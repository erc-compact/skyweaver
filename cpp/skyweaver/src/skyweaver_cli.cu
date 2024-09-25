#include "boost/program_options.hpp"
#include "errno.h"
#include "psrdada_cpp/cli_utils.hpp"
#include "skyweaver/BeamformerPipeline.cuh"
#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/IncoherentDedispersionPipeline.cuh"
#include "skyweaver/MultiFileReader.cuh"
#include "skyweaver/MultiFileWriter.cuh"
#include "skyweaver/PipelineConfig.hpp"
#include "skyweaver/StatisticsCalculator.cuh"
#include "skyweaver/Timer.hpp"
#include "skyweaver/logging.hpp"
#include "skyweaver/nvtx_utils.h"
#include "skyweaver/skyweaver_constants.hpp"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <memory>
#include <omp.h>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <thread>
#include <vector>

#define BOOST_LOG_DYN_LINK 1

namespace
{

std::string skyweaver_splash = R"(
   ____ __                                           
  / __// /__ __ __ _    __ ___  ___ _ _  __ ___  ____
 _\ \ /  '_// // /| |/|/ // -_)/ _ `/| |/ // -_)/ __/
/___//_/\_\ \_, / |__,__/ \__/ \_,_/ |___/ \__//_/   
           /___/                                     

)";

const size_t ERROR_IN_COMMAND_LINE     = 1;
const size_t SUCCESS                   = 0;
const size_t ERROR_UNHANDLED_EXCEPTION = 2;

class NullHandler
{
  public:
    template <typename... Args>
    void init(Args... args){};

    template <typename... Args>
    bool operator()(Args... args)
    {
        return false;
    };
};

const char* build_time = __DATE__ " " __TIME__;

void display_constants()
{
    const int boxWidth = 53;
    std::string border(boxWidth, '-');

    auto print_str = [&](const std::string& label, std::string const& value) {
        std::cout << "[ " << std::setw(boxWidth - 4 - value.size()) << std::left
                  << label << value << " ]" << std::endl;
    };

    auto print_int = [&](const std::string& label, int value) {
        std::string value_str = std::to_string(value);
        print_str(label, value_str);
    };

    std::cout << border << std::endl;
    print_str("Build date: ", build_time);
    print_int("SKYWEAVER_NANTENNAS: ", SKYWEAVER_NANTENNAS);
    print_int("SKYWEAVER_NCHANS: ", SKYWEAVER_NCHANS);
    print_int("SKYWEAVER_NBEAMS: ", SKYWEAVER_NBEAMS);
    print_int("SKYWEAVER_CB_TSCRUNCH: ", SKYWEAVER_CB_TSCRUNCH);
    print_int("SKYWEAVER_CB_FSCRUNCH: ", SKYWEAVER_CB_FSCRUNCH);
    print_int("SKYWEAVER_IB_TSCRUNCH: ", SKYWEAVER_IB_TSCRUNCH);
    print_int("SKYWEAVER_IB_FSCRUNCH: ", SKYWEAVER_IB_FSCRUNCH);
    print_int("SKYWEAVER_IB_SUBTRACTION: ", SKYWEAVER_IB_SUBTRACTION);
    std::cout << border << std::endl;
    std::cout << std::endl;
}

} // namespace

// This patching of the << operator is required to allow
// for float vector arguments to boost program options
namespace std
{
std::ostream& operator<<(std::ostream& os, const std::vector<float>& vec)
{
    for(auto item: vec) { os << item << " "; }
    return os;
}
} // namespace std

template <class Pipeline>
void run_pipeline(Pipeline& pipeline,
                  skyweaver::PipelineConfig& config,
                  skyweaver::MultiFileReader& file_reader,
                  skyweaver::ObservationHeader const& header)
{
    using VoltageType = typename Pipeline::VoltageVectorTypeH;

    BOOST_LOG_NAMED_SCOPE("run_pipeline");
    BOOST_LOG_TRIVIAL(debug) << "Executing pipeline";
    std::size_t input_elements = header.nantennas * config.nchans() *
                                 config.npol() * config.gulp_length_samps();
    BOOST_LOG_TRIVIAL(debug)
        << "Allocating "
        << input_elements * sizeof(typename VoltageType::value_type)
        << " byte input buffer";
    double tsamp = header.obs_nchans / header.obs_bandwidth;
    NVTX_RANGE_PUSH("Input buffer initialisation");
    std::unique_ptr<VoltageType> taftp_input_voltage_a =
        std::make_unique<VoltageType>();
    taftp_input_voltage_a->resize(
        {config.gulp_length_samps() / config.nsamples_per_heap(), // T
         header.nantennas,                                        // A
         config.nchans(),                                         // F
         config.nsamples_per_heap(),                              // T
         config.npol()});                                         // P
    taftp_input_voltage_a->frequencies(config.channel_frequencies());
    taftp_input_voltage_a->tsamp(tsamp);
    taftp_input_voltage_a->dms({0.0});
    std::unique_ptr<VoltageType> taftp_input_voltage_b =
        std::make_unique<VoltageType>();
    taftp_input_voltage_b->like(*taftp_input_voltage_a);
    NVTX_RANGE_POP();

    BOOST_LOG_TRIVIAL(debug)
        << "Input buffer: " << taftp_input_voltage_a->describe();
    std::size_t input_bytes = taftp_input_voltage_a->size() *
                              sizeof(typename VoltageType::value_type);

    NVTX_RANGE_PUSH("Getting total file size");
    std::size_t total_bytes = file_reader.get_total_size();
    BOOST_LOG_TRIVIAL(info) << "Total input size (bytes): " << total_bytes;
    NVTX_RANGE_POP();

    // Set the start offsets and adjust the total bytes
    std::size_t bytes_per_sample =
        header.nantennas * config.nchans() * config.npol() * sizeof(char2);
    std::size_t bytes_per_second = (1.0f / tsamp) * bytes_per_sample;
    std::size_t offset_nsamps =
        static_cast<std::size_t>(config.start_time() / tsamp);
    offset_nsamps = (offset_nsamps / config.nsamples_per_heap()) *
                    config.nsamples_per_heap(); // floor
    pipeline.init(header, offset_nsamps * tsamp);
    std::size_t offset_nbytes = offset_nsamps * bytes_per_sample;
    BOOST_LOG_TRIVIAL(info) << "Starting at " << config.start_time()
                            << " seconds into the observation";
    BOOST_LOG_TRIVIAL(debug)
        << "Offsetting to byte " << offset_nbytes << " of the input data";
    file_reader.seekg(offset_nbytes, std::ios::beg);

    float total_duration     = total_bytes / bytes_per_second;
    float remaining_duration = total_duration - config.start_time();
    if((config.duration() < std::numeric_limits<float>::infinity()) &&
       (config.duration() > remaining_duration)) {
        BOOST_LOG_TRIVIAL(warning)
            << "Requested duration is longer than the remaining input length";
    }
    remaining_duration = std::min(remaining_duration, config.duration());

    skyweaver::Timer stopwatch;
    stopwatch.start("processing_loop");
    std::size_t processed_bytes = 0;
    float data_time_elapsed     = 0.0f;
    float wall_time_elapsed     = 0.0f;
    float real_time_fraction    = 0.0;
    float percentage            = 0.0f;
    std::streamsize nbytes_read = 0;

    // Populate buffer A
    stopwatch.start("file read");
    nbytes_read =
        file_reader.read(reinterpret_cast<char*>(thrust::raw_pointer_cast(
                             taftp_input_voltage_a->data())),
                         input_bytes);
    stopwatch.stop("file read");
    bool thread_error = false;
    // A is full B is empty
    while(!file_reader.eof()) {
        taftp_input_voltage_a.swap(taftp_input_voltage_b);
        // B is full A is empty

        // Here spawn a thread to read the next block to process
        // Thread must write to buffer A
        std::thread reader_thread([&]() {
            try {
                nbytes_read = file_reader.read(
                    reinterpret_cast<char*>(thrust::raw_pointer_cast(
                        taftp_input_voltage_a->data())),
                    input_bytes);
                BOOST_LOG_TRIVIAL(debug)
                    << "read " << nbytes_read << " bytes from file";
            } catch(std::runtime_error& e) {
                BOOST_LOG_TRIVIAL(error) << "Error on input read: " << e.what();
                thread_error = true;
            }
        });
        // Buffer B is full from the previous read and so is now ready to be
        // processed
        pipeline(*taftp_input_voltage_b);
        // Buffer B is now finished processing and we can wait on the A read
        reader_thread.join();
        if(thread_error) {
            throw std::runtime_error("Error in input file read");
        }
        data_time_elapsed +=
            config.gulp_length_samps() * taftp_input_voltage_a->tsamp();
        percentage =
            std::min(100.0f * data_time_elapsed / remaining_duration, 100.0f);
        wall_time_elapsed  = stopwatch.elapsed("processing_loop") / 1e6;
        real_time_fraction = data_time_elapsed / wall_time_elapsed;
        processed_bytes += input_bytes;
        BOOST_LOG_TRIVIAL(info)
            << "Progress: " << std::setprecision(6) << percentage
            << "%, Data time: " << data_time_elapsed
            << " s, Wall time: " << wall_time_elapsed << ", "
            << "Realtime fraction: " << real_time_fraction;
        if(data_time_elapsed >= config.duration()) {
            break;
        }
    }
    stopwatch.stop("processing_loop");
    stopwatch.show_all_timings();
}

auto pre_write_callback = [] (std::size_t const size, skyweaver::MultiFileWriterConfig const& config)
{
    if (!config.pre_write.is_enabled) return;
    std::filesystem::space_info space = std::filesystem::space(config.output_dir);
    size_t limit = std::min(config.pre_write.wait.min_free_space, size);
    if(space.available >= limit)
        return;

    BOOST_LOG_TRIVIAL(info)
        << space.available
        << " bytes available space is not enough. Need at least "
        << limit
        << " bytes in " << config.output_dir << ".";
    BOOST_LOG_TRIVIAL(warning) << "Start pausing.";
    int max_iterations = (config.pre_write.wait.iterations == 0) ? INT_MAX : config.pre_write.wait.iterations;
    for (int i = 0; i < max_iterations; i++)
    {
        sleep(config.pre_write.wait.sleep_time);
        space = std::filesystem::space(config.output_dir);
        if (space.available >= limit)
        {
            BOOST_LOG_TRIVIAL(info) << "Space has been freed up. Will proceed.";
            return;
        }
    }
    throw std::runtime_error("Space for writing hasn't been freed up in time.");
};

template <typename BfTraits, bool enable_incoherent_dedispersion>
void setup_pipeline(skyweaver::PipelineConfig& config)
{
    BOOST_LOG_NAMED_SCOPE("setup_pipeline");
    BOOST_LOG_TRIVIAL(debug) << "Setting up the pipeline";
    // Update the config
    NVTX_RANGE_PUSH("File reader initialisation and header fetch");
    skyweaver::MultiFileReader file_reader(config);
    auto const& header = file_reader.get_header();
    NVTX_RANGE_POP();
    BOOST_LOG_TRIVIAL(debug) << "Validating headers and updating configuration";
    validate_header(header, config);
    update_config(config, header);
    BOOST_LOG_TRIVIAL(debug) << "Creating pipeline handlers";
    using OutputType = typename BfTraits::QuantisedPowerType;

    using IBWriterType =
        skyweaver::MultiFileWriter<skyweaver::BTFPowersH<OutputType>>;
    typename IBWriterType::CreateStreamCallBackType create_stream_callback_ib =
        skyweaver::detail::create_dada_file_stream<
            skyweaver::BTFPowersH<OutputType>>;

    std::unique_ptr<IBWriterType> ib_handler;

    using StatsWriterType =
        skyweaver::MultiFileWriter<skyweaver::FPAStatsD<skyweaver::Statistics>>;
    std::unique_ptr<StatsWriterType> stats_handler;

    typename StatsWriterType::CreateStreamCallBackType
        create_stream_callback_stats =
            skyweaver::detail::create_dada_file_stream<
                skyweaver::FPAStatsD<skyweaver::Statistics>>;
    if constexpr(enable_incoherent_dedispersion) {
        using CBWriterType = skyweaver::MultiFileWriter<skyweaver::TDBPowersH<OutputType>>;   
        typename CBWriterType::CreateStreamCallBackType
        create_stream_callback_cb =
            skyweaver::detail::create_dada_file_stream<skyweaver::TDBPowersH<OutputType>>;
        std::unique_ptr<CBWriterType> cb_file_writer;
        if (config.pre_write_config().is_enabled)
        {
            ib_handler.reset(new IBWriterType(config, "ib", create_stream_callback_ib, pre_write_callback));
            stats_handler.reset(new StatsWriterType(config, "stats", create_stream_callback_stats, pre_write_callback));
            cb_file_writer.reset(new CBWriterType(config, "cb", create_stream_callback_cb, pre_write_callback));
        }else{
            ib_handler.reset(new IBWriterType(config, "ib", create_stream_callback_ib));
            stats_handler.reset(new StatsWriterType(config, "stats", create_stream_callback_stats));
            cb_file_writer.reset(new CBWriterType(config, "cb", create_stream_callback_cb));
        }

        skyweaver::IncoherentDedispersionPipeline<OutputType,
                                                  OutputType,
                                                  decltype(* cb_file_writer.get())>
        incoherent_dispersion_pipeline(config, * cb_file_writer.get());
        skyweaver::BeamformerPipeline<decltype(incoherent_dispersion_pipeline),
                                      decltype(* ib_handler.get()),
                                      decltype(* stats_handler.get()),
                                      BfTraits>
            pipeline(config,
                     incoherent_dispersion_pipeline,
                     * ib_handler.get(),
                     * stats_handler.get());
        run_pipeline(pipeline, config, file_reader, header);
    } else {
        using CBWriterType = skyweaver::MultiFileWriter<skyweaver::TFBPowersD<OutputType>>;   
        std::unique_ptr<CBWriterType> cb_file_writer;
        typename CBWriterType::CreateStreamCallBackType
        create_stream_callback_cb =
            skyweaver::detail::create_dada_file_stream<skyweaver::TFBPowersD<OutputType>>;           
        if (config.pre_write_config().is_enabled)
        {
            ib_handler.reset(new IBWriterType(config, "ib", create_stream_callback_ib, pre_write_callback));
            cb_file_writer.reset(new CBWriterType(config, "cb", create_stream_callback_cb, pre_write_callback));
            stats_handler.reset(new StatsWriterType(config, "stats", create_stream_callback_stats, pre_write_callback));
        }else{
            ib_handler.reset(new IBWriterType(config, "ib", create_stream_callback_ib));
            cb_file_writer.reset(new CBWriterType(config, "cb", create_stream_callback_cb));
            stats_handler.reset(new StatsWriterType(config, "stats", create_stream_callback_stats));
        }

        skyweaver::BeamformerPipeline<decltype(*cb_file_writer.get()),
                                      decltype(*ib_handler.get()),
                                      decltype(* stats_handler.get()),
                                      BfTraits>
            pipeline(config, * cb_file_writer.get(), * ib_handler.get(), * stats_handler.get());
        run_pipeline(pipeline, config, file_reader, header);
    }
}

int main(int argc, char** argv)
{
    std::cout << skyweaver_splash;
    display_constants();

    try {
        skyweaver::PipelineConfig config;
        skyweaver::init_logging("warning");
        /**
         * Define and parse the program options
         */
        namespace po = boost::program_options;

        // Generic options group here to contain the configuration file name
        // The config file will be parsed AFTER the the command line options
        // this gives command line options precedence over config file options.
        // Options here are only settable via the command line.
        po::options_description generic("Generic options");
        generic.add_options()("cfg,c",
                              po::value<std::string>()->default_value(""),
                              "Skyweaver configuration file");

        // Main option group that contains parameters settable via both the
        // command line and the config file
        po::options_description main_options("Main options");
        main_options.add_options()

            // Help menu
            ("help,h", "Display help messages")

            // Input file containing list of DADA files to process
            ("input-file",
             po::value<std::string>()->required()->notifier(
                 [&config](std::string key) {
                     config.read_input_file_list(key);
                 }),
             "File containing list of DADA files to process")

            // Input file for delay solutions
            // This can contain any number of beams but only beams
            // up to the maximum configured for generation will be
            // produced. Antenna ordering in the file must match
            // the antenna order of the input data.
            ("delay-file",
             po::value<std::string>()->required()->notifier(
                 [&config](std::string key) { config.delay_file(key); }),
             "File containing delay solutions")

            // Output directory where all results will be written
            ("output-dir",
             po::value<std::string>()
                 ->default_value(config.output_dir())
                 ->notifier(
                     [&config](std::string key) { config.output_dir(key); }),
             "The output directory for all results")

            // Output file for block statistics
            ("output-level",
             po::value<float>()
                 ->default_value(config.output_level())
                 ->notifier([&config](float key) { config.output_level(key); }),
             "The desired standard deviation for output data")

            /**
             * Defines a dedispersion plan to be executed.
             * Argument is colon-separated with no spaces.
             * Parameters:
             * <coherent_dm>:<start_incoherent_dm>:<end_incoherent_dm>:<dm_step>:<tscrunch>
             * The tscrunch is defined relative to the beamformer output.
             *
             * --ddplan 5.0:0.0:10.0:0.1:1
             * or
             * --ddplan 5.0:1
             *
             * In the configuration file it can only be specified with:
             *
             * ddplan=5.0:0.0:10.0:0.1:1
             * ddplan=5.0:1
             */
            ("ddplan",
             po::value<std::vector<std::string>>()
                 ->multitoken()
                 ->required()
                 ->notifier(
                     [&config](std::vector<std::string> const& descriptors) {
                         for(auto const& descriptor: descriptors) {
                             config.ddplan().add_block(descriptor);
                         }
                     }),
             "A dispersion plan definition string "
             "(<coherent_dm>:<start_incoherent_dm>:"
             "<end_incoherent_dm>:<dm_step>:<tscrunch>) or "
             "(<coherent_dm>:<tscrunch>) "
             "or (<coherent_dm>)")

                ("enable-incoherent-dedispersion",
                 po::value<bool>()->default_value(true)->notifier(
                     [&config](bool const& enable) {
                         config.enable_incoherent_dedispersion(enable);
                     }),
                 "Turn on/off incoherent dedispersion after beamforming")

                    ("start-time",
                     po::value<float>()->default_value(0.0f)->notifier(
                         [&config](float const& start_time) {
                             config.start_time(start_time);
                         }),
                     "Time since start of the data stream from which to start "
                     "processing (seconds)")

                        ("duration",
                         po::value<float>()
                             ->default_value(
                                 std::numeric_limits<float>::infinity())
                             ->notifier([&config](float const& duration) {
                                 config.duration(duration);
                             }),
                         "Number of seconds of data to process")

            // Number of samples to read in each gulp
            ("gulp-size",
             po::value<std::size_t>()
                 ->default_value(config.gulp_length_samps())
                 ->notifier([&config](std::size_t const& gulp_size) {
                     // Round off to next multiple of 256
                     if(gulp_size % config.nsamples_per_heap() != 0 || // 256
                        gulp_size % config.nsamples_per_block() !=
                            0) { // 32* tscrunch

                         std::size_t larger_value =
                             std::max(config.nsamples_per_heap(),
                                      config.nsamples_per_block());

                         BOOST_LOG_TRIVIAL(warning)
                             << "Rounding up gulp-size to next multiple of "
                                "NSAMPLES PER HEAP and NSAMPLES PER BLOCK";

                         config.gulp_length_samps((gulp_size / larger_value) *
                                                  larger_value);
                     } else {
                         config.gulp_length_samps(gulp_size);
                     }
                 }),
             "The number of samples to read in each gulp ")

            // Stokes mode I, Q, U, V or IQUV
            ("stokes-mode",
             po::value<std::string>()
                 ->default_value(config.stokes_mode())
                 ->notifier([&config](std::string stokes) {
                     for(auto& c: stokes) c = (char)toupper(c);
                     config.stokes_mode(stokes);
                 }),
             "The Stokes mode to use, can be either I, Q, U, V or IQUV")

            // Logging options
            ("nthreads",
             po::value<std::size_t>()->default_value(16)->notifier(
                 [](std::size_t nthreads) { omp_set_num_threads(nthreads); }),
             "The number of threads to use for incoherent dedispersion")

            ("wait-for-space",
              po::value<std::string>()
                 ->notifier(
                     [&config](std::string key) { config.configure_wait(key); }),
             "Wait for enough disk space for the output. "
             "<nr_of_cycles (0 for infinity)>:<sleep_time_per_cycle [s]>:<minimum_size(e.g. 200M or 32G)>")

          // Logging options
            ("log-level",
             po::value<std::string>()->default_value("info")->notifier(
                 [](std::string level) { skyweaver::set_log_level(level); }),
             "The logging level to use (debug, info, warning, error)");

        // set options allowed on command line
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(main_options);

        // set options allowed in config file
        po::options_description config_file_options;
        config_file_options.add(main_options);

        po::variables_map variable_map;
        try {
            po::store(po::command_line_parser(argc, argv)
                          .options(cmdline_options)
                          .run(),
                      variable_map);
            if(variable_map.count("help")) {
                std::cout << "skyweavercpp -- C++/CUDA beamformer pipeline for "
                             "COMPACT-ERC"
                          << std::endl
                          << cmdline_options << std::endl;
                return SUCCESS;
            }
        } catch(po::error& e) {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

        auto config_file = variable_map.at("cfg").as<std::string>();

        if(config_file != "") {
            std::ifstream config_fs(config_file.c_str());
            if(!config_fs.is_open()) {
                std::cerr << "Unable to open configuration file: "
                          << config_file << " (" << std::strerror(errno)
                          << ")\n";
                return ERROR_UNHANDLED_EXCEPTION;
            } else {
                po::store(po::parse_config_file(config_fs, config_file_options),
                          variable_map);
            }
        }
        po::notify(variable_map);

        /**
         * All the application code goes here
         */
        NVTX_MARKER("Application start");
        BOOST_LOG_NAMED_SCOPE("skyweaver_cli")
        BOOST_LOG_TRIVIAL(info)
            << "Initialising the skyweaver beamforming pipeline";
        if(config_file != "") {
            BOOST_LOG_TRIVIAL(info) << "Configuration file: " << config_file;
        }
        BOOST_LOG_TRIVIAL(info)
            << "Input file count: " << config.input_files().size();
        BOOST_LOG_TRIVIAL(info) << "Delay file: " << config.delay_file();
        BOOST_LOG_TRIVIAL(info) << "Output dir: " << config.output_dir();
        BOOST_LOG_TRIVIAL(info) << "Output level: " << config.output_level();
        BOOST_LOG_TRIVIAL(info) << "Gulp size: " << config.gulp_length_samps();
        BOOST_LOG_TRIVIAL(info) << "Stokes mode: " << config.stokes_mode();
        BOOST_LOG_TRIVIAL(info) << config.ddplan();
        if(config.enable_incoherent_dedispersion()) {
            BOOST_LOG_TRIVIAL(info) << "Incoherent dedispersion enabled";
            if(config.stokes_mode() == "I") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::I>,
                               true>(config);
            } else if(config.stokes_mode() == "Q") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::Q>,
                               true>(config);
            } else if(config.stokes_mode() == "U") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::U>,
                               true>(config);
            } else if(config.stokes_mode() == "V") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::V>,
                               true>(config);
            } else if(config.stokes_mode() == "QU") {
                setup_pipeline<skyweaver::StokesTraits<
                                   skyweaver::StokesParameter::Q, skyweaver::StokesParameter::U>,
                               true>(config);
            } else if(config.stokes_mode() == "IV") {
                setup_pipeline<skyweaver::StokesTraits<
                                   skyweaver::StokesParameter::I, skyweaver::StokesParameter::V>,
                               true>(config);             
            } else if(config.stokes_mode() == "IQUV") {
                setup_pipeline<skyweaver::FullStokesBeamformerTraits, true>(
                    config);
            } else {
                throw std::runtime_error(
                    "Invalid Stokes mode passed, must be one "
                    "of I, Q, U, V, QU, IV or IQUV");
            }
        } else {
            BOOST_LOG_TRIVIAL(info) << "Incoherent dedispersion disabled";
            if(config.stokes_mode() == "I") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::I>,
                               false>(config);
            } else if(config.stokes_mode() == "Q") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::Q>,
                               false>(config);
            } else if(config.stokes_mode() == "U") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::U>,
                               false>(config);
            } else if(config.stokes_mode() == "V") {
                setup_pipeline<skyweaver::SingleStokesBeamformerTraits<
                                   skyweaver::StokesParameter::V>,
                               false>(config);
            } else if(config.stokes_mode() == "QU") {
                setup_pipeline<skyweaver::StokesTraits<
                                   skyweaver::StokesParameter::Q, skyweaver::StokesParameter::U>,
                               false>(config);
            } else if(config.stokes_mode() == "IV") {
                setup_pipeline<skyweaver::StokesTraits<
                                   skyweaver::StokesParameter::I, skyweaver::StokesParameter::V>,
                               false>(config); 
            } else if(config.stokes_mode() == "IQUV") {
                setup_pipeline<skyweaver::FullStokesBeamformerTraits, false>(
                    config);
            } else {
                throw std::runtime_error(
                    "Invalid Stokes mode passed, must be one "
                    "of I, Q, U, V, QU, IV or IQUV");
            }
        }
    } catch(std::exception& e) {
        std::cerr << "Unhandled Exception reached the top of main: " << e.what()
                  << ", application will now exit" << std::endl;
        return ERROR_UNHANDLED_EXCEPTION;
    }
    return SUCCESS;
}
