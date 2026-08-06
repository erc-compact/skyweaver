#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/SkyCleaver.hpp"
#include "skyweaver/SkyCleaverConfig.hpp"
#include "skyweaver/logging.hpp"
#include "skyweaver/skycleaver_utils.hpp"

#include <algorithm>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <cerrno>
#include <cmath>
#include <cstdlib>
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
#include <type_traits>
#include <vector>

#define BOOST_LOG_DYN_LINK 1

namespace
{

std::string skycleaver_splash          = R"(
        __                __                               
.-----.|  |--.--.--.----.|  |.-----.---.-.--.--.-----.----.
|__ --||    <|  |  |  __||  ||  -__|  _  |  |  |  -__|   _|
|_____||__|__|___  |____||__||_____|___._|\___/|_____|__|  
             |_____|                                       

)";
const size_t ERROR_IN_COMMAND_LINE     = 1;
const size_t SUCCESS                   = 0;
const size_t ERROR_UNHANDLED_EXCEPTION = 2;

const char* build_time = __DATE__ " " __TIME__;



} // namespace

namespace std
{
std::ostream& operator<<(std::ostream& os, const std::vector<float>& vec)
{
    for(auto item: vec) { os << item << " "; }
    return os;
}
} // namespace std
template <typename InputVectorType, typename OutputVectorType>
void run_pipeline(skyweaver::SkyCleaverConfig& config)
{
    skyweaver::SkyCleaver<InputVectorType, OutputVectorType> skycleaver(config);
    skycleaver.cleave();
}

int main(int argc, char** argv)
{
    std::cout << skycleaver_splash;
    std::cout << "Build time: " << build_time << std::endl;
    // skyweaver::init_logging("warning");

    skyweaver::SkyCleaverConfig config;

    namespace po = boost::program_options;

    po::options_description generic("Generic options");
    generic.add_options()("cfg,c",
                          po::value<std::string>()->default_value(""),
                          "Skycleaver configuration file");

    po::options_description main_options("Main options");
    main_options.add_options()("help,h", "Produce help message")(
        "root-dir,r",
        po::value<std::string>()->required()->notifier(
            [&config](std::string key) { config.root_dir(key); }),
        "The output directory for all results")(
        "output-dir",
        po::value<std::string>()
            ->default_value(config.output_dir())
            ->notifier([&config](std::string key) { config.output_dir(key); }),
        "The output directory for all results")(
        "out-prefix",
        po::value<std::string>()
            ->default_value(config.out_prefix())
            ->notifier([&config](std::string key) { config.out_prefix(key); }),
        "The prefix for all output files")(
        "nthreads",
        po::value<unsigned int>()
            ->default_value(config.nthreads())
            ->notifier([&config](unsigned int key) { config.nthreads(key); }),
        "The number of threads to use for processing")(
        "nsamples-per-block",
        po::value<std::size_t>()
            ->default_value(config.nsamples_per_block())
            ->notifier(
                [&config](std::size_t key) { config.nsamples_per_block(key); }),
        "The number of samples per block")(
        "nbridges",
        po::value<std::size_t>()
            ->default_value(config.nbridges())
            ->notifier([&config](std::size_t key) { config.nbridges(key); }),
        "The number of bridges")(
        "stream-id",
        po::value<std::size_t>()
            ->default_value(config.stream_id())
            ->notifier([&config](std::size_t key) { config.stream_id(key); }),
        "The stream id")(
        "max-ram-gb",
        po::value<std::size_t>()
            ->default_value(config.max_ram_gb())
            ->notifier([&config](std::size_t key) { config.max_ram_gb(key); }),
        "The maximum amount of RAM to use in GB")(
        "max-output-filesize",
        po::value<std::size_t>()
            ->default_value(config.max_output_filesize())
            ->notifier([&config](std::size_t key) {
                config.max_output_filesize(key);
            }),
        "The maximum output file size in bytes")(
        "dada-header-size",
        po::value<std::size_t>()
            ->default_value(config.dada_header_size())
            ->notifier(
                [&config](std::size_t key) { config.dada_header_size(key); }),
        "The size of the DADA header")(
        "log-level",
        po::value<std::string>()->default_value("info")->notifier(
            [](std::string level) { skyweaver::init_logging(level); }),
        "The logging level to use (debug, info, warning, error)")(
        "start_sample",
        po::value<std::size_t>()
            ->default_value(config.start_sample())
            ->notifier(
                [&config](std::size_t key) { config.start_sample(key); }),
        "Start from this sample")(
        "nsamples_to_read",
        po::value<std::size_t>()
            ->default_value(config.start_sample())
            ->notifier(
                [&config](std::size_t key) { config.nsamples_to_read(key); }),
        "total number of samples to read from start_sample")(
            "targets_file",
            po::value<std::string>()->default_value("")->notifier(
                [&config](std::string key) {
                    config.targets_file(key);
                }), "update beam names and positions from this file"
        )(
        "required_beams",
        po::value<std::string>()->default_value("")->notifier(
            [&config](std::string key) {
                config.required_beams(skyweaver::get_list_from_string<int>(key));
            }),
        "Comma separated list of beams to process. Syntax - beam1, beam2, "
        "beam3:beam4:step, beam5 etc..")(
        "required_dms",
        po::value<std::string>()->default_value("")->notifier(
            [&config](std::string key) {
                config.required_dms(skyweaver::get_list_from_string<double>(key));
            }),
        "Comma separated list of DMs to process. Syntax - dm1, dm2, "
        "dm1:dm2:step, etc..")(
        "out-stokes",
        po::value<std::string>()
            ->default_value(config.out_stokes())
            ->notifier([&config](std::string key) {
                if(key.find_first_not_of("IQUVL") != std::string::npos) {
                    throw std::runtime_error("Invalid Stokes mode: " + key);
                }
                config.out_stokes(key);
            }),
        "The list of stokes needed - these will be separte output files");

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(main_options);

    // set options allowed in config file
    po::options_description config_file_options;
    config_file_options.add(main_options);
    po::variables_map variable_map;
    try {
        po::store(
            po::command_line_parser(argc, argv).options(cmdline_options).run(),
            variable_map);
        if(variable_map.count("help")) {
            std::cout
                << "skycleaver -- A pipeline that cleaves input TDB files, "
                   "and cleaves them to form output Sigproc Filterbank files."
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
            std::cerr << "Unable to open configuration file: " << config_file
                      << " (" << std::strerror(errno) << ")\n";
            return ERROR_UNHANDLED_EXCEPTION;
        } else {
            po::store(po::parse_config_file(config_fs, config_file_options),
                      variable_map);
        }
    }
    po::notify(variable_map);

    BOOST_LOG_NAMED_SCOPE("skycleaver_cli");
    BOOST_LOG_TRIVIAL(info) << "Configuration: " << config_file;
    BOOST_LOG_TRIVIAL(info) << "root_dir: " << config.root_dir();
    BOOST_LOG_TRIVIAL(info) << "output_dir: " << config.output_dir();
    BOOST_LOG_TRIVIAL(info) << "out_prefix: " << config.out_prefix();
    BOOST_LOG_TRIVIAL(info) << "nthreads: " << config.nthreads();
    BOOST_LOG_TRIVIAL(info)
        << "nsamples_per_block: " << config.nsamples_per_block();
    BOOST_LOG_TRIVIAL(info) << "nchans: " << config.nchans();
    BOOST_LOG_TRIVIAL(info) << "nbeams: " << config.nbeams();
    BOOST_LOG_TRIVIAL(info) << "nbridges: " << config.nbridges();
    BOOST_LOG_TRIVIAL(info) << "ndms: " << config.ndms();
    BOOST_LOG_TRIVIAL(info) << "max_ram_gb: " << config.max_ram_gb();
    BOOST_LOG_TRIVIAL(info)
        << "max_output_filesize: " << config.max_output_filesize();
    BOOST_LOG_TRIVIAL(info)
        << "dada_header_size: " << config.dada_header_size();
    BOOST_LOG_TRIVIAL(info) << "start_sample: " << config.start_sample();
    BOOST_LOG_TRIVIAL(info)
        << "nsamples_to_read: " << config.nsamples_to_read();
    BOOST_LOG_TRIVIAL(info) << "out_stokes: " << config.out_stokes();

    if(config.required_beams().size() > 0) {
        for(auto beam: config.required_beams()) {
            BOOST_LOG_TRIVIAL(info) << "required_beam: " << beam;
        }
    }
    if(config.required_dms().size() > 0) {
        for(auto dm: config.required_dms()) {
            BOOST_LOG_TRIVIAL(info) << "required_dm: " << dm;
        }
    }
    

    run_pipeline<skyweaver::TDBPowersStdH<char>,
                     skyweaver::TFPowersStdH<uint8_t>>(config);


    // skyweaver::SkyCleaver skycleaver(config);
    // skycleaver.cleave();

    // if(config.stokes_mode() == "I" || config.stokes_mode() == "Q" ||
    //    config.stokes_mode() == "U" || config.stokes_mode() == "V") {
    //     run_pipeline<skyweaver::TDBPowersStdH<std::vector<char>>,
    //                  skyweaver::TFPowersStdH<uint8_t>>(config);
    // } else if(config.stokes_mode() == "IV" || config.stokes_mode() == "QU") {
    //     run_pipeline<skyweaver::TDBPowersStdH<char2>,
    //                  skyweaver::TFPowersStdH<uint8_t>>(config);

    // } else if(config.stokes_mode() == "IQUV") {
    //     run_pipeline<skyweaver::TDBPowersStdH<char4>,
    //                  skyweaver::TFPowersStdH<uint8_t>>(config);
    // } else {
    //     throw std::runtime_error("Invalid Stokes mode: " +
    //                              config.stokes_mode());
    // }

}
