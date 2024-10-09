
#include <cmath>
#include <filesystem>
#include <iostream>
#include <omp.h>
#include <regex>
#include <ranges>
#include <type_traits>
#include "skyweaver/types.cuh"
#include "skyweaver/SkyCleaver.cuh"

namespace fs = std::filesystem;
using SkyCleaver       = skyweaver::SkyCleaver;


using FreqType         = skyweaver::SkyCleaver::FreqType;
using BridgeReader     = skyweaver::BridgeReader;
using MultiFileReader  = skyweaver::MultiFileReader;
using OutputVectorType = skyweaver::SkyCleaver::OutputVectorType;
using InputVectorType  = skyweaver::SkyCleaver::InputVectorType;

namespace
{
template <typename T>
std::string to_string_with_padding(T num, int width, int precision = -1)
{
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0');
    if(precision >=
       0) { // Check if precision is specified for floating-point numbers
        oss << std::fixed << std::setprecision(precision);
    }
    oss << num;
    return oss.str();
}
std::vector<std::string>
get_subdirs(std::string directory_path,
            std::regex numeric_regex = std::regex("^[0-9]+$"))
{
    std::vector<std::string> subdirs;
    try {
        if(fs::exists(directory_path) && fs::is_directory(directory_path)) {
            for(const auto& entry: fs::directory_iterator(directory_path)) {
                if(fs::is_directory(entry.status())) {
                    std::string folder_name = entry.path().filename().string();
                    if(std::regex_match(folder_name, numeric_regex)) {
                        BOOST_LOG_TRIVIAL(debug)
                            << "Found subdirectory: " << folder_name;
                        subdirs.push_back(folder_name);
                    }
                }
            }
        } else {
            std::runtime_error(
                "Root directory does not exist or is not a directory.");
        }
    } catch(const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        std::runtime_error("Error reading subdirectories in root directory: " +
                           directory_path);
    }

    return subdirs;
}

std::vector<std::string> get_files(std::string directory_path,
                                   std::string extension)
{
    std::vector<std::string> files;
    try {
        if(fs::exists(directory_path) && fs::is_directory(directory_path)) {
            for(const auto& entry: fs::directory_iterator(directory_path)) {
                if(fs::is_regular_file(entry.status())) {
                    std::string file_name = entry.path().string();
                    if(file_name.find(extension) != std::string::npos) {
                        files.push_back(file_name);
                    }
                }
            }
        } else {
            std::runtime_error("No files in bridge directory: " +
                               directory_path);
        }
    } catch(const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        std::runtime_error("Error reading files in bridge directory: " +
                           directory_path);
    }

    return files;
}

} // namespace

void SkyCleaver::init_readers()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::init_readers")

    std::string root_dir    = _config.root_dir();
    std::string root_prefix = _config.root_prefix();
    std::size_t stream_id   = _config.stream_id();

    // get the list of directories in root/stream_id(for the nex)
    std::vector<std::string> freq_dirs =
        get_subdirs(root_dir + "/" + std::to_string(stream_id));

    BOOST_LOG_TRIVIAL(info)
        << "Found " << freq_dirs.size()
        << " frequency directories in root directory: " << root_dir;

    std::map<FreqType, long double> bridge_timestamps;
    long double latest_timestamp = 0.0;

    for(const auto& freq_dir: freq_dirs) {
        std::vector<std::string> tdb_files = get_files(
            root_dir + "/" + std::to_string(stream_id) + "/" + freq_dir,
            ".tdb");
        BOOST_LOG_TRIVIAL(info) << "Found " << tdb_files.size()
                                << " TDB files for frequency: " << freq_dir;
        if(tdb_files.empty()) {
            BOOST_LOG_TRIVIAL(warning)
                << "No TDB files found for frequency: " << freq_dir;
            continue;
        }

        std::size_t freq = static_cast<std::size_t>(std::stoul(freq_dir));

        _bridge_readers[freq] =
            std::make_unique<MultiFileReader>(tdb_files,
                                              _config.dada_header_size(),
                                              false);
        long double timestamp = _bridge_readers[freq]->get_header().utc_start;
        bridge_timestamps.insert({freq, timestamp});
        if(timestamp > latest_timestamp) {
            latest_timestamp = timestamp;
        }
        _available_freqs.push_back(freq);

        BOOST_LOG_TRIVIAL(debug)
            << "Added bridge reader for frequency: " << freq_dir;
    }

    int nbridges = _config.nbridges();

    std::size_t gulp_size =
        _config.nsamples_per_block() * _config.ndms() * _config.nbeams();

    ObservationHeader const& header =
        (*_bridge_readers.begin()).second->get_header();
    BOOST_LOG_TRIVIAL(info)
        << "Read header from first file: " << header.to_string();

    long double obs_centre_freq = header.obs_frequency;
    long double obs_bandwidth = header.obs_bandwidth;

    long double start_freq = obs_centre_freq - obs_bandwidth / 2;



    for(int i = 0; i < nbridges; i++) {
       
        int ifreq = std::lround(std::floor(start_freq + (i + 0.5) * obs_bandwidth / nbridges));
        _expected_freqs.push_back(ifreq);
        BOOST_LOG_TRIVIAL(info)
            << "Expected frequency [" << i << "]: " << ifreq;

        if(_bridge_readers.find(ifreq) == _bridge_readers.end()) {
            BOOST_LOG_TRIVIAL(warning)
                << "Frequency " << ifreq
                << " not found in bridge readers, will write zeros";
        }
        _bridge_data[ifreq] = std::make_unique<InputVectorType>(
            std::initializer_list{_config.nsamples_per_block(),
                                  _config.ndms(),
                                  _config.nbeams()},
            0);
    }

    // print _expected_freqs

    std::size_t smallest_data_size = std::numeric_limits<std::size_t>::max();

    for(const auto& [freq, reader]: _bridge_readers) {
        // at this point, all non-existed frequencies have been added with zero
        // data now check if there are any unexpected frequencies in the bridge
        // readers.
        if(std::find(_expected_freqs.begin(), _expected_freqs.end(), freq) ==
           _expected_freqs.end()) {
            throw std::runtime_error("Frequency " + std::to_string(freq) +
                                     " not found in expected frequencies");
        }

        // now time align all the bridges to the latest timestamp
        long double timestamp = bridge_timestamps[freq];
        long double time_diff = latest_timestamp - timestamp;
        long double tsamp =
            reader->get_header().tsamp *
            1e-6; // Header has it in microseconds, converting to seconds
        std::size_t nsamples = std::floor(time_diff / tsamp);

        BOOST_LOG_TRIVIAL(debug)
            << "Frequency: " << freq << " Timestamp: " << timestamp
            << "tsamp: " << tsamp << " Latest timestamp: " << latest_timestamp
            << " Time difference: " << time_diff
            << " Number of samples to skip: " << nsamples;

        BOOST_LOG_TRIVIAL(debug)
            << "Seeking " << nsamples * _config.ndms() * _config.nbeams()
            << " bytes in bridge reader for frequency: " << freq;

        std::size_t bytes_seeking = (nsamples * _config.ndms() *
                                    _config.nbeams() *
                                    sizeof(InputVectorType::value_type));

        _bridge_readers[freq]->seekg(bytes_seeking,
                                     std::ios_base::beg);

        std::size_t data_size =
            _bridge_readers[freq]->get_total_size() - bytes_seeking;
        BOOST_LOG_TRIVIAL(debug)
            << "Data size for frequency: " << freq << " is " << data_size;
        if(data_size < smallest_data_size) {
            smallest_data_size = data_size;
        }
    }


    if(smallest_data_size % (_config.ndms() * _config.nbeams()) != 0) {
        std::runtime_error("Data size is not a multiple of ndms * nbeams");
    }

    std::size_t smallest_nsamples =
        std::floor(smallest_data_size / _config.ndms() / _config.nbeams());

    if (smallest_nsamples < _config.start_sample()) {
        std::runtime_error("start_sample is greater than the smallest_nsamples in the data.");
    }

    smallest_nsamples = smallest_nsamples - _config.start_sample();

    BOOST_LOG_TRIVIAL(info)
        << "Smallest data size: " << smallest_data_size
        << " Smallest number of samples: " << smallest_nsamples;

    if(smallest_nsamples < _config.nsamples_per_block()) {
        std::runtime_error(
            "Smallest data size is less than nsamples_per_block");
    }

    if(_config.nsamples_to_read() > 0) {

        if(smallest_nsamples < _config.nsamples_to_read()) {
            std::runtime_error(
                "Smallest data size is less than nsamples_to_read");
        }
        

        _nsamples_to_read = _config.nsamples_to_read();
    }
    else{
        _nsamples_to_read = smallest_nsamples;
    }

    std::size_t bytes_seeking =
        (_config.start_sample() * _config.ndms() * _config.nbeams() *
         sizeof(InputVectorType::value_type));

    if(bytes_seeking > 0) {
        BOOST_LOG_TRIVIAL(info) << "Seeking " << bytes_seeking
                                << " bytes in bridge readers to start sample: "
                                << _config.start_sample();
        for(const auto& [freq, reader]: _bridge_readers) {
            _bridge_readers[freq]->seekg(bytes_seeking, std::ios_base::cur);
        }
    }

    BOOST_LOG_TRIVIAL(info)
        << "Added " << _bridge_data.size() << " bridge readers to SkyCleaver";

    _header = _bridge_readers[_available_freqs[0]]->get_header();
    BOOST_LOG_TRIVIAL(info) << "Adding first header to SkyCleaver: " << _header.to_string();
    _header.nchans = _header.nchans * _config.nbridges();
    _header.nbeams = _config.nbeams();
}

void SkyCleaver::init_writers()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::init_writers")
    BOOST_LOG_TRIVIAL(debug)
        << "_config.output_dir(); " << _config.output_dir();


    if(!fs::exists(_config.output_dir())) {
        fs::create_directories(_config.output_dir());
    }
    std::string out_prefix = _config.out_prefix().empty()
                                 ? ""
                                 : _config.out_prefix() + "_";
    std::string output_dir = _config.output_dir();

    

    for(int idm = 0; idm < _config.ndms(); idm++) {

        std::string prefix = _config.ndms() > 1 ? out_prefix + "idm_" +
                         to_string_with_padding(_header.dms[idm], 9, 3) + "_": out_prefix;

        if(! _config.required_dms().empty()) {
            const auto& required_dms = _config.required_dms();
            if(std::ranges::find(required_dms, _header.dms[idm]) == required_dms.end()) {
                BOOST_LOG_TRIVIAL(info)
                    << "DM " << _header.dms[idm] << " is not required, skipping from writing";
                continue;
            }
        }
        for(int ibeam = 0; ibeam < _config.nbeams(); ibeam++) {

            // skip if beam is not used
            if (! _config.required_beams().empty()) {
                const auto& required_beams = _config.required_beams();
                std::cerr << "required_beams: " << required_beams.size() << "ibeam: " << ibeam << std::endl;
                if(std::ranges::find(required_beams, ibeam) == required_beams.end()) {
                    BOOST_LOG_TRIVIAL(info)
                        << "Beam " << ibeam << " is not required, skipping from writing";
                    continue;
                }
            }
            MultiFileWriterConfig writer_config;
            writer_config.header_size   = _config.dada_header_size();
            writer_config.max_file_size = _config.max_output_filesize();
            writer_config.stokes_mode   = _config.stokes_mode();
            writer_config.base_output_dir  = output_dir;
            writer_config.prefix        = prefix + "cb_" + to_string_with_padding(ibeam, 5);;
            writer_config.extension     = ".fil";

            BOOST_LOG_TRIVIAL(info)
                << "Writer config: " << writer_config.to_string();

            typename MultiFileWriter<OutputVectorType>::CreateStreamCallBackType
                create_stream_callback_sigproc =
                    skyweaver::detail::create_sigproc_file_stream<
                        OutputVectorType>;
            _beam_writers[idm][ibeam] =
                std::make_unique<MultiFileWriter<OutputVectorType>>(
                    writer_config,
                    "",
                    create_stream_callback_sigproc);
            _header.ibeam = ibeam;
            _beam_writers[idm][ibeam]->init(_header);

            _beam_data[idm][ibeam] = std::make_shared<OutputVectorType>(
                std::initializer_list{_config.nsamples_per_block(),
                                      _config.nbridges()},
                0);

            _beam_data[idm][ibeam]->reference_dm(_header.refdm);

            _total_beam_writers++;
        }
    }

    BOOST_LOG_TRIVIAL(info)
        << "Added " << _total_beam_writers << " beam writers to SkyCleaver";
}

SkyCleaver::SkyCleaver(SkyCleaverConfig const& config): _config(config)
{
    _timer.start("skycleaver::init_readers");
    init_readers();
    _timer.stop("skycleaver::init_readers");
    _timer.start("skycleaver::init_writers");
    init_writers();
    _timer.stop("skycleaver::init_writers");
}

void SkyCleaver::cleave()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::cleave")

    /** read data **/

    for(std::size_t nsamples_read = 0; nsamples_read < _nsamples_to_read;
        nsamples_read += _config.nsamples_per_block()) {
        std::size_t gulp_samples =
            _nsamples_to_read - nsamples_read < _config.nsamples_per_block()
                ? _nsamples_to_read - nsamples_read
                : _config.nsamples_per_block();

        BOOST_LOG_TRIVIAL(info) << "Cleaving samples: " << nsamples_read
                                << " to " << nsamples_read + gulp_samples;

        std::size_t gulp_size =
            gulp_samples * _config.ndms() * _config.nbeams();

        int nthreads_read = _config.nthreads() > _config.nbridges()
                                ? _config.nbridges()
                                : _config.nthreads();

        omp_set_num_threads(nthreads_read);

        _timer.start("skyweaver::read_data");

        std::vector<bool> read_status(
            _available_freqs.size(),
            false); // since we cannot throw exceptions in parallel regions
#pragma omp parallel for
        for(std::size_t i = 0; i < _available_freqs.size(); i++) {
            FreqType freq = _available_freqs[i];
            if(_bridge_readers.find(freq) == _bridge_readers.end()) {
                read_status[i] = true;
            }
            const auto& reader = _bridge_readers[freq];
            if(reader->eof()) {
                read_status[i] = true;
            }

            std::streamsize read_size =
                reader->read(reinterpret_cast<char*>(thrust::raw_pointer_cast(
                                 _bridge_data[freq]->data())),
                             gulp_size); // read a big chunk of data
            BOOST_LOG_TRIVIAL(debug)
                << "Read " << read_size << " bytes from bridge" << freq;
            if(read_size < gulp_size * sizeof(InputVectorType::value_type)) {
                BOOST_LOG_TRIVIAL(warning)
                    << "Read less data than expected from bridge " << freq;
                read_status[i] = true;
            }
        }

        if(std::any_of(read_status.begin(), read_status.end(), [](bool status) {
               return status;
           })) {
            std::runtime_error("Some bridges have had unexpected reads");
        }

        BOOST_LOG_TRIVIAL(info) << "Read data from bridge readers";

        _timer.stop("skyweaver::read_data");
        _timer.start("skyweaver::process_data");

    /** Process data **/

        omp_set_num_threads(_config.nthreads());

        std::size_t nbridges           = _config.nbridges();
        std::size_t nsamples_per_block = _config.nsamples_per_block();
        std::size_t ndms               = _config.ndms();
        std::size_t nbeams             = _config.nbeams();

#pragma omp parallel for schedule(static) collapse(2)
        for(std::size_t ibeam = 0; ibeam < nbeams; ibeam++) {
            for(std::size_t idm = 0; idm < ndms; idm++) {
                // _beam_data[idm][ibeam] is not found, skip
                if(_beam_data.find(idm) == _beam_data.end()) {
                    continue;
                }
                if(_beam_data[idm].find(ibeam) == _beam_data[idm].end()) {
                    continue;
                }

#pragma omp simd
                for(std::size_t isample = 0; isample < nsamples_per_block; isample++) {
                    const std::size_t base_index =
                        isample * ndms * nbeams + idm * nbeams + ibeam;

                    std::size_t ifreq      = 0;
                    const std::size_t out_offset = isample * nbridges;
                    for(const auto& [freq, ifreq_data]:
                        _bridge_data) { // for each frequency
                        _beam_data[idm][ibeam]->at(out_offset + nbridges - 1 -
                                                   ifreq) = clamp<uint8_t>(127 + ifreq_data->at(base_index));
                        ++ifreq;
                    }
                }
            }
        }
        BOOST_LOG_TRIVIAL(info) << "Processed data";
        _timer.stop("skyweaver::process_data");

    /** Write data **/


        _timer.start("skyweaver::write_data");

#pragma omp parallel for schedule(static) collapse(2)
        for(int idm = 0; idm < _config.ndms(); idm++) {
            for(int ibeam = 0; ibeam < _config.nbeams(); ibeam++) {

                if(_beam_data.find(idm) == _beam_data.end()) {
                    continue;
                }
                if(_beam_data[idm].find(ibeam) == _beam_data[idm].end()) {
                    continue;
                }
                _beam_writers[idm][ibeam]->write(*_beam_data[idm][ibeam],
                                                 _config.stream_id());
            }
        }

        _timer.stop("skyweaver::write_data");
    }
    _timer.show_all_timings();
}
