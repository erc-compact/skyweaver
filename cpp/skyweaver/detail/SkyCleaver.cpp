
#include "skyweaver/SkyCleaver.hpp"

#include "skyweaver/beamformer_utils.cuh"
#include "skyweaver/types.cuh"
#include "skyweaver/skycleaver_utils.hpp"


#include <cmath>
#include <filesystem>
#include <iostream>
#include <omp.h>
#include <ranges>
#include <regex>
#include <type_traits>

namespace fs = std::filesystem;

using BridgeReader    = skyweaver::BridgeReader;
using MultiFileReader = skyweaver::MultiFileReader;
using BeamInfo        = skyweaver::BeamInfo;

namespace
{

std::string trim(const std::string& str) {
    auto start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) {
        return ""; // String is all whitespace
    }
    auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

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
                        //check if .tmp not in filename
                        if(file_name.find(".tmp") != std::string::npos) {
                            continue;
                        }
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

    std::sort(files.begin(), files.end());

    return files;
}


} // namespace

void parse_target_file(std::string file_name, std::vector<BeamInfo>& beam_infos){
    std::ifstream targets_file(file_name);
    if(!targets_file.is_open()) {
        std::runtime_error("Error opening target file: " + file_name);
    }
    // the file is in csv format. First read the header to know the positions of name, ra and dec
    std::string header;
    do{
        std::getline(targets_file, header);
    } while(header.empty() || header.find("#") != std::string::npos);

    std::vector<std::string> header_tokens;
    std::stringstream header_stream(header);
    std::string token;
    while(std::getline(header_stream, token, ',')) {
        header_tokens.push_back(token);
    }
    std::size_t name_pos = std::distance(header_tokens.begin(), std::find(header_tokens.begin(), header_tokens.end(), "name"));
    std::size_t ra_pos = std::distance(header_tokens.begin(), std::find(header_tokens.begin(), header_tokens.end(), "ra"));
    std::size_t dec_pos = std::distance(header_tokens.begin(), std::find(header_tokens.begin(), header_tokens.end(), "dec"));

    if(name_pos == header_tokens.size() || ra_pos == header_tokens.size() || dec_pos == header_tokens.size()) {
        std::runtime_error("Invalid header in target file: " + file_name);
    }

    std::string line;
    while(std::getline(targets_file, line)) {

        line = trim(line);

        //if empty line or # anywhere in line, continue
        if(line.empty() || line.find("#") != std::string::npos) {
            BOOST_LOG_TRIVIAL(debug) << "Skipping line: " << line;
            continue;
        }
   
        std::vector<std::string> tokens;
        std::stringstream line_stream(line);
        std::string token;
        while(std::getline(line_stream, token, ',')) {
            tokens.push_back(token);
        }
        if(tokens.size() != header_tokens.size()) {
            std::runtime_error("Invalid number of columns in target file: " + file_name);
        }
        BeamInfo beam_info;
        beam_info.beam_name = tokens[name_pos];
        beam_info.beam_ra = tokens[ra_pos];
        beam_info.beam_dec = tokens[dec_pos];
        beam_infos.push_back(beam_info);
    }
}

void compare_bridge_headers(const skyweaver::ObservationHeader& first,
                            const skyweaver::ObservationHeader& second)
{
    if(first.nchans != second.nchans) {
        throw std::runtime_error("Number of channels in bridge readers do not "
                                 "match. Expected: " +
                                 std::to_string(first.nchans) +
                                 " Found: " + std::to_string(second.nchans));
    }
    if(first.nbeams != second.nbeams) {
        throw std::runtime_error(
            "Number of beams in bridge readers do not match. "
            "Expected: " +
            std::to_string(first.nbeams) +
            " Found: " + std::to_string(second.nbeams));
    }
    if(first.nbits != second.nbits) {
        throw std::runtime_error(
            "Number of bits in bridge readers do not match. "
            "Expected: " +
            std::to_string(first.nbits) +
            " Found: " + std::to_string(second.nbits));
    }
    if(first.tsamp != second.tsamp) {
        throw std::runtime_error(
            "Sampling time in bridge readers do not match. "
            "Expected: " +
            std::to_string(first.tsamp) +
            " Found: " + std::to_string(second.tsamp));
    }
    if(first.stokes_mode != second.stokes_mode) {
        throw std::runtime_error("Stokes mode in bridge readers do not match. "
                                 "Expected: " +
                                 first.stokes_mode +
                                 " Found: " + second.stokes_mode);
    }
}

template <typename InputVectorType, typename OutputVectorType>
void skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::init_readers()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::init_readers")

    std::string root_dir  = _config.root_dir();
    std::size_t stream_id = _config.stream_id();

    // get the list of directories in root/stream_id(for the nex)
    std::vector<std::string> freq_dirs =
        get_subdirs(root_dir + "/" + std::to_string(stream_id));

    BOOST_LOG_TRIVIAL(info)
        << "Found " << freq_dirs.size()
        << " frequency directories in root directory: " << root_dir;

    std::map<skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::FreqType,
             long double>
        bridge_timestamps;
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

    _header = _bridge_readers[_available_freqs[0]]->get_header();
    for(const auto& [freq, reader]: _bridge_readers) {
        compare_bridge_headers(_header, reader->get_header());
    }
    BOOST_LOG_TRIVIAL(info)
        << "Number of beams: " << _header.nbeams
        << " Number of DMS: " << _header.ndms
        << " Stokes mode: " << _header.stokes_mode
        << " Number of channels: " << _header.nchans;


    _config.nbeams(_header.nbeams);
    _config.ndms(_header.ndms);
    _config.stokes_mode(_header.stokes_mode);
    _config.nchans(_header.nchans);

    BOOST_LOG_TRIVIAL(info)
        << "Number of beams: " << _config.nbeams()
        << " Number of DMS: " << _config.ndms()
        << " Stokes mode: " << _config.stokes_mode()
        << " Number of channels: " << _config.nchans();

    std::vector<std::vector<std::size_t>> stokes_positions;
    for(const auto stokes: _config.out_stokes()) {
        std::size_t pos = _config.stokes_mode().find(stokes);
        if(pos == std::string::npos) {

            if(stokes == 'L') {
                std::size_t pos1 = _config.stokes_mode().find("Q");
                std::size_t pos2 = _config.stokes_mode().find("U");
                if(pos1 == std::string::npos || pos2 == std::string::npos) {
                    throw std::runtime_error("Asked for L, but beamformed data does not have Q and/or U");
                }
                stokes_positions.push_back({pos1, pos2});
                continue;
            }
            else {
                throw std::runtime_error("Requested stokes not found in beamformed data: " + stokes);
            }
        }
        stokes_positions.push_back({pos});
    }

    _config.stokes_positions(stokes_positions);


    long double obs_centre_freq = _header.obs_frequency;
    long double obs_bandwidth   = _header.obs_bandwidth;

    long double start_freq = obs_centre_freq - obs_bandwidth / 2;

    for(int i = 0; i < _config.nbridges(); i++) {
      std::size_t ifreq = std::lround(std::floor(
            start_freq + (i + 0.5) * obs_bandwidth / _config.nbridges()));
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
                                  _config.nbeams(),
                                  _config.stokes_mode().size()},
            0);
    }

    std::size_t smallest_data_size = std::numeric_limits<std::size_t>::max();
    std::size_t dbp_factor =
        _config.ndms() * _config.nbeams() * _config.stokes_mode().size();

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
            << "Seeking " << nsamples * dbp_factor
            << " bytes in bridge reader for frequency: " << freq;

        std::size_t bytes_seeking =
            (nsamples * dbp_factor *
             sizeof(typename InputVectorType::value_type));

        _bridge_readers[freq]->seekg(bytes_seeking, std::ios_base::beg);

        std::size_t data_size =
            _bridge_readers[freq]->get_total_size() - bytes_seeking;
        BOOST_LOG_TRIVIAL(debug)
            << "Data size for frequency: " << freq << " is " << data_size;
        if(data_size < smallest_data_size) {
            smallest_data_size = data_size;
        }
    }

    if(smallest_data_size % dbp_factor != 0) {
        std::runtime_error(
            "Data size is not a multiple of ndms * nbeams * nstokes");
    }

    std::size_t smallest_nsamples = std::floor(smallest_data_size / dbp_factor);

    if(smallest_nsamples < _config.start_sample()) {
        std::runtime_error(
            "start_sample is greater than the smallest_nsamples in the data.");
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
    } else {
        _nsamples_to_read = smallest_nsamples;
    }

    std::size_t bytes_seeking = (_config.start_sample() * dbp_factor *
                                 sizeof(typename InputVectorType::value_type));

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
    BOOST_LOG_TRIVIAL(info)
        << "Adding first header to SkyCleaver: " << _header.to_string();
    _header.nchans = _header.nchans * _config.nbridges();
    _header.nbeams = _config.nbeams();

    _nthreads_read = _config.nthreads() > _config.nbridges()
                         ? _config.nbridges()
                         : _config.nthreads();
}
template <typename InputVectorType, typename OutputVectorType>
void skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::init_writers()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::init_writers")

    BOOST_LOG_TRIVIAL(debug)
        << "_config.output_dir(); " << _config.output_dir();
    std::string output_dir = _config.output_dir();

    if(!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
    }

    for(std::size_t istokes = 0; istokes < _config.out_stokes().size();
        istokes++) {
        for(int idm = 0; idm < _config.ndms(); idm++) {
            if(!_config.required_dms().empty()) {
                const auto& required_dms = _config.required_dms();
                if(std::ranges::find(required_dms, _header.dms[idm]) ==
                   required_dms.end()) {
                    BOOST_LOG_TRIVIAL(info)
                        << "DM " << _header.dms[idm]
                        << " is not required, skipping from writing";
                    continue;
                }
            }
            for(int ibeam = 0; ibeam < _config.nbeams(); ibeam++) {
                // skip if beam is not used
                if(!_config.required_beams().empty()) {
                    const auto& required_beams = _config.required_beams();
                    std::cerr << "required_beams: " << required_beams.size()
                              << "ibeam: " << ibeam << std::endl;
                    if(std::ranges::find(required_beams, ibeam) ==
                       required_beams.end()) {
                        BOOST_LOG_TRIVIAL(info)
                            << "Beam " << ibeam
                            << " is not required, skipping from writing";
                        continue;
                    }
                }
                BeamInfo beam_info = _beam_infos[ibeam];
                MultiFileWriterConfig writer_config;
                writer_config.header_size   = _config.dada_header_size();
                writer_config.max_file_size = _config.max_output_filesize();
                writer_config.stokes_mode   = _config.out_stokes().at(istokes);
                writer_config.base_output_dir = output_dir;
                writer_config.inner_dir = beam_info.beam_name;
                writer_config.prefix          = _config.out_prefix();
                std::string suffix = "idm_" +
                              to_string_with_padding(_header.dms[idm], 9, 3);
                    
                writer_config.suffix = suffix + "_" +
                                       beam_info.beam_name + "_" +
                                       _config.out_stokes().at(istokes);
                writer_config.extension = ".fil";

                BOOST_LOG_TRIVIAL(info)
                    << "Writer config: " << writer_config.to_string();

                typename MultiFileWriter<OutputVectorType>::
                    CreateStreamCallBackType create_stream_callback_sigproc =
                        skyweaver::detail::create_sigproc_file_stream<
                            OutputVectorType>;
                _beam_writers[istokes][idm][ibeam] =
                    std::make_unique<MultiFileWriter<OutputVectorType>>(
                        writer_config,
                        "",
                        create_stream_callback_sigproc);
                _header.ibeam = ibeam;

                _header.ra = beam_info.beam_ra;
                _header.dec = beam_info.beam_dec;
                _beam_writers[istokes][idm][ibeam]->init(_header);

                _beam_data[istokes][idm][ibeam] =
                    std::make_shared<OutputVectorType>(
                        std::initializer_list{_config.nsamples_per_block(),
                                              _config.nbridges()},
                        0);

                _beam_data[istokes][idm][ibeam]->reference_dm(_header.refdm);
                _total_beam_writers++;
            }
        }
    }

    BOOST_LOG_TRIVIAL(info)
        << "Added " << _total_beam_writers << " beam writers to SkyCleaver";
}
template <typename InputVectorType, typename OutputVectorType>
skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::SkyCleaver(
    SkyCleaverConfig& config)
    : _config(config)
{
    BOOST_LOG_TRIVIAL(info) << "Reading and initialising beam details from file: "
                        << _config.targets_file();
    
    parse_target_file(_config.targets_file(), _beam_infos);

    BOOST_LOG_TRIVIAL(info) << "Number of beams in target file: " << _beam_infos.size();

    _timer.start("skycleaver::init_readers");
    init_readers();
    _timer.stop("skycleaver::init_readers");
    

    if(_beam_infos.size() < _config.nbeams()){ // there are some null beams with zeros, do not create filterbanks for them. 
        std::string required_beams = "0:" + std::to_string(_beam_infos.size()-1);
        std::vector<int> required_beam_numbers = skyweaver::get_list_from_string<int>(required_beams);

        if(_config.required_beams().empty()) { // if nothing given, set the valid beams to required beams
        BOOST_LOG_TRIVIAL(warning) << "Number of beams in target file is less than the number of beams in the header. "
                                << "Setting required beams to: " << required_beams;            
        _config.required_beams(skyweaver::get_list_from_string<int>(required_beams));
        }
        else{
            for(auto beam_num: _config.required_beams()){ // if given, check if all requested beams are valid beams
                if(std::find(required_beam_numbers.begin(), required_beam_numbers.end(), beam_num) == required_beam_numbers.end()){
                    std::runtime_error("Beam number " + std::to_string(beam_num) + " not found in target file.");
                }
            }
        }
    }
    

    _timer.start("skycleaver::init_writers");
    init_writers();
    _timer.stop("skycleaver::init_writers");
}
template <typename InputVectorType, typename OutputVectorType>
void skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::read(
    std::size_t gulp_samples)
{
    std::size_t gulp_size = gulp_samples * _config.ndms() * _config.nbeams() *
                            _config.stokes_mode().size();
    BOOST_LOG_TRIVIAL(info) << "Reading gulp samples: " << gulp_samples
                            << " with size: " << gulp_size;

    omp_set_num_threads(_nthreads_read);

    std::vector<bool> read_failures(
        _available_freqs.size(),
        false); // since we cannot throw exceptions in parallel regions
#pragma omp parallel for
    for(std::size_t i = 0; i < _available_freqs.size(); i++) {
        skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::FreqType
            freq = _available_freqs[i];
        if(_bridge_readers.find(freq) == _bridge_readers.end()) {
            read_failures[i] = true;
        }
        const auto& reader = _bridge_readers[freq];
        if(reader->eof()) {
            BOOST_LOG_TRIVIAL(warning)
                << "End of file reached for bridge " << freq;
            read_failures[i] = true;
        }

        std::streamsize read_size =
            reader->read(reinterpret_cast<char*>(thrust::raw_pointer_cast(
                             _bridge_data[freq]->data())),
                         gulp_size); // read a big chunk of data
        BOOST_LOG_TRIVIAL(debug)
            << "Read " << read_size << " bytes from bridge" << freq;
        if(read_size <
           gulp_size * sizeof(typename InputVectorType::value_type)) {
            BOOST_LOG_TRIVIAL(warning)
                << "Read less data than expected from bridge " << freq;
            read_failures[i] = true;
        }
    }
    bool failed = false;
    for(int i = 0; i < read_failures.size(); i++) {
        if(read_failures[i]) {
            BOOST_LOG_TRIVIAL(error)
                << "Reading bridge [" << i << "]: failed " << std::endl;
            failed = true;
        }
    }
    if(failed) {
        std::runtime_error("Failed to read data from bridge readers.");
    }

    BOOST_LOG_TRIVIAL(info) << "Read data from bridge readers";
}

template <typename InputVectorType, typename OutputVectorType>
void skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::cleave()
{
    BOOST_LOG_NAMED_SCOPE("SkyCleaver::cleave")

    for(std::size_t nsamples_read = 0; nsamples_read < _nsamples_to_read;
        nsamples_read += _config.nsamples_per_block()) {
        std::size_t gulp_samples =
            _nsamples_to_read - nsamples_read < _config.nsamples_per_block()
                ? _nsamples_to_read - nsamples_read
                : _config.nsamples_per_block();

        BOOST_LOG_TRIVIAL(info) << "Cleaving samples: " << nsamples_read
                                << " to " << nsamples_read + gulp_samples;

        _timer.start("skyweaver::read_data");
        read(gulp_samples);
        _timer.stop("skyweaver::read_data");

        _timer.start("skyweaver::process_data");
        omp_set_num_threads(_config.nthreads());

        std::size_t nbridges    = _config.nbridges();
        std::size_t ndms        = _config.ndms();
        std::size_t nbeams      = _config.nbeams();
        std::size_t nstokes_out = _config.out_stokes().size();
        std::size_t nstokes_in  = _config.stokes_mode().size();

#pragma omp parallel for schedule(static) collapse(3)
        for(std::size_t istokes = 0; istokes < nstokes_out; istokes++) {
            for(std::size_t ibeam = 0; ibeam < nbeams; ibeam++) {
                for(std::size_t idm = 0; idm < ndms;
                    idm++) { // cannot separate loops, so do checks later

                    if(_beam_data.find(istokes) == _beam_data.end()) {
                        continue;
                    }
                    if(_beam_data[istokes].find(idm) ==
                       _beam_data[istokes].end()) {
                        continue;
                    }
                    if(_beam_data[istokes][idm].find(ibeam) ==
                       _beam_data[istokes][idm].end()) {
                        continue;
                    }


                    const std::vector<std::size_t> stokes_positions =
                        _config.stokes_positions()[istokes];
                
#pragma omp simd
                    for(std::size_t isample = 0; isample < gulp_samples;
                        isample++) {
                        const std::size_t out_offset = isample * nbridges;

                        // This is stupid but preferred over a more elegant solution that is not fast, this can be easily vectorised
                        if(stokes_positions.size() == 1) {
                            const std::size_t base_index =
                            isample * ndms * nbeams * nstokes_in +
                            idm * nbeams * nstokes_in + ibeam * nstokes_in +
                            stokes_positions[0];
 
                            std::size_t ifreq            = 0;
                            for(const auto& [freq, ifreq_data]:
                            _bridge_data) { // for each frequency
                            _beam_data[istokes][idm][ibeam]->at(
                                out_offset + nbridges - 1 - ifreq) =
                                clamp<uint8_t>(127 +
                                               ifreq_data->at(base_index));
                            ++ifreq;
                            }
                        }
                        else{
                            std::size_t ifreq            = 0;
                            for(const auto& [freq, ifreq_data]: _bridge_data) { 
                                float value = 0;
                                for(int stokes_position=0; stokes_position<stokes_positions.size(); stokes_position++) {
                                    const std::size_t base_index =
                                        isample * ndms * nbeams * nstokes_in +
                                        idm * nbeams * nstokes_in + ibeam * nstokes_in +
                                        stokes_positions[stokes_position];
                                    value += (ifreq_data->at(base_index) * ifreq_data->at(base_index));
                                }
                                // for each frequency
                                _beam_data[istokes][idm][ibeam]->at(
                                    out_offset + nbridges - 1 - ifreq) =
                                    clamp<uint8_t>(127 +
                                                sqrt(value));
                                ++ifreq;
                            }
                        }

                    }
                        
                }
            }
        }
        BOOST_LOG_TRIVIAL(info) << "Processed data";
        _timer.stop("skyweaver::process_data");
        _timer.start("skyweaver::write_data");
        write();
        _timer.stop("skyweaver::write_data");
    }
    _timer.show_all_timings();

}

template <typename InputVectorType, typename OutputVectorType>
void skyweaver::SkyCleaver<InputVectorType, OutputVectorType>::write()
{
    omp_set_num_threads(_config.nthreads());
#pragma omp parallel for schedule(static) collapse(3)
    for(std::size_t istokes = 0; istokes < _config.out_stokes().size();
        istokes++) {
        for(std::size_t idm = 0; idm < _config.ndms(); idm++) {
            for(std::size_t ibeam = 0; ibeam < _config.nbeams(); ibeam++) {
                if(_beam_data.find(istokes) == _beam_data.end()) {
                    continue;
                }
                if(_beam_data[istokes].find(idm) == _beam_data[istokes].end()) {
                    continue;
                }
                if(_beam_data[istokes][idm].find(ibeam) ==
                   _beam_data[istokes][idm].end()) {
                    continue;
                }
                _beam_writers[istokes][idm][ibeam]->write(
                    *_beam_data[istokes][idm][ibeam],
                    _config.stream_id());
            }
        }
    }
}
