#include "skyweaver/PipelineConfig.hpp"

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <limits>

namespace skyweaver
{

PipelineConfig::PipelineConfig()
    : _input_files({}), _check_input_contiguity(false), _dada_header_size(4096),
      _delay_file(""), _output_dir("./"), _max_output_filesize(10000000000000),
      _output_file_prefix(""), 
      _enable_incoherent_dedispersion(true), _cfreq(1284000000.0),
      _bw(13375000.0), _channel_frequencies_stale(true),
      _gulp_length_samps(4096), _start_time(0.0f),
      _duration(std::numeric_limits<float>::infinity()), _total_nchans(4096),
      _stokes_mode("I"), _output_level(24.0f),
      _nbeams_per_file(SKYWEAVER_NBEAMS)
{
}

PipelineConfig::~PipelineConfig()
{
}

void PipelineConfig::delay_file(std::string const& fpath)
{
    _delay_file = fpath;
}

std::string const& PipelineConfig::delay_file() const
{
    return _delay_file;
}

void PipelineConfig::input_files(std::vector<std::string> const& files)
{
    _input_files = files;
}

std::vector<std::string> const& PipelineConfig::input_files() const
{
    return _input_files;
}

void PipelineConfig::check_input_contiguity(bool check)
{
    _check_input_contiguity = check;
}

bool PipelineConfig::check_input_contiguity() const
{
    return _check_input_contiguity;
}

void PipelineConfig::dada_header_size(std::size_t size)
{
    _dada_header_size = size;
}

std::size_t PipelineConfig::dada_header_size() const
{
    return _dada_header_size;
}

void PipelineConfig::read_input_file_list(std::string filename)
{
    BOOST_LOG_TRIVIAL(debug) << "Reading input file list from " << filename;
    std::string line;
    std::ifstream ifs(filename.c_str());
    if(!ifs.is_open()) {
        std::cerr << "Unable to open input file list: " << filename << " ("
                  << std::strerror(errno) << ")\n";
        throw std::runtime_error(std::strerror(errno));
    }
    _input_files.resize(0);
    // TODO: Check that the files from the input list exist
    while(std::getline(ifs, line)) {
        BOOST_LOG_TRIVIAL(debug) << line;
        boost::algorithm::trim(line);
        if(line[0] == '#' || line.empty()) {
            // Line is a comment
            BOOST_LOG_TRIVIAL(debug) << "is a comment or empty";
            continue;
        } else {
            _input_files.push_back(line);
        }
        BOOST_LOG_TRIVIAL(debug) << "trimmed: " << line;
    }
    ifs.close();
}

void PipelineConfig::output_dir(std::string const& path)
{
    _output_dir = path;
}

std::string const& PipelineConfig::output_dir() const
{
    return _output_dir;
}

std::size_t PipelineConfig::max_output_filesize() const
{
    return _max_output_filesize;
}

void PipelineConfig::max_output_filesize(std::size_t nbytes)
{
    _max_output_filesize = nbytes;
}

std::string const& PipelineConfig::output_file_prefix() const
{
    return _output_file_prefix;
}

void PipelineConfig::output_file_prefix(std::string const& prefix)
{
    _output_file_prefix = prefix;
}

double PipelineConfig::centre_frequency() const
{
    return _cfreq;
}

void PipelineConfig::centre_frequency(double cfreq)
{
    _cfreq                     = cfreq;
    _channel_frequencies_stale = true;
}

double PipelineConfig::bandwidth() const
{
    return _bw;
}

void PipelineConfig::bandwidth(double bw)
{
    _bw                        = bw;
    _channel_frequencies_stale = true;
}

std::vector<float> const& PipelineConfig::coherent_dms() const
{
    return _ddplan.coherent_dms();
}

DedispersionPlan& PipelineConfig::ddplan()
{
    return _ddplan;
}

DedispersionPlan const& PipelineConfig::ddplan() const
{
    return _ddplan;
}



void PipelineConfig::enable_incoherent_dedispersion(bool enable)
{
    _enable_incoherent_dedispersion = enable;
}

bool PipelineConfig::enable_incoherent_dedispersion() const
{
    return _enable_incoherent_dedispersion;
}

std::vector<double> const& PipelineConfig::channel_frequencies() const
{
    if(_channel_frequencies_stale) {
        calculate_channel_frequencies();
    }
    return _channel_frequencies;
}

void PipelineConfig::calculate_channel_frequencies() const
{
    /**
     * Need to revisit this implementation as it is not clear how the
     * frequencies are labeled for the data out of the F-engine. Either
     * way is a roughly correct place-holder.
     */
    double chbw    = bandwidth() / nchans();
    double fbottom = centre_frequency() - bandwidth() / 2.0;
    _channel_frequencies.clear();
    for(std::size_t chan_idx = 0; chan_idx < nchans(); ++chan_idx) {
        _channel_frequencies.push_back(fbottom + chbw / 2.0 +
                                       (chbw * chan_idx));
    }
    _channel_frequencies_stale = false;
}

std::size_t PipelineConfig::gulp_length_samps() const
{
    return _gulp_length_samps;
}

void PipelineConfig::gulp_length_samps(std::size_t nsamps)
{
    _gulp_length_samps = nsamps;
}

float PipelineConfig::start_time() const
{
    return _start_time;
}

void PipelineConfig::start_time(float start_time_)
{
    _start_time = start_time_;
}

float PipelineConfig::duration() const
{
    return _duration;
}

void PipelineConfig::duration(float duration_)
{
    _duration = duration_;
}

void PipelineConfig::output_level(float level)
{
    _output_level = level;
}

float PipelineConfig::output_level() const
{
    return _output_level;
}

std::size_t PipelineConfig::total_nchans() const
{
    return _total_nchans;
}

void PipelineConfig::total_nchans(std::size_t nchans)
{
    _total_nchans = nchans;
}

std::size_t PipelineConfig::nbeams_per_file() const
{
    return _nbeams_per_file;
}

void PipelineConfig::nbeams_per_file(std::size_t nbeams_per_file)
{
    _nbeams_per_file = nbeams_per_file;
}

} // namespace skyweaver
