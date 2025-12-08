#ifndef SKYWEAVER_OBSERVATIONHEADER_HPP
#define SKYWEAVER_OBSERVATIONHEADER_HPP

#include "psrdada_cpp/raw_bytes.hpp"
#include "skyweaver/Header.hpp"
#include "skyweaver/PipelineConfig.hpp"

#include <cmath>

namespace skyweaver
{

struct ObservationHeader {
    std::size_t nchans    = 0; // Number of frequency channels in the subband
    std::size_t npol      = 0; // Number of polarisations
    std::size_t nbits     = 0; // Number of bits per sample
    std::size_t nantennas = 0; // Number of antennas
    std::size_t sample_clock_start = 0; // The start epoch in sampler ticks
    std::size_t chan0_idx          = 0; // The index of the first channel
    std::size_t beam0_idx          = 0; // The index of the first beam
    std::size_t obs_nchans =
        0; // The total number of channels in the observation
    long double bandwidth = 0.0; // Bandwidth in Hz of the subband
    long double obs_bandwidth =
        0.0;                     // The full bandwidth in Hz of the observation
    long double frequency = 0.0; // Centre frequency in Hz of the subband
    long double obs_frequency =
        0.0;                        // Centre frequency in Hz of the observation
    long double tsamp        = 0.0; // Sampling interval in microseconds
    long double sample_clock = 0.0; // The sampling rate in Hz
    long double sync_time    = 0.0; // The UNIX epoch of the sampler zero
    long double utc_start    = 0.0; // The UTC start time of the data
    long double mjd_start    = 0.0; // The MJD start time of the data
    std::size_t obs_offset =
        0; // The offset of the current file from UTC_START in bytesß
    long double refdm  = 0.0; // Reference DM
    std::size_t ibeam  = 0.0; // Beam number
    std::size_t nbeams = 0;   // Number of beams
    std::string source_name;  // Name of observation target
    std::string ra;           // Right ascension
    std::string dec;          // Declination
    std::string telescope;    // Telescope name
    std::string instrument;   // Name of the recording instrument
    std::string order;        // Order of the dimensions in the data
    std::size_t ndms;         // Number of DMs
    std::vector<float> dms;   // DMs
    std::string stokes_mode;  // Stokes mode

    std::string to_string() const; // Convert the header to a string

    long double az;              // Azimuth
    long double za;              // Zenith angle
    std::size_t machineid   = 0; // Machine ID
    std::size_t nifs        = 0; // Number of IFs
    std::size_t telescopeid = 0; // Telescope ID
    std::size_t datatype    = 0; // Data type
    std::size_t barycentric = 0; // Barycentric correction
    std::string rawfile;         // Raw file name
    double fch1 = 0.0;           // Centre frequency of the first channel
    double foff = 0.0;           // Frequency offset between channels

    bool sigproc_params = false; // Whether to include sigproc parameters
    ObservationHeader() = default;
    ObservationHeader(ObservationHeader const&)            = default;
    ObservationHeader& operator=(ObservationHeader const&) = default;
};

// template for comparing two floating point objects

template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, bool>::type
is_close(T a, T b, T tolerance = 1e-12)
{
    return std::fabs(a - b) < tolerance;
}

/**
 * @brief Parse header information for a DADA header block
 *
 * @param raw_header A RawBytes object containing the DADA header
 * @param header An ObservationHeader instance
 */
void read_dada_header(psrdada_cpp::RawBytes& raw_header,
                      ObservationHeader& header);

void validate_header(ObservationHeader const& header,
                     PipelineConfig const& config);

void update_config(PipelineConfig& config, ObservationHeader const& header);

bool are_headers_similar(ObservationHeader const& header1,
                         ObservationHeader const& header2);

} // namespace skyweaver

#endif // SKYWEAVER_OBSERVATIONHEADER_HPP
