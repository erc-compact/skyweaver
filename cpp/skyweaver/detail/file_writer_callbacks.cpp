

#include "skyweaver/FileOutputStream.hpp"
#include "skyweaver/MultiFileWriter.cuh"
#include "skyweaver/ObservationHeader.hpp"
#include "skyweaver/SigprocHeader.hpp"

#include <boost/log/trivial.hpp>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
/**
 * The expectation is that each file contains data in
 * TBTF order. Need to explicitly update:
 * INNER_T --> the number of timesamples per block that was processed
 * NBEAMS --> the total number of beams in the file
 * DM --> the dispersion measure of this data
 * Freq --> the centre frequency of this subband (???)
 * BW --> the bandwidth of the subband (???)
 * TSAMP --> the sampling interval of the data
 * NPOL --> normally 1 but could be 4 for full stokes
 * CHAN0_IDX --> this comes from the obs header and uniquely identifies the
 * bridge
 */
std::string default_dada_header = R"(
HEADER       DADA
HDR_VERSION  1.0
HDR_SIZE     4096
DADA_VERSION 1.0

FILE_SIZE    100000000000
FILE_NUMBER  0

UTC_START    1708082229.000020336 
MJD_START    60356.47024305579093

SOURCE       J1644-4559
RA           16:44:49.27
DEC          -45:59:09.7
TELESCOPE    MeerKAT
INSTRUMENT   CBF-Feng
RECEIVER     L-band
FREQ         1284000000.000000
BW           856000000.000000
TSAMP        4.7850467290
STOKES       I

NBIT         8
NDIM         1
NPOL         1
NCHAN        64
NBEAM       800
ORDER        TFB

CHAN0_IDX 2688
  )";

#define MJD_UNIX_EPOCH 40587.0
} // namespace
namespace skyweaver
{
namespace detail
{
template <typename VectorType>
inline std::unique_ptr<FileOutputStream>
create_dada_file_stream(MultiFileWriterConfig const& config,
                        ObservationHeader const& header,
                        VectorType const& stream_data,
                        std::size_t stream_idx)
{
    BOOST_LOG_TRIVIAL(debug) << "Creating stream based on stream prototype: "
                             << stream_data.describe();
    // Here we round the file size to a multiple of the stream prototype
    std::size_t filesize =
        std::max(1ul,
                 config.max_file_size / stream_data.size() /
                     sizeof(typename VectorType::value_type)) *
        stream_data.size();

    BOOST_LOG_TRIVIAL(debug)
        << "Maximum allowed file size = " << filesize << " bytes (+header)";

    std::stringstream output_dir;
    output_dir << config.output_dir << "/" << std::fixed << std::setfill('0')
               << std::setw(9) << static_cast<int>(header.frequency);

    std::stringstream output_basename;
    output_basename << config.output_basename << "_" << std::fixed
                    << std::setfill('0') << std::setw(9)
                    << static_cast<int>(header.frequency);

    std::unique_ptr<FileOutputStream> file_stream =
        std::make_unique<FileOutputStream>(
            output_dir.str(),
            output_basename.str(),
            config.extension,
            filesize,
            [&, header, stream_data, stream_idx, filesize](
                std::size_t& header_size,
                std::size_t bytes_written,
                std::size_t file_idx) -> std::shared_ptr<char const> {
                header_size       = config.header_size;
                char* temp_header = new char[header_size];
                std::fill(temp_header, temp_header + header_size, 0);
                std::memcpy(temp_header,
                            default_dada_header.c_str(),
                            default_dada_header.size());
                psrdada_cpp::RawBytes bytes(temp_header,
                                            header_size,
                                            header_size,
                                            false);
                Header header_writer(bytes);
                header_writer.set<std::string>("SOURCE", header.source_name);
                header_writer.set<std::string>("RA", header.ra);
                header_writer.set<std::string>("DEC", header.dec);
                header_writer.set<std::size_t>("NBEAM", stream_data.nbeams());
                header_writer.set<std::size_t>("NCHAN",
                                               stream_data.nchannels());
                header_writer.set<std::size_t>("OBS_NCHAN", header.obs_nchans);
                header_writer.set<long double>("OBS_FREQUENCY",
                                               header.obs_frequency);
                header_writer.set<long double>("OBS_BW", header.obs_bandwidth);
                header_writer.set<std::size_t>("NSAMP", stream_data.nsamples());
                if(stream_data.ndms()) {
                    header_writer.set<std::size_t>("NDMS", stream_data.ndms());
                    header_writer.set("DMS", stream_data.dms(), 7);
                }
                header_writer.set<long double>(
                    "COHERENT_DM",
                    static_cast<long double>(stream_data.reference_dm()));
                try {
                    header_writer.set<long double>(
                        "FREQ",
                        std::accumulate(stream_data.frequencies().begin(),
                                        stream_data.frequencies().end(),
                                        0.0) /
                            stream_data.frequencies().size());
                } catch(std::runtime_error&) {
                    BOOST_LOG_TRIVIAL(warning)
                        << "Warning: Frequencies array was stale, using the "
                           "centre frequency from the header";
                    header_writer.set<long double>("FREQ", header.frequency);
                }

                header_writer.set<long double>("BW", header.bandwidth);
                header_writer.set<long double>("TSAMP",
                                               stream_data.tsamp() * 1e6);
                if(config.stokes_mode == "IQUV") {
                    header_writer.set<std::size_t>("NPOL", 4);
                } else {
                    header_writer.set<std::size_t>("NPOL", 1);
                }
                header_writer.set<std::string>("STOKES_MODE",
                                               config.stokes_mode);
                header_writer.set<std::string>("ORDER",
                                               stream_data.dims_as_string());
                header_writer.set<std::size_t>("CHAN0_IDX", header.chan0_idx);
                header_writer.set<std::size_t>("FILE_SIZE", filesize);
                header_writer.set<std::size_t>("FILE_NUMBER", file_idx);
                header_writer.set<std::size_t>("OBS_OFFSET", bytes_written);
                header_writer.set<std::size_t>("OBS_OVERLAP", 0);

                long double tstart =
                    header.utc_start + stream_data.utc_offset();

                header_writer.set<long double>("UTC_START", tstart);
                header_writer.set<long double>("MJD_START",
                                               MJD_UNIX_EPOCH +
                                                   tstart / 86400.0);
                std::shared_ptr<char const> header_ptr(
                    temp_header,
                    std::default_delete<char[]>());

                return header_ptr;
            });
    return file_stream;
}

template <typename VectorType>
inline std::unique_ptr<FileOutputStream>
create_sigproc_file_stream(MultiFileWriterConfig const& config,
                           ObservationHeader const& obs_header,
                           VectorType const& stream_data,
                           std::size_t stream_idx)
{
    BOOST_LOG_TRIVIAL(debug) << "Creating stream based on stream prototype: "
                             << stream_data.describe();

    ObservationHeader header = obs_header;

    BOOST_LOG_TRIVIAL(info) << "Header: " << header.to_string();

    // Here we round the file size to a multiple of the stream prototype
    std::size_t filesize =
        std::max(1ul,
                 config.max_file_size / stream_data.size() /
                     sizeof(typename VectorType::value_type)) *
        stream_data.size();
    BOOST_LOG_TRIVIAL(debug)
        << "Maximum allowed file size = " << filesize << " bytes (+header)";

    double foff = -1 *
                  static_cast<double>(header.obs_bandwidth / header.nchans) /
                  1e6; // MHz
    // 1* foff instead of 0.5* foff below because the dedispersion causes all
    // the frequencies to change by half the bandwidth to refer to the bottom of
    // the channel
    double fch1 =
        static_cast<double>(header.obs_frequency + header.obs_bandwidth / 2.0) /
            1e6 +
        foff; // MHz

    double utc_start = static_cast<double>(header.utc_start);
    header.mjd_start = (utc_start / 86400.0) + MJD_UNIX_EPOCH;

    uint32_t datatype    = 0;
    uint32_t barycentric = 0;
    // uint32_t ibeam = 0;
    double az     = 0.0;
    double za     = 0.0;
    uint32_t nifs = header.npol;

    header.sigproc_params = true;
    header.rawfile        = std::string("unset");
    header.fch1           = fch1;
    header.foff           = foff;
    header.tsamp          = header.tsamp / 1e6;
    header.az             = az;
    header.za             = za;
    header.datatype       = datatype;
    header.barycentric    = barycentric;
    header.nifs           = nifs;
    header.telescopeid    = 64;

    BOOST_LOG_TRIVIAL(info) << "Creating Sigproc file stream";

    // // Here we should update the tstart of the default header to be the
    // // start of the stream
    // // reset the total bytes counter to keep the time tracked correctly
    // // Generate the new base filename in <utc>_<tag> format
    // std::stringstream base_filename;
    // // Get UTC time string
    // std::time_t unix_time =
    //     static_cast<std::time_t>(( header.mjd_start  - 40587.0) * 86400.0);
    // struct std::tm* ptm = std::gmtime(&unix_time);

    // // Swapped out put_time call for strftime due to put_time
    // // causing compiler bugs prior to g++ 5.x
    // char formatted_time[80];
    // strftime(formatted_time, 80, "%Y-%m-%d-%H:%M:%S", ptm);
    // base_filename << formatted_time;

    // make config.output_dir if it does not exist

    std::unique_ptr<FileOutputStream> file_stream =
        std::make_unique<FileOutputStream>(
            config.output_dir,
            config.output_basename,
            config.extension,
            filesize,
            [header](std::size_t& header_size,
                     std::size_t bytes_written,
                     std::size_t file_idx) -> std::shared_ptr<char const> {
                // We do not explicitly delete[] this array
                // Cleanup is handled by the shared pointer
                // created below
                std::ostringstream header_stream;
                // get ostream from temp_header

                SigprocHeader sigproc_header(header);
                double mjd_offset = (((bytes_written / (header.nbits / 8.0)) /
                                      (header.nchans)) *
                                     header.tsamp) /
                                    (86400.0);
                sigproc_header.add_time_offset(mjd_offset);
                sigproc_header.write_header(header_stream);
                std::string header_str = header_stream.str();
                header_size            = header_str.size();
                char* header_cstr      = new char[header_size];
                std::copy(header_str.begin(), header_str.end(), header_cstr);
                std::shared_ptr<char const> header_ptr(
                    header_cstr,
                    std::default_delete<char[]>());
                return header_ptr;
            });
    return file_stream;
}

} // namespace detail
} // namespace skyweaver