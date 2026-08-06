#include "skyweaver/Header.hpp"
#include "skyweaver/MultiFileWriter.cuh"

#include <cstring>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

/**
 * Now write a DADA file per DM
 * with optional time splitting
 */

namespace skyweaver
{

namespace
{

std::string get_formatted_time(long double unix_timestamp)
{
    char formatted_time[80];
    std::time_t unix_time = unix_timestamp;
    struct std::tm* ptm   = std::gmtime(&unix_time);
    strftime(formatted_time, sizeof(formatted_time), "%Y-%m-%d-%H:%M:%S", ptm);
    return std::string(formatted_time);
}

} // namespace

// template <typename VectorType>
// MultiFileWriter<VectorType>::MultiFileWriter(PipelineConfig const& config,
//                                              std::string tag)
//     : _config(config), _tag(tag)
// {
// }

template <typename VectorType>
MultiFileWriter<VectorType>::MultiFileWriter(
    PipelineConfig const& config,
    std::string tag,
    CreateStreamCallBackType create_stream_callback)
    : _tag(tag), _create_stream_callback(create_stream_callback)
{
    _config.header_size     = config.dada_header_size();
    _config.max_file_size   = config.max_output_filesize();
    _config.stokes_mode     = config.stokes_mode();
    _config.base_output_dir = config.output_dir();
    _config.inner_dir = "";

}

template <typename VectorType>
MultiFileWriter<VectorType>::MultiFileWriter(
    MultiFileWriterConfig config,
    std::string tag,
    CreateStreamCallBackType create_stream_callback)
    : _config(config), _tag(tag),
      _create_stream_callback(create_stream_callback)
{
}

template <typename VectorType>
MultiFileWriter<VectorType>::~MultiFileWriter(){};

template <typename VectorType>
void MultiFileWriter<VectorType>::init(ObservationHeader const& header)
{
    _header = header; // Make copy of the header
}

template <typename VectorType>
bool MultiFileWriter<VectorType>::has_stream(std::size_t stream_idx)
{
    auto it = _file_streams.find(stream_idx);
    return (it != _file_streams.end());
}

template <typename VectorType>
FileOutputStream&
MultiFileWriter<VectorType>::create_stream(VectorType const& stream_data,
                                           std::size_t stream_idx)
{
    _config.output_dir = get_output_dir(stream_data, stream_idx);

    if(_config.extension.empty()) {
        _config.extension = get_extension(stream_data);
    }

    _config.output_basename = get_basefilename(stream_data, stream_idx);

    _file_streams[stream_idx] =
        _create_stream_callback(_config, _header, stream_data, stream_idx);

    return *_file_streams[stream_idx];
}

template <typename VectorType>
std::string
MultiFileWriter<VectorType>::get_output_dir(VectorType const& stream_data,
                                            std::size_t stream_idx)
{
    // Output directory format
    // <utcstart>/<freq:%f.02>/<stream_id>
    std::stringstream output_dir;
    output_dir << _config.base_output_dir << "/"
               << get_formatted_time(_header.utc_start) << "/" << stream_idx;

    if(!_config.inner_dir.empty()) {
        output_dir << "/" << _config.inner_dir;
    }

    return output_dir.str();
}

template <typename VectorType>
std::string
MultiFileWriter<VectorType>::get_basefilename(VectorType const& stream_data,
                                              std::size_t stream_idx)
{
    // Output file format
    // <prefix>_<utcstart>_<dm:%f.03>_<byte_offset>.<extension>
    std::stringstream base_filename;
    if(!_config.prefix.empty()) {
        base_filename << _config.prefix << "_";
    }
    base_filename << get_formatted_time(_header.utc_start) << "_" << stream_idx
                  << "_cdm_" << std::fixed << std::setprecision(3)
                  << std::setfill('0') << std::setw(9)
                  << stream_data.reference_dm(); 

    if(!_config.suffix.empty()) {
        base_filename << "_" << _config.suffix;
    }

    if(!_tag.empty()) {
        base_filename << "_" << _tag;
    }
    return base_filename.str();
}

template <typename VectorType>
std::string
MultiFileWriter<VectorType>::get_extension(VectorType const& stream_data)
{
    std::string dims = stream_data.dims_as_string();
    for(auto& c: dims) { c = std::tolower(static_cast<unsigned char>(c)); }
    if(dims == "t") {
        return ".dat";
    } else if(dims == "tf") {
        return ".fil";
    }
    return "." + dims;
}

template <typename VectorType>
bool MultiFileWriter<VectorType>::operator()(VectorType const& stream_data,
                                             std::size_t stream_idx)
{
    if(!has_stream(stream_idx)) {
        create_stream(stream_data, stream_idx);
    }
    if constexpr(is_device_vector<VectorType>::value) {
        thrust::host_vector<typename VectorType::value_type> stream_data_h =
            stream_data.vector();
        _file_streams.at(stream_idx)
            ->write(reinterpret_cast<char const*>(stream_data_h.data()),
                    stream_data_h.size() *
                        sizeof(typename VectorType::value_type));
    } else {
        _file_streams.at(stream_idx)
            ->write(reinterpret_cast<char const*>(
                        thrust::raw_pointer_cast(stream_data.data())),
                    stream_data.size() *
                        sizeof(typename VectorType::value_type));
    }
    return false;
}

template <typename VectorType>
bool MultiFileWriter<VectorType>::write(VectorType const& stream_data,
                                        std::size_t stream_idx)
{
    return this->operator()(stream_data, stream_idx);
}

} // namespace skyweaver