#include "skyweaver/MultiFileReader.cuh"

#include <boost/log/trivial.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thrust/host_vector.h>
using namespace skyweaver;

std::size_t MultiFileReader::safe_read(std::ifstream& input_stream,
                                       char* buffer,
                                       std::size_t nbytes)
{
    BOOST_LOG_TRIVIAL(debug)
        << "At byte " << input_stream.tellg() << " of the input file";
    BOOST_LOG_TRIVIAL(debug) << "Safe reading " << nbytes << " bytes";
    input_stream.read(buffer, nbytes);
    std::size_t nbytes_read = input_stream.gcount();
    if(input_stream.eof()) {
        // Reached the end of the delay model file
        // TODO: Decide what the behaviour should be here.
        throw std::runtime_error("Reached end of delay model file");
    } else if(input_stream.fail() || input_stream.bad()) {
        std::ostringstream error_msg;
        error_msg << "Error: Unable to read from delay file: "
                  << std::strerror(errno);
        throw std::runtime_error(error_msg.str().c_str());
    }
    BOOST_LOG_TRIVIAL(debug) << "Safe Read complete";
    return nbytes_read;
}

MultiFileReader::MultiFileReader(PipelineConfig const& config)
    : MultiFileReader(config.input_files(),
                      config.dada_header_size(),
                      config.check_input_contiguity()) { }

MultiFileReader::MultiFileReader(std::vector<std::string> dada_files,std::size_t dada_header_size, bool check_input_contiguity)
    : _current_file_idx(0), _current_position(0), _eof_flag(false) {

    static_assert(sizeof(char2) == 2 * sizeof(char),
                  "Size of char2 is not as expected");

    this->_total_size       = 0;
    this->_dada_files       = dada_files;
    this->_dada_header_size = dada_header_size;
    std::vector<char> header_bytes(_dada_header_size);



    if(_dada_files.size() == 0) {
        throw std::runtime_error("No input files provided for MultiFileReader");
    }


    for(const auto& file: _dada_files) {
        std::ifstream f(file, std::ifstream::in | std::ifstream::binary);
        if(!f.is_open()) {
            throw std::runtime_error("Could not open file " + file);
        }
        safe_read(f, header_bytes.data(), _dada_header_size);
        read_header(header_bytes);
        f.seekg(0, std::ios::end); // Move to the end of the file
        std::size_t size = f.tellg();

        if (size < _dada_header_size) {
            throw std::runtime_error("File " + file + " is too small to contain a header");
        }

        _total_size += (size - _dada_header_size); // skip header
        _sizes.push_back(size - _dada_header_size);
        BOOST_LOG_TRIVIAL(debug)
            << "File: " << file << " Binary size:" << size - _dada_header_size;
        f.close();
    }
    BOOST_LOG_TRIVIAL(debug) << "Total size of data: " << _total_size;

    if(check_input_contiguity) {
        check_contiguity();
    }

    open(); // open the first file
}

void MultiFileReader::check_contiguity()
{
    ObservationHeader first_header = _headers[0];
    std::size_t size_so_far = _sizes[0] + first_header.obs_offset;

    for(int i=1; i < _sizes.size(); i++) {
        are_headers_similar(first_header, _headers[i]);
        BOOST_LOG_TRIVIAL(debug) << "Header of " << _dada_files[i-1] << " and " << _dada_files[i] << " are similar";
        BOOST_LOG_TRIVIAL(debug) << "Checking contiguity between them";
        BOOST_LOG_TRIVIAL(debug) << "Offset of " << _dada_files[i-1] << " is " << _headers[i-1].obs_offset;
        BOOST_LOG_TRIVIAL(debug) << "Offset of " << _dada_files[i] << " is " << _headers[i].obs_offset;
        BOOST_LOG_TRIVIAL(debug) << "Size so far: " << size_so_far;
        if(_headers[i].obs_offset != size_so_far) {
            throw std::runtime_error("These are not contiguous:" + _dada_files[i-1] + " and " + _dada_files[i]);
        }
        size_so_far += _sizes[i];
       
    }
    if(size_so_far != _total_size) {
        throw std::runtime_error("Total size does not match the sum of individual sizes");
    }
}

MultiFileReader::~MultiFileReader()
{
    this->close();
}

std::vector<std::string> const& MultiFileReader::get_dada_files() const
{
    return _dada_files;
}

void MultiFileReader::read_header(std::vector<char>& header_bytes)
{
    psrdada_cpp::RawBytes block(header_bytes.data(),
                                _dada_header_size,
                                _dada_header_size,
                                false);

    ObservationHeader header;
    read_dada_header(block, header);
    _headers.push_back(header);

    BOOST_LOG_TRIVIAL(debug) << "Header parameters: " << this->header;
}

skyweaver::ObservationHeader const& MultiFileReader::get_header() const
{
    return _headers[0];
}

skyweaver::ObservationHeader const& MultiFileReader::get_current_header() const
{
    return _headers[_current_file_idx];
}

void MultiFileReader::open()
{
    BOOST_LOG_TRIVIAL(debug)
        << "Opening file: " << _dada_files[_current_file_idx]
        << "of size: " << _sizes[_current_file_idx];
    _current_stream.open(_dada_files[_current_file_idx],
                         std::ifstream::in | std::ifstream::binary);
    if(_current_stream.good()) {
        _current_stream.seekg(
            _dada_header_size,
            std::ios::beg); // Skip header for subsequent _dada_files
    } else {
        throw std::runtime_error("Could not open file " +
                                 _dada_files[_current_file_idx]);
    }
    _current_position = _dada_header_size;
    _eof_flag         = false;
    BOOST_LOG_TRIVIAL(debug)
        << "Current file: " << _dada_files[_current_file_idx];
}

void MultiFileReader::open_next()
{
    BOOST_LOG_TRIVIAL(debug)
        << "Opening next file: " << _dada_files[_current_file_idx]
        << "of size: " << _sizes[_current_file_idx];
    _current_file_idx++;
    if(_current_file_idx < _dada_files.size()) {
        _current_stream.close();
        open();
    } else {
        _eof_flag = true;
        BOOST_LOG_TRIVIAL(debug) << "End of files reached";
    }
}

void MultiFileReader::open_previous()
{
    _current_file_idx--;
    if(_current_file_idx >= 0) {
        _current_stream.close();
        open();
    } else {
        _eof_flag = true;
    }
}

std::size_t MultiFileReader::get_total_size() const
{
    return _total_size;
}

bool MultiFileReader::can_read(std::size_t bytes) const
{
    std::size_t remainingSize = _total_size - _current_position;
    return remainingSize >= bytes;
}

void MultiFileReader::seekg(long pos, std::ios_base::seekdir dir)
{
    BOOST_LOG_TRIVIAL(debug) << "Seeking to " << pos << " bytes";
    std::size_t targetPos = (dir == std::ios_base::cur)
                                ? _current_position + pos
                                : pos;       // current : begin;
    if(dir == std::ios_base::end && pos < 0) // end
        targetPos = this->_total_size + pos;

    if(targetPos > _total_size) {
        throw std::runtime_error("Invalid seek position: " +
                                 std::to_string(targetPos));
    }

    std::size_t cumulativeSize = 0;
    for(size_t i = 0; i < _dada_files.size(); i++) {
        cumulativeSize += _sizes[i];
        if(targetPos < cumulativeSize) {
            std::size_t filePos =
                targetPos - cumulativeSize + _sizes[i]; // position in file
            if(_current_file_idx != i) {
                _current_stream.close();
                _current_file_idx = i;
                open();
            }

            _current_stream.seekg(_dada_header_size + filePos, std::ios::beg);
            BOOST_LOG_TRIVIAL(debug)
                << "Actually seeking to " << filePos << " bytes from "
                << _dada_files[_current_file_idx];
            BOOST_LOG_TRIVIAL(debug)
                << "Current position: " << _current_stream.tellg();
            _current_position = targetPos;
            _eof_flag         = false;
            break;
        } else {
            BOOST_LOG_TRIVIAL(debug)
                << "Seek is after file: " << _dada_files[i];
        }
    }
    if(_current_position >= _total_size || targetPos >= _total_size) {
        _eof_flag = true;
    }
}

void MultiFileReader::close()
{
    if(_current_stream.is_open()) {
        _current_stream.close();
    }
    _eof_flag = true;
}

std::streamsize MultiFileReader::read(char* raw_ptr, std::streamsize bytes)
{
    std::streamsize totalRead = 0;

    while(bytes > 0) {
        BOOST_LOG_TRIVIAL(debug)
            << "Attempting to read " << bytes << " bytes starting from "
            << _current_stream.tellg() << " of "
            << _dada_files[_current_file_idx];

        if(!_current_stream.is_open()) {
            std::stringstream ss;
            ss << "File not open: " << _dada_files[_current_file_idx];
            throw std::runtime_error(ss.str());
        }

        _current_stream.read(raw_ptr + totalRead, bytes);

        std::streamsize bytesRead = _current_stream.gcount();
        BOOST_LOG_TRIVIAL(debug) << "Read " << bytesRead << " bytes from "
                                 << _dada_files[_current_file_idx];

        if(bytesRead < bytes) {
            BOOST_LOG_TRIVIAL(debug)
                << "This is less than required, opening next file"
                << _dada_files[_current_file_idx];
            open_next();
        }

        bytes -= bytesRead;
        totalRead += bytesRead;
        _current_position += bytesRead;

        if(_eof_flag)
            break;
    }
    if(_current_position >= _total_size) {
        _eof_flag = true;
    }
    BOOST_LOG_TRIVIAL(debug) << "Total read: " << totalRead;
    BOOST_LOG_TRIVIAL(debug) << "Current position: " << _current_position;
    return totalRead;
}

bool MultiFileReader::eof() const
{
    return _eof_flag;
}

bool MultiFileReader::good() const
{
    return _current_stream.good() && !_eof_flag && !_current_stream.bad() &&
           !_current_stream.fail();
}

std::size_t MultiFileReader::tellg() const
{
    return _current_position;
}

bool MultiFileReader::is_open() const
{
    return _current_stream.is_open();
}

template <typename T>
MultiFileReader& operator>>(MultiFileReader& reader, T& value)
{
    if(!(reader._current_stream >> value)) {
        reader._eof_flag = true;
    }
    return reader;
}
