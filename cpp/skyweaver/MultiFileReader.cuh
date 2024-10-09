#ifndef MULTIFILE_READER_HPP
#define MULTIFILE_READER_HPP

#include "cuda.h"
#include "psrdada_cpp/psrdadaheader.hpp"
#include "psrdada_cpp/raw_bytes.hpp"
#include "skyweaver/ObservationHeader.hpp"
#include "skyweaver/PipelineConfig.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>
namespace skyweaver
{
class MultiFileReader;
}

//TODO: Document this interface

class skyweaver::MultiFileReader
{
  private:
    std::vector<std::string> _dada_files;
    std::vector<std::size_t> _sizes;
    std::vector<ObservationHeader> _headers;
    std::ifstream _current_stream;
    int _current_file_idx;
    std::size_t _current_position;
    bool _eof_flag;
    std::size_t _dada_header_size;
    std::size_t _total_size;
    bool _is_open;
    

    void read_header(std::vector<char>& headerBytes);
    void check_contiguity();
    std::shared_ptr<ObservationHeader> header;

  public:
    MultiFileReader(PipelineConfig const& config);
    MultiFileReader(std::vector<std::string> dada_files,std::size_t dada_header_size, bool check_contiguity = false);
    void open();
    void open_next();
    void open_previous();
    void close();
    void seekg(long pos, std::ios_base::seekdir dir = std::ios_base::beg);
    
    std::size_t tellg() const;
    bool eof() const;
    bool good() const;
    bool can_read(std::size_t bytes) const;
    std::streamsize read(char* raw_ptr, std::streamsize bytes);
    

    bool is_open() const;
    std::size_t
    safe_read(std::ifstream& input_stream, char* buffer, std::size_t nbytes);
    std::vector<std::string> const& get_dada_files() const;

    std::size_t get_total_size() const;

    skyweaver::ObservationHeader const& get_header() const;
    skyweaver::ObservationHeader const& get_current_header() const;

    template <typename T>
    friend MultiFileReader& operator>>(MultiFileReader& reader, T& value);

    ~MultiFileReader();
};

#endif // MULTIFILE_READER_HPP
