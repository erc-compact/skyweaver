#ifndef MULTIFILE_READER_HPP
#define MULTIFILE_READER_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mutex>
#include <memory>
#include "PsrDadaHeader.hpp"


class MultiFileReader {
private:
    std::vector<std::string> files;
    std::vector<long> sizes;
    std::ifstream currentStream;
    int currentFileIndex;
    std::size_t currentPosition;
    std::mutex mtx;
    bool eofFlag;
    std::size_t headerSize;
    std::size_t totalSize;

    void read_header();
    std::unique_ptr<psrdada_cpp::PsrDadaHeader> header;

public:
    MultiFileReader(const std::vector<std::string>& fileNames, long headerSize = 0);
    void open();
    void open_next();
    void open_previous();
    void close();
    void seekg(long pos, std::ios_base::seekdir dir = std::ios_base::beg);
    long tellg() const;
    bool eof() const;
    bool good() const;
    std::streamsize read(std::unique_ptr<std::vector<char>> buffer, std::streamsize bytes);

    template <typename T>
    friend MultiFileReader& operator>>(MultiFileReader& reader, T& value);
};

#endif // MULTIFILE_READER_HPP