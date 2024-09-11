#include "skyweaver/FileOutputStream.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

namespace skyweaver
{

namespace fs = std::filesystem;

void create_directories(const fs::path& path)
{
    // Check if the directory already exists
    if(!fs::exists(path)) {
        // Directory does not exist, attempt to create it
        if(!fs::create_directories(path)) {
            throw std::runtime_error("Failed to create directory: " +
                                     path.string());
        }
    } else if(!fs::is_directory(path)) {
        // Path exists but is not a directory
        throw std::runtime_error("Path exists but is not a directory: " +
                                 path.string());
    }
}

FileOutputStream::File::File(std::string const& fname, std::size_t bytes)
    : _full_path(fname), _bytes_requested(bytes), _bytes_written(0)
{
    _stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    _stream.open(_full_path, std::ofstream::out | std::ofstream::binary);
    if(_stream.is_open()) {
        BOOST_LOG_TRIVIAL(info) << "Opened output file " << _full_path;
    } else {
        std::stringstream error_message;
        error_message << "Could not open file " << _full_path;
        BOOST_LOG_TRIVIAL(error) << error_message.str();
        throw std::runtime_error(error_message.str());
    }
}

FileOutputStream::File::~File()
{
    if(_stream.is_open()) {
        BOOST_LOG_TRIVIAL(info) << "Closing file " << _full_path;
        _stream.close();
    }
}

std::size_t FileOutputStream::File::write(PipelineConfig const& config, char const* ptr, std::size_t bytes)
{
    BOOST_LOG_TRIVIAL(debug)
        << "Writing " << bytes << " bytes to " << _full_path;
    std::size_t bytes_remaining = _bytes_requested - _bytes_written;
    try {
        if(bytes > bytes_remaining) {
            wait_for_space(bytes_remaining);
            _stream.write(ptr, bytes_remaining);
            _bytes_written += bytes_remaining;
            BOOST_LOG_TRIVIAL(debug)
                << "Partial write of " << bytes_remaining << " bytes";
            return bytes_remaining;
        } else {
            wait_for_space(bytes);
            _stream.write(ptr, bytes);
            _bytes_written += bytes;
            BOOST_LOG_TRIVIAL(debug) << "Completed write";
            return bytes;
        }
    } catch(const std::ofstream::failure& e) {
        std::string reason;

        if(_stream.bad()) {
            reason = "badbit set.";
        } else if(_stream.fail()) {
            reason = "failbit set.";
        } else if(_stream.eof()) {
            reason = "eofbit set.";
        }

        BOOST_LOG_TRIVIAL(error) << "Error while writing to " << _full_path
                                 << " (" << e.what() << ") because of reason: " << reason;
       
        throw;
    }
}

void FileOutputStream::File::wait_for_space(PipelineConfig const& config, size_t requested_bytes)
{
    std::filesystem::space_info space = std::filesystem::space(_full_path);
    if(space.available >= requested_bytes)
        return;

    BOOST_LOG_TRIVIAL(warning)
        << space.available
        << " bytes available space is not enough for writing "
        << requested_bytes << " bytes to " << _full_path << ".";
    BOOST_LOG_TRIVIAL(warning) << "Start pausing.";

    while(space.available < requested_bytes) {
        sleep(10);
        space = std::filesystem::space(_full_path);
    }
    BOOST_LOG_TRIVIAL(warning) << "Space has been freed up. Will proceed.";
}

FileOutputStream::FileOutputStream(PipelineConfig const& config,
                       std::string const& directory,
                       std::string const& base_filename,
                       std::string const& extension,
                       std::size_t bytes_per_file,
                       HeaderUpdateCallback header_update_callback)
    : _directory(directory), _base_filename(base_filename),
      _extension(extension), _bytes_per_file(bytes_per_file),
      _header_update_callback(header_update_callback), _total_bytes_written(0),
      _current_file(nullptr), _file_count(0), _config(config)
{
    if(_bytes_per_file == 0) {
        throw std::runtime_error(
            "The number of bytes per file must be greater than zero");
    }
    BOOST_LOG_TRIVIAL(debug) << "Creating output file stream with parameters,\n"
                             << "Output directory: " << _directory << "\n"
                             << "Base filename: " << _base_filename << "\n"
                             << "Extension: " << _extension << "\n"
                             << "Number of bytes per file: " << _bytes_per_file;
    // Can make this a cli arg if needed
    umask(0022);
    create_directories(directory);
}

FileOutputStream::~FileOutputStream()
{
    if(_current_file) {
        _current_file.reset(nullptr);
    }
}

void FileOutputStream::write(char const* ptr, std::size_t bytes)
{
    BOOST_LOG_TRIVIAL(debug) << "Writing " << bytes << " bytes to file stream";
    if(_current_file) {
        std::size_t bytes_written = _current_file->write(ptr, bytes);
        _total_bytes_written += bytes_written;
        if(bytes_written < bytes) {
            new_file();
            write(ptr + bytes_written, bytes - bytes_written);
        }
    } else {
        new_file();
        write(ptr, bytes);
    }
}

void FileOutputStream::new_file()
{
    std::stringstream full_path;
    full_path << _directory << "/" << _base_filename << "_" << std::setfill('0')
              << std::setw(16) << _total_bytes_written << std::setfill(' ')
              << _extension;
    std::size_t header_bytes;
    BOOST_LOG_TRIVIAL(debug) << "Retrieving updated header";
    // The callback needs to guarantee the lifetime of the returned pointer here
    std::shared_ptr<char const> header_ptr =
        _header_update_callback(header_bytes,
                                _total_bytes_written,
                                _file_count);
    _current_file.reset(
        new File(full_path.str(), _bytes_per_file + header_bytes));
    // Here we are directly invoking the write method on the File object
    // to avoid potential bugs when the header is not completely written
    BOOST_LOG_TRIVIAL(debug) << "Writing updated header";
    if(_current_file->write(header_ptr.get(), header_bytes) != header_bytes) {
        throw std::runtime_error("Unable to write header to File instance");
    }
    ++_file_count;
}

} // namespace skyweaver
