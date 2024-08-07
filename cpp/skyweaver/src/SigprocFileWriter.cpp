#include "skyweaver/SigprocFileWriter.hpp"

#include <ctime>
#include <iomanip>

namespace skyweaver
{

SigprocFileWriter::SigprocFileWriter()
    : _max_filesize(1 << 20), _state(DISABLED), _tag(""), _outdir("./"),
      _extension(".fil"), _total_bytes(0), _new_stream_required(false),
      _stream(nullptr)
{
}

SigprocFileWriter::~SigprocFileWriter()
{
    _stream.reset(nullptr);
}

void SigprocFileWriter::tag(std::string const& tag_)
{
    _tag = tag_;
}

std::string const& SigprocFileWriter::tag() const
{
    return _tag;
}

void SigprocFileWriter::directory(std::string const& dir)
{
    _outdir = dir;
}

std::string const& SigprocFileWriter::directory() const
{
    return _outdir;
}

void SigprocFileWriter::max_filesize(std::size_t size)
{
    _max_filesize = size;
}

std::size_t SigprocFileWriter::max_filesize() const
{
    return _max_filesize;
}

void SigprocFileWriter::init(psrdada_cpp::RawBytes& _header)
{
    SigprocHeader parser;
    parser.read_header(block, _header);
}

bool SigprocFileWriter::operator()(psrdada_cpp::RawBytes& block)
{
    if(_state == ENABLED) {
        BOOST_LOG_TRIVIAL(debug) << "(" << _tag << ") Writer state is ENABLED";
        if((_stream == nullptr) || _new_stream_required) {
            new_stream();
        }
        // Also debugging print here the block index
        BOOST_LOG_TRIVIAL(debug) << "Writing DADA block at: " << block.ptr();
        _stream->write(block.ptr(), block.used_bytes());
    } else if(_state == DISABLED) {
        BOOST_LOG_TRIVIAL(debug) << "(" << _tag << ") Writer state is DISABLED";
        if(_stream != nullptr) {
            _stream.reset(nullptr);
        }
    } else {
        BOOST_LOG_TRIVIAL(error)
            << "(" << _tag << ") Writer state is neither ENABLED or DISABLED";
    }
    _total_bytes += block.used_bytes();
    // For debugging purposes we memset the buffer to zero here to
    // avoid "ghost" pulses
    std::memset(block.ptr(), 0, block.used_bytes());
    return false;
}

void SigprocFileWriter::new_stream()
{
    BOOST_LOG_TRIVIAL(info)
        << "Header parameters for new file with tag '" << _tag << "':";
    BOOST_LOG_TRIVIAL(info) << "rawfile: " << _header.rawfile;
    BOOST_LOG_TRIVIAL(info) << "source: " << _header.source;
    BOOST_LOG_TRIVIAL(info) << "az: " << _header.az;
    BOOST_LOG_TRIVIAL(info) << "dec: " << _header.dec;
    BOOST_LOG_TRIVIAL(info) << "fch1: " << _header.fch1;
    BOOST_LOG_TRIVIAL(info) << "foff: " << _header.foff;
    BOOST_LOG_TRIVIAL(info) << "ra: " << _header.ra;
    BOOST_LOG_TRIVIAL(info) << "rdm: " << _header.rdm;
    BOOST_LOG_TRIVIAL(info) << "tsamp: " << _header.tsamp;
    BOOST_LOG_TRIVIAL(info) << "tstart: " << _header.tstart;
    BOOST_LOG_TRIVIAL(info) << "za: " << _header.za;
    BOOST_LOG_TRIVIAL(info) << "datatype: " << _header.datatype;
    BOOST_LOG_TRIVIAL(info) << "barycentric: " << _header.barycentric;
    BOOST_LOG_TRIVIAL(info) << "ibeam: " << _header.ibeam;
    BOOST_LOG_TRIVIAL(info) << "machineid: " << _header.machineid;
    BOOST_LOG_TRIVIAL(info) << "nbeams: " << _header.nbeams;
    BOOST_LOG_TRIVIAL(info) << "nbits: " << _header.nbits;
    BOOST_LOG_TRIVIAL(info) << "nchans: " << _header.nchans;
    BOOST_LOG_TRIVIAL(info) << "nifs: " << _header.nifs;
    BOOST_LOG_TRIVIAL(info) << "telescopeid: " << _header.telescopeid;
    // Here we should update the tstart of the default header to be the
    // start of the stream
    _header.tstart =
        _header.tstart +
        (((_total_bytes / (_header.nbits / 8.0)) / (_header.nchans)) *
         _header.tsamp) /
            (86400.0);
    // reset the total bytes counter to keep the time tracked correctly
    _total_bytes = 0;
    // Generate the new base filename in <utc>_<tag> format
    std::stringstream base_filename;
    // Get UTC time string
    std::time_t unix_time =
        static_cast<std::time_t>((_header.tstart - 40587.0) * 86400.0);
    struct std::tm* ptm = std::gmtime(&unix_time);

    // Swapped out put_time call for strftime due to put_time
    // causing compiler bugs prior to g++ 5.x
    char formatted_time[80];
    strftime(formatted_time, 80, "%Y-%m-%d-%H:%M:%S", ptm);
    base_filename << formatted_time;

    if(_tag != "") {
        base_filename << "_" << _tag;
    }

    _stream.reset(new FileStream(
        _outdir,
        base_filename.str(),
        _extension,
        _max_filesize,
        [&](std::size_t& header_size,
            std::size_t bytes_written,
            std::size_t file_count) -> std::shared_ptr<char const> {
            // We do not explicitly delete[] this array
            // Cleanup is handled by the shared pointer
            // created below
            char* temp_header = new char[4096];
            SigprocHeader parser;
            // Make a copy of the header to edit
            FilHead fh = header();
            // Here we do any required updates to the header
            fh.tstart = fh.tstart +
                        (((bytes_written / (fh.nbits / 8.0)) / (fh.nchans)) *
                         fh.tsamp) /
                            (86400.0);
            header_size = parser.write_header(temp_header, fh);
            std::shared_ptr<char const> header_ptr(
                temp_header,
                std::default_delete<char[]>());
            return header_ptr;
        }));
    _new_stream_required = false;
}

void SigprocFileWriter::enable()
{
    if(_state == DISABLED) {
        _new_stream_required = true;
    }
    _state = ENABLED;
}

void SigprocFileWriter::disable()
{
    _state               = DISABLED;
    _new_stream_required = false;
}

FilHead& SigprocFileWriter::header()
{
    return _header;
}

} // namespace skyweaver
