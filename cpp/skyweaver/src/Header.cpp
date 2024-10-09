#include "skyweaver/Header.hpp"

#include "ascii_header.h"
#include "inttypes.h"

#include <cstring>
#include <iomanip>

namespace skyweaver
{

Header::Header(psrdada_cpp::RawBytes& header): _header(header)
{
}

Header::~Header()
{
}

void Header::purge()
{
    std::memset(static_cast<void*>(_header.ptr()), 0, _header.total_bytes());
}

void Header::fetch_header_string(char const* key) const
{
    if(!has_key(key)) {
        throw std::runtime_error(std::string("Could not find ") + key +
                                 " key in DADA header.");
    }
    if(strcmp("unset", _buffer) == 0) {
        throw std::runtime_error(std::string("The header key ") + key +
                                 " was unset");
    }
}

bool Header::has_key(char const* key) const
{
    return ascii_header_get(_header.ptr(), key, "%s", _buffer) != -1;
}



template <>
long double Header::get<long double>(char const* key) const
{
    fetch_header_string(key);
    long double value = std::strtold(_buffer, NULL);
    BOOST_LOG_TRIVIAL(debug)
        << "Header get: " << key << " = " << std::setprecision(15) << value;
    return value;
}

template <>
std::size_t Header::get<std::size_t>(char const* key) const
{
    fetch_header_string(key);
    std::size_t value = std::strtoul(_buffer, NULL, 0);
    BOOST_LOG_TRIVIAL(debug) << "Header get: " << key << " = " << value;
    return value;
}

template <>
std::string Header::get<std::string>(char const* key) const
{
    fetch_header_string(key);
    std::string value = std::string(_buffer);
    BOOST_LOG_TRIVIAL(debug) << "Header get: " << key << " = " << value;
    return value;
}

template <>
void Header::set<long double>(char const* key, long double value)
{
    BOOST_LOG_TRIVIAL(debug)
        << "Header set: " << key << " = " << std::setprecision(15) << value;
    ascii_header_set(this->_header.ptr(), key, "%Lf", value);
}

template <>
void Header::set<std::size_t>(char const* key, std::size_t value)
{
    BOOST_LOG_TRIVIAL(debug) << "Header set: " << key << " = " << value;
    ascii_header_set(this->_header.ptr(), key, "%" PRId64, value);
}

template <>
void Header::set<std::string>(char const* key, std::string value)
{
    BOOST_LOG_TRIVIAL(debug) << "Header set: " << key << " = " << value;
    ascii_header_set(this->_header.ptr(), key, "%s", value.c_str());
}

} // namespace skyweaver
