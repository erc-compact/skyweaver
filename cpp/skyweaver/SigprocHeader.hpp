
#ifndef SKYWEAVER_SIGPROCHEADER_HPP
#define SKYWEAVER_SIGPROCHEADER_HPP

#include "psrdada_cpp/raw_bytes.hpp"
#include "skyweaver/ObservationHeader.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace skyweaver
{


class SigprocHeader
{
  public:
    SigprocHeader();

    SigprocHeader(ObservationHeader const& obs_header);

    SigprocHeader(SigprocHeader const&) = delete; 

    ~SigprocHeader();

    void write_header(std::ostream& stream);
    void add_time_offset(double offset_mjd);

  private:
    std::string rawfile  = "unset";
    std::string source   = "unset";
    double az            = 0.0; // azimuth angle in deg
    double dec           = 0.0; // source declination
    double fch1          = 0.0; // frequency of the top channel in MHz
    double foff          = 0.0; // channel bandwidth in MHz
    double ra            = 0.0; // source right ascension
    double rdm           = 0.0; // reference DM
    double tsamp         = 0.0; // sampling time in seconds
    double tstart        = 0.0; // observation start time in MJD format
    double za            = 0.0; // zenith angle in deg
    uint32_t datatype    = 0;   // data type ID
    uint32_t barycentric = 0;   // barycentric flag
    uint32_t ibeam       = 0;   // beam number
    uint32_t machineid   = 0;   // machine ID 0 = FAKE
    uint32_t nbeams      = 0;   // number of beams in the observation
    uint32_t nbits       = 0;   // number of bits per sample
    uint32_t nchans      = 0;   // number of frequency channels
    uint32_t nifs        = 0;   // number of ifs (pols)
    uint32_t telescopeid = 0;   // telescope ID

    /*
     * @brief write string to the header
     */
    void header_write(std::ostream& stream, std::string const& str);

    void header_write(std::ostream& stream,
                      std::string const& str,
                      std::string const& name);


    /*
     * @brief write a value to the stream
     */
    template <typename T>
    void header_write(std::ostream& stream, std::string const& name, T val);
};

} // namespace skyweaver

#include "skyweaver/detail/SigprocHeader.cpp"

#endif // SKYWEAVER_SIGPROCHEADER_HPP
