#include "skyweaver/SigprocHeader.hpp"

#include <boost/algorithm/string.hpp>
#include <chrono>
#include <iterator>
#include <string>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>




namespace skyweaver
{

SigprocHeader::SigprocHeader()
{
}

SigprocHeader::SigprocHeader(ObservationHeader const& obs_header)
{
    source = obs_header.source_name;
    auto ra_val = obs_header.ra;
    auto dec_val = obs_header.dec;
   	std::vector<std::string> ra_s;
	std::vector<std::string> dec_s;      
    boost::split(ra_s, ra_val, boost::is_any_of(":"));
    boost::split(dec_s, dec_val, boost::is_any_of(":"));
    ra = stod(boost::join(ra_s, ""));
    dec = stod(boost::join(dec_s, ""));
    fch1 = obs_header.fch1;
    foff = obs_header.foff;
    tsamp = obs_header.tsamp;
    tstart = obs_header.mjd_start;
    nbits = obs_header.nbits;
    nifs = obs_header.nifs;
    nchans = obs_header.nchans;
    ibeam = obs_header.ibeam;
    rawfile = obs_header.rawfile;
    telescopeid = obs_header.telescopeid;
    machineid = obs_header.machineid;
    datatype = obs_header.datatype;
    barycentric = obs_header.barycentric;
    nbeams = obs_header.nbeams;
    az = obs_header.az;
    za = obs_header.za;
}


void SigprocHeader::add_time_offset(double offset_mjd)
{
    tstart += offset_mjd;
}

SigprocHeader::~SigprocHeader()
{
}

void SigprocHeader::header_write(std::ostream& stream, std::string const& str)
{
    int len = str.size();
    stream.write(reinterpret_cast<char*>(&len), sizeof(int));
    stream << str;
}

void SigprocHeader::header_write(std::ostream& stream,
                                 std::string const& str,
                                 std::string const& val)
{
    header_write(stream, str);
    header_write(stream, val);
}

void SigprocHeader::write_header(std::ostream& stream)
{
    header_write(stream, "HEADER_START");
    header_write<std::uint32_t>(stream, "telescope_id",telescopeid);
    header_write<std::uint32_t>(stream, "machine_id",machineid);
    header_write<std::uint32_t>(stream, "data_type",datatype);
    header_write<std::uint32_t>(stream, "barycentric",barycentric);
    header_write(stream, "source_name",source);
    header_write<double>(stream, "src_raj",ra);
    header_write<double>(stream, "src_dej",dec);
    header_write<std::uint32_t>(stream, "nbits",nbits);
    header_write<std::uint32_t>(stream, "nifs",nifs);
    header_write<std::uint32_t>(stream, "nchans",nchans);
    header_write<std::uint32_t>(stream, "ibeam",ibeam);
    header_write<double>(stream, "fch1",fch1);
    header_write<double>(stream, "foff",foff);
    header_write<double>(stream, "tstart",tstart);
    header_write<double>(stream, "tsamp",tsamp);
    header_write(stream, "HEADER_END");
}

} // namespace skyweaver
