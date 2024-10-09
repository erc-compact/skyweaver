#include "skyweaver/SigprocHeader.hpp"

#include <string>

namespace skyweaver
{

template <typename T>
void SigprocHeader::header_write(std::ostream& stream,
                                 std::string const& name,
                                 T val)
{
    header_write(stream, name);
    stream.write(reinterpret_cast<char*>(&val), sizeof(val));
}

} // namespace skyweaver
