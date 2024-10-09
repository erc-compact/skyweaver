#ifndef SKYWEAVER_SKYCLEAVER_HPP
#define SKYWEAVER_SKYCLEAVER_HPP
#include "boost/log/trivial.hpp"
#include "psrdada_cpp/psrdadaheader.hpp"
#include "psrdada_cpp/raw_bytes.hpp"
#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/MultiFileReader.cuh"
#include "skyweaver/MultiFileWriter.cuh"
#include "skyweaver/ObservationHeader.hpp"
#include "skyweaver/SkyCleaverConfig.hpp"
#include "skyweaver/Timer.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace skyweaver
{

struct BridgeReader {
  public:
    std::vector<std::string> _tdb_filenames;
    std::unique_ptr<MultiFileReader> _tdb_reader;
    std::string freq;

}; // BridgeReader


class SkyCleaver
{
  public:
    using InputType        = int8_t;
    using OutputType       = uint8_t;
    using InputVectorType  = TDBPowersH<InputType>;
    using OutputVectorType = TFPowersH<OutputType>;
    using FreqType         = std::size_t; // up to the nearest Hz
    using BeamNumberType   = std::size_t;
    using DMNumberType     = std::size_t;

  private:
    SkyCleaverConfig const& _config;
    std::map<FreqType, std::unique_ptr<MultiFileReader>> _bridge_readers;
    std::map<FreqType, std::unique_ptr<InputVectorType>> _bridge_data;

    std::vector<FreqType> _expected_freqs;
    std::vector<FreqType> _available_freqs;
    std::size_t _nsamples_to_read;
    ObservationHeader _header;

    std::vector<std::string> _beam_filenames;
    std::map<DMNumberType,
             std::map<BeamNumberType,
                      std::unique_ptr<MultiFileWriter<OutputVectorType>>>>
        _beam_writers;
    std::map<DMNumberType,
             std::map<BeamNumberType, std::shared_ptr<OutputVectorType>>>
        _beam_data;
    std::size_t _total_beam_writers;

    void init_readers();
    void init_writers();

    Timer _timer;

  public:
    SkyCleaver(SkyCleaverConfig const& config);
    SkyCleaver(SkyCleaver const&)     = delete;
    void operator=(SkyCleaver const&) = delete;

    void cleave();

}; // class SkyCleaver

} // namespace skyweaver

#endif // SKYWEAVER_SKYCLEAVER_HPP