#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/test/MultiFileWriterTester.cuh"
#include "skyweaver/MultiFileWriter.cuh"

#include <cstdio>
#include <filesystem>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace skyweaver
{
namespace test
{

MultiFileWriterTester::MultiFileWriterTester(): ::testing::Test()
{
}

MultiFileWriterTester::~MultiFileWriterTester()
{
}

void MultiFileWriterTester::SetUp()
{
    char template_dirname[] = "/tmp/skyweaver_test_XXXXXX";
    char* directory_path    = mkdtemp(template_dirname);
    _config.output_dir(std::string(directory_path));
    auto& plan = _config.ddplan();
    for(float dm = 0.0f; dm < 5; dm += 1.2345f) { plan.add_block(dm); }
}

void MultiFileWriterTester::TearDown()
{
    fs::remove_all(_config.output_dir());
}

TEST_F(MultiFileWriterTester, simple_updating_write)
{
    using InputType = TFBPowersH<int8_t>;

    BOOST_LOG_TRIVIAL(debug)
        << "Testing in tmp directory: " << _config.output_dir();

    using WriterType = MultiFileWriter<InputType>;
    typename WriterType::CreateStreamCallBackType create_stream_callback = detail::create_dada_file_stream<InputType>;
       

    WriterType mfw(_config, "test", create_stream_callback);
    _config.max_output_filesize(1000);
    InputType powers({_config.nsamples_per_block(), 64, _config.nbeams()});
    powers.dms({0.0});
    ObservationHeader header;
    header.nchans = _config.nchans();
    header.tsamp  = 0.000064;
    BOOST_LOG_TRIVIAL(debug) << "HEADER: \n" << header.to_string();
    mfw.init(header);
    for(std::size_t file_idx = 0; file_idx < 4; ++file_idx) {
        for(std::size_t dm_idx = 0; dm_idx < _config.coherent_dms().size();
            ++dm_idx) {
            mfw(powers, dm_idx);
        }
    }
}

} // namespace test
} // namespace skyweaver