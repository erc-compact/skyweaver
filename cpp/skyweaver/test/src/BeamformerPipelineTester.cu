#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/IncoherentDedispersionPipeline.cuh"
#include "skyweaver/MultiFileWriter.cuh"
#include "skyweaver/PipelineConfig.hpp"
#include "skyweaver/beamformer_utils.cuh"
#include "skyweaver/test/BeamformerPipelineTester.cuh"

#include <filesystem>
namespace fs = std::filesystem;

namespace skyweaver
{
namespace test
{
namespace
{
std::string default_dada_header = R"(
HEADER       DADA
HDR_VERSION  1.0
HDR_SIZE     4096
DADA_VERSION 1.0

# DADA parameters
OBS_ID       3

FILE_SIZE    1073745920
FILE_NUMBER  0

# time of the rising edge of the first time sample
# UTC_START    1708082482.092000
UTC_START    1708082170.0
MJD_START    60356.473172

OBS_OFFSET   0
OBS_OVERLAP  0

SOURCE       J1644-4559_Offset1
RA           16:44:26.25
DEC          -45:59:09.6
TELESCOPE    MeerKAT
INSTRUMENT   CBF-Feng
RECEIVER     L-band
FREQ         1430414062.5
BW           13515625.0
OBS_FREQ     1284000000.000000
OBS_BW       856000000.0
TSAMP        0.0000047850467290

BYTES_PER_SECOND 3424000000.0

NBIT         8
NDIM         2
NPOL         2
NCHAN        64
OBS_NCHAN    4096
NANT         57
ORDER        TAFTP
INNER_T      256


#MeerKAT specifics
SYNC_TIME    1708039531.000000
SAMPLE_CLOCK 1712000000.0
SAMPLE_CLOCK_START 0.0
CHAN0_IDX 2688
)";
}

class NullHandler
{
  public:
    template <typename... Args>
    void init(Args... args){};

    template <typename... Args>
    bool operator()(Args... args)
    {
        return false;
    };
};

template <typename BfTraits>
BeamformerPipelineTester<BfTraits>::BeamformerPipelineTester()
    : ::testing::Test()
{
}

template <typename BfTraits>
BeamformerPipelineTester<BfTraits>::~BeamformerPipelineTester()
{
}

template <typename BfTraits>
void BeamformerPipelineTester<BfTraits>::SetUp()
{
    if(_config.nantennas() < 57) {
        GTEST_SKIP();
    }
    if(_config.nbeams() < 67) {
        GTEST_SKIP();
    }
}

template <typename BfTraits>
void BeamformerPipelineTester<BfTraits>::TearDown()
{
}

typedef ::testing::Types<SingleStokesBeamformerTraits<StokesParameter::I>,
                         SingleStokesBeamformerTraits<StokesParameter::Q>,
                         SingleStokesBeamformerTraits<StokesParameter::U>,
                         SingleStokesBeamformerTraits<StokesParameter::V>,
                         StokesTraits<StokesParameter::Q, StokesParameter::U>,
                         StokesTraits<StokesParameter::I, StokesParameter::V>,
                         FullStokesBeamformerTraits>
    StokesTypes;
TYPED_TEST_SUITE(BeamformerPipelineTester, StokesTypes);

TYPED_TEST(BeamformerPipelineTester, instantiate)
{
    using BfTraits = typename TestFixture::BfTraitsType;
    this->_config.ddplan().add_block(0.0f, 1);
    this->_config.output_dir("/tmp/");
    this->_config.delay_file("data/test_delays.bin");
    NullHandler cb_handler;
    NullHandler ib_handler;
    NullHandler stats_handler;
    BeamformerPipeline<decltype(cb_handler),
                       decltype(ib_handler),
                       decltype(stats_handler),
                       BfTraits>(this->_config,
                                 cb_handler,
                                 ib_handler,
                                 stats_handler);
}

TYPED_TEST(BeamformerPipelineTester, full_pipeline_test)
{
    using BfTraits = typename TestFixture::BfTraitsType;
    this->_config.ddplan().add_block(0.0f, 1);
    this->_config.output_dir("/tmp/");
    this->_config.delay_file("data/test_delays.bin");
    ObservationHeader header;
    std::vector<char> header_bytes(default_dada_header.size() + 1);
    std::strcpy(header_bytes.data(), default_dada_header.c_str());
    psrdada_cpp::RawBytes raw_header(header_bytes.data(),
                                     default_dada_header.size(),
                                     default_dada_header.size(),
                                     false);
    read_dada_header(raw_header, header);
    validate_header(header, this->_config);
    update_config(this->_config, header);
    using WriterType =
        MultiFileWriter<TDBPowersH<typename BfTraits::QuantisedPowerType>>;
    typename WriterType::CreateStreamCallBackType create_stream_callback =
        detail::create_dada_file_stream<
            TDBPowersH<typename BfTraits::QuantisedPowerType>>;
    WriterType cb_handler(this->_config, "cb", create_stream_callback);
    NullHandler ib_handler;
    NullHandler stats_handler;
    using IDPipelineType =
        IncoherentDedispersionPipeline<typename BfTraits::QuantisedPowerType,
                                       typename BfTraits::QuantisedPowerType,
                                       decltype(cb_handler)>;
    using BPipelineType    = BeamformerPipeline<IDPipelineType,
                                             decltype(ib_handler),
                                             decltype(stats_handler),
                                             BfTraits>;
    using InputVectorTypeH = typename BPipelineType::VoltageVectorTypeH;

    IDPipelineType dedispersion_pipeline(this->_config, cb_handler);
    BPipelineType pipeline(this->_config,
                           dedispersion_pipeline,
                           ib_handler,
                           stats_handler);

    InputVectorTypeH input({
        this->_config.gulp_length_samps() /
            this->_config.nsamples_per_heap(), // T
        header.nantennas,                      // A
        this->_config.nchans(),                // F
        this->_config.nsamples_per_heap(),     // T
        this->_config.npol()                   // P
    });
    input.frequencies(this->_config.channel_frequencies());
    input.dms({0.0f});
    pipeline.init(header);
    for(int ii = 0; ii < 100; ++ii) { pipeline(input); }
}

} // namespace test
} // namespace skyweaver
