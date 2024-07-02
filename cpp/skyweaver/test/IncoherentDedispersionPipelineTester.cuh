#ifndef SKYWEAVER_TEST_INCOHERENTDEDISPERSIONPIPELINETESTER_CUH
#define SKYWEAVER_TEST_INCOHERENTDEDISPERSIONPIPELINETESTER_CUH

#include "skyweaver/PipelineConfig.hpp"
#include "skyweaver/ObservationHeader.hpp"
#include "thrust/host_vector.h"
#include <gtest/gtest.h>

namespace skyweaver
{
namespace test
{

template <typename InputType_, typename OutputType_>
struct IDPipelineTraits
{
    using InputType = InputType_;
    using OutputType = OutputType_;
};

template <typename Traits>
class IncoherentDedispersionPipelineTester: public ::testing::Test
{
  public:
    using OutputVectorType = thrust::host_vector<typename Traits::OutputType>;

  protected:
    void SetUp() override;
    void TearDown() override;

  public:
    IncoherentDedispersionPipelineTester();
    ~IncoherentDedispersionPipelineTester();

  void init(ObservationHeader const& header);

  void operator()(OutputVectorType const&, std::size_t dm_idx);

  protected:
    bool _init_called;
    int _operator_call_count;
    PipelineConfig _config;
    ObservationHeader _init_arg;
    OutputVectorType _operator_arg_1;
    std::size_t _operator_arg_2;
};

} // namespace test
} // namespace skyweaver

#endif // SKYWEAVER_TEST_INCOHERENTDEDISPERSIONPIPELINETESTER_CUH