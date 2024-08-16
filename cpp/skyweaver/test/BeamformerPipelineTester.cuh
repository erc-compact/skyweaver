#ifndef SKYWEAVER_TEST_BEAMFORMERPIPELINETESTER_CUH
#define SKYWEAVER_TEST_BEAMFORMERPIPELINETESTER_CUH

#include "skyweaver/BeamformerPipeline.cuh"
#include "skyweaver/PipelineConfig.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace skyweaver
{
namespace test
{

template <typename BfTraits>
class BeamformerPipelineTester: public ::testing::Test
{
  
  protected:
    void SetUp() override;
    void TearDown() override;

  public:
    using BfTraitsType = BfTraits;
    BeamformerPipelineTester();
    ~BeamformerPipelineTester();
  protected:
    PipelineConfig _config;
};

} // namespace test
} // namespace skyweaver

#endif // SKYWEAVER_TEST_BEAMFORMERPIPELINETESTER_CUH
