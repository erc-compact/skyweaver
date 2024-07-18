#ifndef SKYWEAVER_TEST_DELAYMANAGERTESTER_CUH
#define SKYWEAVER_TEST_DELAYMANAGERTESTER_CUH

#include "skyweaver/DelayManager.cuh"
#include "skyweaver/PipelineConfig.hpp"

#include <gtest/gtest.h>

namespace skyweaver
{
namespace test
{

class DelayManagerTester: public ::testing::Test
{
  protected:
    void SetUp() override;
    void TearDown() override;

  public:
    DelayManagerTester();
    ~DelayManagerTester();

  protected:
    PipelineConfig _config;
    cudaStream_t _stream;
};

} // namespace test
} // namespace skyweaver

#endif // SKYWEAVER_TEST_DELAYMANAGERTESTER_CUH