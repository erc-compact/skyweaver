#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/test/DelayManagerTester.cuh"
#include "thrust/host_vector.h"

namespace skyweaver
{
namespace test
{

DelayManagerTester::DelayManagerTester(): ::testing::Test(), _stream(0)
{
}

DelayManagerTester::~DelayManagerTester()
{
}

void DelayManagerTester::SetUp()
{
  CUDA_ERROR_CHECK(cudaStreamCreate(&_stream));
}

void DelayManagerTester::TearDown()
{
  CUDA_ERROR_CHECK(cudaStreamDestroy(_stream));
}

TEST_F(DelayManagerTester, test_valid_read)
{
  // TODO
}

TEST_F(DelayManagerTester, test_invalid_read)
{
  // TODO
}

} // namespace test
} // namespace skyweaver