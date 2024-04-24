#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/test/DelayManagerTester.cuh"
#include "thrust/host_vector.h"

#include <cstdlib>
#include <iostream>

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

/**
 * The data/ directory contains a test delay file called
 * test_delays.bin.
 *
 * It contains 5 delay models for 67 beams and 57 antennas.
 *
 * The first valid epoch is   1708082168.957
 * The last valid epoch is 1708082468.957
 *
 * delay offsets are in the range 1e-7 to 1e-10
 * delay rates are in the range 1e-12 to 1e-14
 * weights are either 1 or 0
 */

TEST_F(DelayManagerTester, test_valid_read)
{
    DelayManager delay_manager("data/test_delays.bin", _stream);
}

TEST_F(DelayManagerTester, test_invalid_read)
{
    // TODO
}

} // namespace test
} // namespace skyweaver