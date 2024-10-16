#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/skyweaver_constants.hpp"
#include "skyweaver/test/WeightsManagerTester.cuh"
#include "thrust/host_vector.h"

#define TWOPI 6.283185307179586

namespace skyweaver
{
namespace test
{

WeightsManagerTester::WeightsManagerTester(): ::testing::Test(), _stream(0)
{
}

WeightsManagerTester::~WeightsManagerTester()
{
}

void WeightsManagerTester::SetUp()
{
    _config.centre_frequency(1.4e9);
    _config.bandwidth(56.0e6);
    CUDA_ERROR_CHECK(cudaStreamCreate(&_stream));
}

void WeightsManagerTester::TearDown()
{
    CUDA_ERROR_CHECK(cudaStreamDestroy(_stream));
}

void WeightsManagerTester::calc_weights_c_reference(
    thrust::host_vector<float3> const& delay_models,
    thrust::host_vector<char2>& weights,
    std::vector<double> const& channel_frequencies,
    int nantennas,
    int nbeams,
    int nchans,
    double current_epoch,
    double delay_epoch,
    double tstep,
    int ntsteps)
{
    double2 weight;
    char2 compressed_weight;
    for(int antenna_idx = 0; antenna_idx < nantennas; ++antenna_idx) {
        for(int beam_idx = 0; beam_idx < nbeams; ++beam_idx) {
            float3 delay_model =
                delay_models[beam_idx * nantennas + antenna_idx];
            for(int chan_idx = 0; chan_idx < nchans; ++chan_idx) {
                double frequency = channel_frequencies[chan_idx];
                for(int time_idx = 0; time_idx < ntsteps; ++time_idx) {
                    double t = (current_epoch - delay_epoch) + time_idx * tstep;
                    double phase =
                        (t * delay_model.z + delay_model.y) * frequency;
                    sincos(TWOPI * phase, &weight.y, &weight.x);
                    compressed_weight.x =
                        (char)round(weight.x * 127.0 * delay_model.x);
                    compressed_weight.y =
                        (char)round(-1.0 * weight.y * 127.0 * delay_model.x);
                    int output_idx =
                        nantennas * (nbeams * (time_idx * nchans + chan_idx) +
                                     beam_idx) +
                        antenna_idx;
                    weights[output_idx] = compressed_weight;
                }
            }
        }
    }
}

void WeightsManagerTester::compare_against_host(
    DelayVectorTypeD const& delays,
    WeightsVectorTypeD const& weights,
    TimeType current_epoch,
    TimeType delay_epoch)
{
    // Implicit device to host copies
    thrust::host_vector<float3> host_delays = delays;
    thrust::host_vector<char2> cuda_weights = weights;
    thrust::host_vector<char2> c_weights(cuda_weights.size());
    calc_weights_c_reference(host_delays,
                             c_weights,
                             _config.channel_frequencies(),
                             _config.nantennas(),
                             _config.nbeams(),
                             _config.channel_frequencies().size(),
                             current_epoch,
                             delay_epoch,
                             0.0,
                             1);
    for(int ii = 0; ii < cuda_weights.size(); ++ii) {
        ASSERT_EQ(c_weights[ii].x, cuda_weights[ii].x);
        ASSERT_EQ(c_weights[ii].y, cuda_weights[ii].y);
    }
}

TEST_F(WeightsManagerTester, test_zero_value)
{
    // This is always the size of the delay array
    std::size_t delays_size = _config.nbeams() * _config.nantennas();
    WeightsManager weights_manager(_config, _stream);
    // This is a thrust::device_vector<float3>
    DelayVectorTypeD delays(delays_size, {1.0f, 0.0f, 0.0f});
    TimeType current_epoch = 10.0;
    TimeType delay_epoch   = 9.0;
    // First try everything with only zeros
    auto const& weights =
        weights_manager.weights(delays, current_epoch, delay_epoch);
    compare_against_host(delays, weights, current_epoch, delay_epoch);
}

TEST_F(WeightsManagerTester, test_real_value)
{
    // This is always the size of the delay array
    std::size_t delays_size = _config.nbeams() * _config.nantennas();
    WeightsManager weights_manager(_config, _stream);
    // This is a thrust::device_vector<float3>
    DelayVectorTypeD delays(delays_size, {1.0f, 1e-11f, 1e-10f});
    TimeType current_epoch = 10.0;
    TimeType delay_epoch   = 9.0;
    // First try everything with only zeros
    auto const& weights =
        weights_manager.weights(delays, current_epoch, delay_epoch);
    compare_against_host(delays, weights, current_epoch, delay_epoch);
}

TEST_F(WeightsManagerTester, test_real_values)
{
    // This is always the size of the delay array
    std::size_t delays_size = _config.nbeams() * _config.nantennas();
    WeightsManager weights_manager(_config, _stream);
    // This is a thrust::device_vector<float3>
    DelayVectorTypeD delays(delays_size, {0.0f, 0.0f});
    for(int ii = 0; ii < delays_size; ++ii) {
        delays[ii] = {ii * 1e-11f, ii * 1e-15f};
    }
    TimeType current_epoch = 10.0;
    TimeType delay_epoch   = 9.0;
    // First try everything with only zeros
    auto const& weights =
        weights_manager.weights(delays, current_epoch, delay_epoch);
    compare_against_host(delays, weights, current_epoch, delay_epoch);
}

} // namespace test
} // namespace skyweaver
