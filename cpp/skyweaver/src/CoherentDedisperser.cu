#include "cuda.h"
#include "cufft.h"
#include "skyweaver/CoherentDedisperser.cuh"
#include "skyweaver/dedispersion_utils.cuh"

#include <cmath>
#include <cuda_runtime.h>
#include <iomanip>
#include <iostream>
#include <thrust/complex.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <vector>

namespace
{

#define NCHANS_PER_BLOCK 128

// Function to check if a number's prime factors are only 2, 3, 5, or
bool has_small_prime_factors(size_t n)
{
    if(n == 0)
        return false;
    while(n % 2 == 0) n /= 2;
    while(n % 3 == 0) n /= 3;
    while(n % 5 == 0) n /= 5;
    while(n % 7 == 0) n /= 7;
    return n == 1;
}

// Function to find the next optimal FFT size
size_t next_optimal_fft_size(size_t N)
{
    size_t n = N;
    while(!has_small_prime_factors(n)) { n++; }
    return n;
}

// Function to compute padding amount
size_t compute_padding(size_t N)
{
    size_t optimal_size = 0;
    size_t n            = N - 1;

    do {
        n            = n + 1;
        optimal_size = next_optimal_fft_size(n);
        BOOST_LOG_TRIVIAL(info) << "Trying optimal size: " << optimal_size;
    } while(optimal_size % NCHANS_PER_BLOCK != 0);
    return optimal_size - N;
}

} // namespace

namespace skyweaver
{

void create_coherent_dedisperser_config(CoherentDedisperserConfig& config,
                                        PipelineConfig const& pipeline_config)
{
    // the centre frequency and bandwidth are for the bridge. This is taken from
    // Observation Header (not from the user)
    float f_low =
        pipeline_config.centre_frequency() - pipeline_config.bandwidth() / 2.0f;

    float f_high =
        f_low + pipeline_config.bandwidth() / pipeline_config.nchans();

    // pipeline_config.centre_frequency() + pipeline_config.bandwidth() / 2.0f;
    float tsamp = pipeline_config.nchans() / pipeline_config.bandwidth();

    if(pipeline_config.coherent_dms().empty()) {
        throw std::runtime_error("No coherent DMs specified");
    }

    auto it      = std::max_element(pipeline_config.coherent_dms().begin(),
                               pipeline_config.coherent_dms().end());
    float max_dm = *it;
    BOOST_LOG_TRIVIAL(debug) << "Constructing coherent dedisperser plan";
    std::size_t max_dm_delay_samps =
        DMSampleDelay(max_dm, f_low, tsamp)(f_high);

    if(max_dm_delay_samps > 2 * pipeline_config.gulp_length_samps()) {
        throw std::runtime_error(
            "Gulp length must be at least 2 times the maximum DM delay");
    }

    if((pipeline_config.gulp_length_samps() + max_dm_delay_samps) % 2 != 0) {
        max_dm_delay_samps++;
    }

    create_coherent_dedisperser_config(config,
                                       pipeline_config.gulp_length_samps(),
                                       max_dm_delay_samps,
                                       pipeline_config.nchans(),
                                       pipeline_config.npol(),
                                       pipeline_config.nantennas(),
                                       tsamp,
                                       f_low,
                                       pipeline_config.bandwidth(),
                                       pipeline_config.coherent_dms());
}
/*
 * @brief      Create a new CoherentDedisperser object, mostly used only for
 * testing
 *
 * @param      config  The config reference
 */
void create_coherent_dedisperser_config(CoherentDedisperserConfig& config,
                                        std::size_t gulp_samps,
                                        std::size_t overlap_samps,
                                        std::size_t num_coarse_chans,
                                        std::size_t npols,
                                        std::size_t nantennas,
                                        double tsamp,
                                        double low_freq,
                                        double bw,
                                        std::vector<float> dms)
{
    config.gulp_samps       = gulp_samps;
    config.num_coarse_chans = num_coarse_chans;
    config.npols            = npols;
    config.nantennas        = nantennas;
    config.tsamp            = tsamp;
    config.low_freq         = low_freq;
    config.bw               = bw;
    config.high_freq        = low_freq + bw;
    config.coarse_chan_bw   = bw / num_coarse_chans;
    config.filter_delay     = tsamp * overlap_samps / 2.0;
    BOOST_LOG_TRIVIAL(debug)
        << "tsamp in create_coherent_dedisperser_config: " << config.tsamp;
    BOOST_LOG_TRIVIAL(info)
        << "overlap_samps just due to dedispersion: " << overlap_samps;
    BOOST_LOG_TRIVIAL(info) << "Filter delay: " << config.filter_delay;

    config.overlap_samps =
        compute_padding(gulp_samps + overlap_samps) + overlap_samps;

    BOOST_LOG_TRIVIAL(info)
        << "overlap_samps after padding for optimal FFT size "
        << config.overlap_samps;

    /* Precompute DM constants */
    config._h_dms = dms;
    config._d_dms = config._h_dms;
    config._d_dm_prefactor.resize(dms.size());
    config._d_ism_responses.resize(dms.size());
    for(int i = 0; i < dms.size(); i++) {
        config._d_ism_responses[i].resize(
            num_coarse_chans * (config.gulp_samps + config.overlap_samps));
    }

    thrust::transform(config._d_dms.begin(),
                      config._d_dms.end(),
                      config._d_dm_prefactor.begin(),
                      DMPrefactor());

    config.fine_chan_bw = config.coarse_chan_bw / config.gulp_samps;

    for(int idx = 0; idx < config._d_dms.size(); idx++) {
        get_dm_responses(config,
                         config._d_dm_prefactor[idx],
                         config._d_ism_responses[idx]);
    }

    // data is FTPA order, we will loop over F, so we are left with TPA order.
    // Let's fuse PA to X, so TX order.
    //  We stride and batch over X and transform T
    std::size_t X        = config.npols * config.nantennas;
    std::size_t fft_size = config.gulp_samps + config.overlap_samps;
    int n[1]             = {static_cast<int>(fft_size)}; // FFT size
    int inembed[1]       = {static_cast<int>(fft_size)};
    int onembed[1]       = {static_cast<int>(fft_size)};
    int istride          = X;
    int ostride          = X;
    int idist            = 1;
    int odist            = 1;
    int batch            = X;

    if(cufftPlanMany(&config._fft_plan,
                     1,
                     n,
                     inembed,
                     istride,
                     idist,
                     onembed,
                     ostride,
                     odist,
                     CUFFT_C2C,
                     batch) != CUFFT_SUCCESS) {
        std::runtime_error("CUFFT error: Plan creation failed");
    }

    BOOST_LOG_TRIVIAL(debug) << "FFT plan created";
}

void CoherentDedisperser::dedisperse(
    TPAVoltagesD<char2> const& d_tpa_voltages_in,
    FTPAVoltagesD<char2>& d_ftpa_voltages_out,
    unsigned int freq_idx,
    unsigned int dm_idx)
{
    BOOST_LOG_NAMED_SCOPE("CoherentDedisperser::dedisperse");
    // d_tpa_voltages_in.size() is with overlap
    _d_fpa_spectra.resize(d_tpa_voltages_in.size(), {0.0f, 0.0f});
    _d_tpa_voltages_in_cufft.resize(d_tpa_voltages_in.size(), {0.0f, 0.0f});
    _d_tpa_voltages_dedispersed.resize(d_tpa_voltages_in.size(), {0.0f, 0.0f});

    BOOST_LOG_TRIVIAL(debug)
        << "Input TPA voltages to dedisperse, d_tpa_voltages_in.size(): "
        << d_tpa_voltages_in.size();
    BOOST_LOG_TRIVIAL(debug)
        << "Output FTPA voltages to write to, d_ftpa_voltages_out.size(): "
        << d_ftpa_voltages_out.size();

    thrust::transform(d_tpa_voltages_in.begin(),
                      d_tpa_voltages_in.end(),
                      _d_tpa_voltages_in_cufft.begin(),
                      [=] __device__(char2 const& val) {
                          cufftComplex complex_val;
                          complex_val.x = val.x;
                          complex_val.y = val.y;
                          return complex_val;
                      });

    BOOST_LOG_TRIVIAL(debug) << "Transformed voltages to cufftComplex";

    cufftExecC2C(_config._fft_plan,
                 thrust::raw_pointer_cast(_d_tpa_voltages_in_cufft.data()),
                 thrust::raw_pointer_cast(_d_fpa_spectra.data()),
                 CUFFT_FORWARD);

    BOOST_LOG_TRIVIAL(debug) << "Executed forward FFT";

    BOOST_LOG_TRIVIAL(debug) << "freq_idx = " << freq_idx;
    BOOST_LOG_TRIVIAL(debug) << "dm_idx = " << dm_idx;

    multiply_by_chirp(_d_fpa_spectra,
                      _d_fpa_spectra,
                      freq_idx,
                      dm_idx); // operating in place..

    BOOST_LOG_TRIVIAL(debug) << "Multiplied by chirp";

    cufftExecC2C(_config._fft_plan,
                 thrust::raw_pointer_cast(_d_fpa_spectra.data()),
                 thrust::raw_pointer_cast(_d_tpa_voltages_dedispersed.data()),
                 CUFFT_INVERSE);

    BOOST_LOG_TRIVIAL(debug) << "Executed inverse FFT";

    std::size_t out_offset =
        freq_idx * _config.nantennas * _config.npols * (_config.gulp_samps);
    std::size_t discard_size =
        _config.nantennas * _config.npols * _config.overlap_samps / 2;

    BOOST_LOG_TRIVIAL(debug) << "Output offset to write from: " << out_offset;
    BOOST_LOG_TRIVIAL(debug) << "discard_size: " << discard_size;
    BOOST_LOG_TRIVIAL(debug)
        << "copying from input from " << discard_size << " to "
        << _d_tpa_voltages_dedispersed.size() - discard_size;
    BOOST_LOG_TRIVIAL(debug)
        << "Total elements copied: "
        << _d_tpa_voltages_dedispersed.size() - 2 * discard_size;
    BOOST_LOG_TRIVIAL(debug) << "Remaining space in output: "
                             << d_ftpa_voltages_out.size() - out_offset;
    BOOST_LOG_TRIVIAL(debug)
        << "copying to output from " << out_offset << " to "
        << out_offset + _d_tpa_voltages_dedispersed.size() - 2 * discard_size;

    std::size_t fft_size = _config.gulp_samps + _config.overlap_samps;

    thrust::transform(
        _d_tpa_voltages_dedispersed.begin() + discard_size,
        _d_tpa_voltages_dedispersed.end() - discard_size,
        d_ftpa_voltages_out.begin() + out_offset,
        [=] __device__(cufftComplex const& val) {
            char2 char2_val;
            char2_val.x = static_cast<char>(
                __float2int_rn(val.x / fft_size)); // scale the data back
            char2_val.y = static_cast<char>(__float2int_rn(val.y / fft_size));
            return char2_val;
        });

    BOOST_LOG_TRIVIAL(debug) << "Transformed dedispersed voltages to char2";
    d_ftpa_voltages_out.reference_dm(_config._h_dms[dm_idx]);
}

void CoherentDedisperser::multiply_by_chirp(

    thrust::device_vector<cufftComplex> const& _d_fpa_spectra_in,
    thrust::device_vector<cufftComplex>& _d_fpa_spectra_out,
    unsigned int freq_idx,
    unsigned int dm_idx)
{
    std::size_t total_chans =
        _config._d_ism_responses[dm_idx].size(); // all coarse + fine chans
    std::size_t fft_size =
        _config.gulp_samps + _config.overlap_samps; // ONLY FINE CHANS
    std::size_t response_offset = freq_idx * fft_size;

    BOOST_LOG_TRIVIAL(debug) << "Freq idx: " << freq_idx;
    BOOST_LOG_TRIVIAL(debug) << "_config.gulp_samps: " << _config.gulp_samps;
    BOOST_LOG_TRIVIAL(debug)
        << "_config.overlap_samps: " << _config.overlap_samps;
    BOOST_LOG_TRIVIAL(debug) << "response_offset: " << response_offset;
    BOOST_LOG_TRIVIAL(debug) << "total_chans: " << total_chans;
    BOOST_LOG_TRIVIAL(debug)
        << "chirp multiply input size: " << _d_fpa_spectra_in.size();
    BOOST_LOG_TRIVIAL(debug)
        << "chirp multiply output size: " << _d_fpa_spectra_out.size();

    dim3 blockSize(_config.nantennas * _config.npols);
    dim3 gridSize(fft_size / NCHANS_PER_BLOCK);
    kernels::dedisperse<<<gridSize, blockSize>>>(
        thrust::raw_pointer_cast(_config._d_ism_responses[dm_idx].data() +
                                 response_offset),
        thrust::raw_pointer_cast(_d_fpa_spectra_in.data()),
        thrust::raw_pointer_cast(_d_fpa_spectra_out.data()),
        total_chans);
}
} // namespace skyweaver
namespace skyweaver
{
namespace kernels
{

__global__ void dedisperse(cufftComplex const* __restrict__ _d_ism_response,
                           cufftComplex const* in,
                           cufftComplex* out,
                           unsigned total_chans)
{
    const unsigned pa_size = blockDim.x; // NANT * NPOL

    volatile __shared__ cufftComplex response[NCHANS_PER_BLOCK];

    const unsigned block_start_chan_idx =
        blockIdx.x * NCHANS_PER_BLOCK; // coarse chan idx

    const unsigned remainder =
        min(total_chans - block_start_chan_idx,
            NCHANS_PER_BLOCK); // how many channels to process

    for(int idx = threadIdx.x; idx < remainder; idx += pa_size) {
        cufftComplex const temp = _d_ism_response[block_start_chan_idx + idx];
        response[idx].x         = temp.x;
        response[idx].y         = temp.y;
    }

    __syncthreads();

    /**
    Each block processes NANT * NPOL in parallel (one per thread)
    Each thread processes NCHANS_PER_BLOCK channels sequentially, for a given
    (iant, ipol)
    **/

    for(int block_ichan_idx = 0; block_ichan_idx < remainder;
        ++block_ichan_idx) {
        const int chan_idx = (block_start_chan_idx + block_ichan_idx) *
                                 pa_size    // get to the correct chan_idx
                             + threadIdx.x; // get to the correct (iant, ipol)
        out[chan_idx] = cuCmulf(response[block_ichan_idx], in[chan_idx]);
    }
}

struct DMResponse {
    int num_coarse_chans;
    double low_freq;
    double coarse_chan_bw;
    double fine_chan_bw;
    double dm_prefactor;
    int num_fine_chans;

    DMResponse(int num_coarse_chans,
               int num_fine_chans,
               double low_freq,
               double coarse_chan_bw,
               double fine_chan_bw,
               double dm_prefactor)
        : num_coarse_chans(num_coarse_chans), num_fine_chans(num_fine_chans),
          low_freq(low_freq), coarse_chan_bw(coarse_chan_bw),
          fine_chan_bw(fine_chan_bw), dm_prefactor(dm_prefactor)
    {
    }

    __device__ inline cufftComplex operator()(int tid) const
    {
        int chan      = tid / num_fine_chans; // Coarse channel
        int fine_chan = tid % num_fine_chans; // fine channel

        double nu_0 = low_freq + chan * coarse_chan_bw;

        double nu = fine_chan * fine_chan_bw; // fine_chan_freq

        double phase_prefactor = nu * nu * dm_prefactor;
        double phase =
            phase_prefactor / ((nu_0 + nu) * nu_0 * nu_0); // precalculate
        cufftDoubleComplex weight;
        sincos(phase,
               &weight.y,
               &weight.x); // TODO: test if it is not approximate
        cufftComplex float_weight;
        float_weight.x = static_cast<float>(weight.x);
        float_weight.y = static_cast<float>(weight.y);

        return float_weight;
    }
};

} // namespace kernels

void get_dm_responses(CoherentDedisperserConfig& config,
                      double dm_prefactor,
                      thrust::device_vector<cufftComplex>& response)
{
    BOOST_LOG_TRIVIAL(debug) << "Generating DM responses";
    std::size_t fft_size = config.gulp_samps + config.overlap_samps;

    thrust::device_vector<int> indices(config.num_coarse_chans * fft_size);
    thrust::sequence(indices.begin(), indices.end());

    // store raw responses in a temporary variable
    thrust::device_vector<cufftComplex> temp_response(response.size());
    BOOST_LOG_TRIVIAL(warning) << "DOING FFTSHIFT" << std::endl;
    // Apply the DMResponse functor using thrust's transform
    thrust::transform(indices.begin(),
                      indices.end(),
                      temp_response.begin(),
                      kernels::DMResponse(config.num_coarse_chans,
                                          fft_size,
                                          config.low_freq,
                                          config.coarse_chan_bw,
                                          config.fine_chan_bw,
                                          dm_prefactor));

    // rotate the response to match cufft output
    
    if(fft_size %2 != 0) {
        throw std::runtime_error("FFT size must be even.");
    }

    std::size_t shift = fft_size / 2;

    for(auto it = temp_response.begin(), it_out = response.begin();
        it < temp_response.end();
        it += fft_size, it_out += fft_size) {
        thrust::copy(it, it + shift, it_out + shift);
        thrust::copy(it + shift, it + fft_size, it_out);
    }
}

} // namespace skyweaver
