#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/IncoherentBeamformer.cuh"
#include "skyweaver/beamformer_utils.cuh"
#include "skyweaver/types.cuh"

#include <cassert>

#define ACC_BUFFER_SIZE 64

namespace skyweaver
{
namespace kernels
{

template <typename BfTraits>
__global__ void icbf_ftpa_general_k(
    char2 const* __restrict__ ftpa_voltages,
    typename BfTraits::RawPowerType* __restrict__ tf_powers_raw,
    typename BfTraits::QuantisedPowerType* __restrict__ tf_powers,
    float const* __restrict__ output_scale,
    float const* __restrict__ output_offset,
    float const* __restrict__ antenna_weights,
    int nsamples,
    int nbeamsets)
{
    // FTPA beamformer
    // Each thread loads the two pols of 1 antenna
    // Detects stokes
    // Sums over tscrunch and fscrunch
    // Total number of threads = N where N is a multiple of Nantennas
    // grid size = nchans/fscrunch, nsamples/tscrunch
    static_assert(SKYWEAVER_NPOL == 2,
                  "icbf_ftpa_general_k only works with dual pol data");

    const int a   = SKYWEAVER_NANTENNAS;
    const int pa  = SKYWEAVER_NPOL * a;
    const int tpa = nsamples * pa;

    // REMOVED volatile
    __shared__ typename BfTraits::RawPowerType acc_buffer[ACC_BUFFER_SIZE];

    for(int ii = threadIdx.x; ii < ACC_BUFFER_SIZE; ii += blockDim.x) {
        acc_buffer[ii] = BfTraits::zero_power;
    }
    __syncthreads();

    const int output_t_idx = blockIdx.x;
    const int output_f_idx = blockIdx.y;
    const int input_t_idx  = output_t_idx * SKYWEAVER_IB_TSCRUNCH;
    const int input_f_idx  = output_f_idx * SKYWEAVER_IB_FSCRUNCH;
    const int a_idx        = threadIdx.x;
    typename BfTraits::RawPowerType power = BfTraits::zero_power;
    for(int f_idx = input_f_idx; f_idx < input_f_idx + SKYWEAVER_CB_FSCRUNCH;
        ++f_idx) {
        for(int t_idx = input_t_idx;
            t_idx < input_t_idx + SKYWEAVER_IB_TSCRUNCH;
            ++t_idx) {
            int input_p0_idx = f_idx * tpa + t_idx * pa + a_idx;
            int input_p1_idx = f_idx * tpa + t_idx * pa + a + a_idx;
            char2 p0_v       = ftpa_voltages[input_p0_idx];
            char2 p1_v       = ftpa_voltages[input_p1_idx];
            cuFloatComplex p0 =
                make_cuFloatComplex((float)p0_v.x, (float)p0_v.y);
            cuFloatComplex p1 =
                make_cuFloatComplex((float)p1_v.x, (float)p1_v.y);
            BfTraits::integrate_stokes(p0, p1, power);
        }
    }
    for(int beamset_idx = 0; beamset_idx < nbeamsets; ++beamset_idx) {
        acc_buffer[a_idx] = power * antenna_weights[beamset_idx * SKYWEAVER_NANTENNAS + a_idx];
        for(unsigned int ii = ACC_BUFFER_SIZE / 2; ii > 0; ii >>= 1) {
            __syncthreads();
            if(a_idx < ii) {
                acc_buffer[a_idx] = acc_buffer[a_idx] + acc_buffer[a_idx + ii];
            }
        }
        __syncthreads();
        if(a_idx == 0) {
            const typename BfTraits::RawPowerType power_f32 = acc_buffer[0];
            const int output_idx = beamset_idx * gridDim.x * gridDim.y +
                                   output_t_idx * gridDim.y + output_f_idx;
            const float scale =
                output_scale[beamset_idx * gridDim.y + output_f_idx];
            const float offset =
                output_offset[beamset_idx * gridDim.y + output_f_idx];
            tf_powers_raw[output_idx] = power_f32;
            tf_powers[output_idx] =
                BfTraits::clamp(BfTraits::rescale(power_f32, offset, scale));
        }
    }
}

} // namespace kernels

template <typename BfTraits>
IncoherentBeamformer<BfTraits>::IncoherentBeamformer(
    PipelineConfig const& config)
    : _config(config)
{
    BOOST_LOG_TRIVIAL(debug) << "Constructing IncoherentBeamformer instance";
}

template <typename BfTraits>
IncoherentBeamformer<BfTraits>::~IncoherentBeamformer()
{
}

template <typename BfTraits>
void IncoherentBeamformer<BfTraits>::beamform(
    VoltageVectorTypeD const& input,
    RawPowerVectorTypeD& output_raw,
    PowerVectorTypeD& output,
    ScalingVectorTypeD const& output_scale,
    ScalingVectorTypeD const& output_offset,
    ScalingVectorTypeD const& antenna_weights,
    int nbeamsets,
    cudaStream_t stream)
{
    // First work out nsamples and resize output if not done already
    BOOST_LOG_TRIVIAL(debug) << "Executing incoherent beamforming";
    
    if(nbeamsets <= 0) {
        throw std::runtime_error(
            "Number of beamsets must be greater than zero");
    }
    std::size_t nsamples = input.nsamples();
    std::size_t output_size =
        (input.size() / _config.nantennas() / _config.npol() /
         _config.ib_tscrunch() / _config.ib_fscrunch()) *
        nbeamsets;
    BOOST_LOG_TRIVIAL(debug) << "Resizing output buffer from " << output.size()
                             << " to " << output_size << " elements";
    output.resize({static_cast<std::size_t>(nbeamsets),
                    nsamples / _config.ib_tscrunch(),
                   _config.nchans() / _config.ib_fscrunch()});
    output.metalike(input);
    output.tsamp(input.tsamp() * _config.ib_tscrunch());
    output_raw.resize({static_cast<std::size_t>(nbeamsets),
                        nsamples / _config.ib_tscrunch(),
                       _config.nchans() / _config.ib_fscrunch()});
    output_raw.metalike(input);
    output_raw.tsamp(input.tsamp() * _config.ib_tscrunch());
    if(output_scale.size() !=
       (_config.nchans() / _config.ib_fscrunch()) * nbeamsets) {
        std::runtime_error("Unexpected number of channels in scaling vector");
    }
    if(output_offset.size() !=
       (_config.nchans() / _config.ib_fscrunch()) * nbeamsets) {
        std::runtime_error("Unexpected number of channels in offset vector");
    }
    if(antenna_weights.size() != _config.nantennas() * nbeamsets) {
        std::runtime_error(
            "Antenna weights vector is not nantennas x nbeamsets in size");
    }
    dim3 block(_config.nantennas());
    dim3 grid(nsamples / _config.ib_tscrunch(),
              _config.nchans() / _config.ib_fscrunch());
    char2 const* ftpa_voltages_ptr = thrust::raw_pointer_cast(input.data());
    float const* output_scale_ptr =
        thrust::raw_pointer_cast(output_scale.data());
    float const* output_offset_ptr =
        thrust::raw_pointer_cast(output_offset.data());
    float const* antenna_weights_ptr =
        thrust::raw_pointer_cast(antenna_weights.data());
    typename PowerVectorTypeD::value_type* tf_powers_ptr =
        thrust::raw_pointer_cast(output.data());
    typename RawPowerVectorTypeD::value_type* tf_powers_raw_ptr =
        thrust::raw_pointer_cast(output_raw.data());
    BOOST_LOG_TRIVIAL(debug) << "Executing incoherent beamforming kernel";
    BOOST_LOG_TRIVIAL(debug) << "Nbeamsets = " << nbeamsets;
    kernels::icbf_ftpa_general_k<BfTraits>
        <<<grid, block, 0, stream>>>(ftpa_voltages_ptr,
                                     tf_powers_raw_ptr,
                                     tf_powers_ptr,
                                     output_scale_ptr,
                                     output_offset_ptr,
                                     antenna_weights_ptr,
                                     static_cast<int>(nsamples),
                                     nbeamsets);
    CUDA_ERROR_CHECK(cudaStreamSynchronize(stream));
    BOOST_LOG_TRIVIAL(debug) << "Incoherent beamforming kernel complete";
}

template class IncoherentBeamformer<
    SingleStokesBeamformerTraits<StokesParameter::I>>;
template class IncoherentBeamformer<
    SingleStokesBeamformerTraits<StokesParameter::Q>>;
template class IncoherentBeamformer<
    SingleStokesBeamformerTraits<StokesParameter::U>>;
template class IncoherentBeamformer<
    SingleStokesBeamformerTraits<StokesParameter::V>>;
template class IncoherentBeamformer<
    StokesTraits<StokesParameter::Q, StokesParameter::U>>;
template class IncoherentBeamformer<
    StokesTraits<StokesParameter::I, StokesParameter::V>>;
template class IncoherentBeamformer<FullStokesBeamformerTraits>;

} // namespace skyweaver
