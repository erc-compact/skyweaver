#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/StatisticsCalculator.cuh"
#include "thrust/host_vector.h"

#include <cstring>
#include <ctime>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#define LOG2_SKYWEAVER_NSAMPLES_PER_HEAP 8

namespace skyweaver
{
namespace kernel
{

__device__ void accumulate(double power,
                           long long& n,
                           double& M1,
                           double& M2,
                           double& M3,
                           double& M4)
{
    long long n1 = n;
    n++;
    double delta    = power - M1;
    double delta_n  = delta / n;
    double delta_n2 = delta_n * delta_n;
    double term1    = delta * delta_n * n1;
    M1 += delta_n;
    M4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * M2 -
          4 * delta_n * M3;
    M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
    M2 += term1;
}

/**
 * @brief Calculate statistics for the given input data
 *
 * @param ftpa_voltages Voltages data in FTPA order (8,8-bit complex)
 * @param results An array of statistics objects to store the outputs
 * @param nsamples The total number of time samples in the T dimension
 *
 * @details The FTPA order input allows for coalesced reads for the input
 *          with the caveat that we are only reading 16 bits at a time.
 *          The output is currently stored as an array of structs. A struct
 *          of arrays would likely be more performant here but is uglier
 *          from a coding persepective.
 *
 */
__global__ void calculate_statistics(char2 const* __restrict__ ftpa_voltages,
                                     Statistics* __restrict__ results,
                                     int nsamples)
{
    // Will operate on FTPA data and calculate statistics for FPA
    const int channel_idx = blockIdx.x;
    const int pol_idx     = blockIdx.y;
    const int antenna_idx = threadIdx.x;
    const int npol        = gridDim.y;
    const int nantennas   = blockDim.x;
    const int tpa_size    = npol * nantennas * nsamples;
    const int offset =
        channel_idx * tpa_size + pol_idx * nantennas + antenna_idx;
    const int stride = npol * nantennas;
    double M1 = 0.0, M2 = 0.0, M3 = 0.0, M4 = 0.0;
    long long n = 0;
    for(int sample_idx = offset; sample_idx < (tpa_size + offset);
        sample_idx += stride) {
        char2 data = ftpa_voltages[sample_idx];
        accumulate(static_cast<double>(data.x), n, M1, M2, M3, M4);
        accumulate(static_cast<double>(data.y), n, M1, M2, M3, M4);
    }

    // Output is ordered in FPA order
    int output_idx =
        channel_idx * npol * nantennas + pol_idx * nantennas + antenna_idx;
    Statistics* output = &results[output_idx];
    output->mean       = M1;
    output->std        = sqrt(M2 / (n - 1.0));
    output->skew       = sqrt((double)n) * M3 / pow(M2, 1.5);
    output->kurtosis   = (double)n * M4 / (M2 * M2) - 3.0;
}

} // namespace kernel

StatisticsCalculator::StatisticsCalculator(PipelineConfig const& config,
                                           cudaStream_t stream)
    : _config(config), _stream(stream)
{
    BOOST_LOG_TRIVIAL(debug)
        << "Constructing new StatisticsCalculator instance";
}

StatisticsCalculator::~StatisticsCalculator()
{
    if(_stats_file.is_open()) {
        _stats_file.close();
    }
}

void StatisticsCalculator::calculate_statistics(
    FTPAVoltagesD<char2> const& ftpa_voltages)
{
    _stats_d.resize({ftpa_voltages.nchannels(),
                     ftpa_voltages.npol(),
                     ftpa_voltages.nantennas()});
    _stats_d.metalike(ftpa_voltages);
    _stats_d.tsamp(ftpa_voltages.tsamp() * ftpa_voltages.nsamples());
    _stats_h.resize({ftpa_voltages.nchannels(),
                     ftpa_voltages.npol(),
                     ftpa_voltages.nantennas()});
    _stats_h.metalike(ftpa_voltages);
    _stats_h.tsamp(ftpa_voltages.tsamp() * ftpa_voltages.nsamples());

    int fpa_size = _stats_h.size();
    if(ftpa_voltages.size() % fpa_size != 0) {
        throw std::runtime_error(
            "FTPA voltages are not a multiple of FPA size");
    }
    int nsamples = ftpa_voltages.size() / fpa_size;
    // call kernel
    char2 const* ftpa_voltages_ptr =
        thrust::raw_pointer_cast(ftpa_voltages.data());
    Statistics* stats_ptr = thrust::raw_pointer_cast(_stats_d.data());
    dim3 dimBlock(_stats_d.nantennas());
    dim3 dimGrid(_stats_d.nchannels(), _stats_d.npol());
    kernel::calculate_statistics<<<dimGrid, dimBlock, 0, _stream>>>(
        ftpa_voltages_ptr,
        stats_ptr,
        nsamples);
    CUDA_ERROR_CHECK(cudaStreamSynchronize(_stream));
    BOOST_LOG_TRIVIAL(debug) << "Finished running calculate_statistics kernel";
    // Copy statistics to host
    _stats_h = _stats_d;
    BOOST_LOG_TRIVIAL(debug) << "Copied input levels to host";
}

void StatisticsCalculator::update_scalings(
    ScalingVectorHType const& beamset_weights,
    int nbeamsets)
{
    // At this stage we have the standard deviations of each channel
    // available on the host (h_input_levels) To support post-fact rescaling
    // of the data it is the scales and offsets that must be preserved to
    // disk.
    const float weights_amp = 127.0f;

    std::size_t reduced_nchans_ib = _config.nchans() / _config.ib_fscrunch();
    std::size_t reduced_nchans_cb = _config.nchans() / _config.cb_fscrunch();

    // Offsets for the coherent beams
    _cb_offsets_d.resize(reduced_nchans_cb * nbeamsets);
    _cb_offsets_h.resize(reduced_nchans_cb * nbeamsets);

    // Scalings for the coherent beams
    _cb_scaling_d.resize(reduced_nchans_cb * nbeamsets);
    _cb_scaling_h.resize(reduced_nchans_cb * nbeamsets);

    // Offsets for the incoherent beam
    _ib_offsets_d.resize(reduced_nchans_ib * nbeamsets);
    _ib_offsets_h.resize(reduced_nchans_ib * nbeamsets);

    // Scalings for the incoherent beam
    _ib_scaling_d.resize(reduced_nchans_ib * nbeamsets);
    _ib_scaling_h.resize(reduced_nchans_ib * nbeamsets);

    for(std::uint32_t beamset_idx = 0; beamset_idx < nbeamsets; ++beamset_idx) {
        const float effective_nantennas = std::accumulate(
            &beamset_weights[_config.nantennas() * beamset_idx],
            &beamset_weights[_config.nantennas() * (beamset_idx + 1)],
            0.0f);

        // define function
        auto get_offset_cb = [&](float x, Statistics const& y) {
            float scale =
                std::pow(weights_amp * y.std *
                             std::sqrt(static_cast<float>(effective_nantennas)),
                         2);
            float dof = 2 * _config.cb_tscrunch() * _config.cb_fscrunch() *
                        _config.npol();
            return x + (scale * dof);
        };

        auto get_scale_cb = [&](float x, Statistics const& y) {
            float scale =
                std::pow(weights_amp * y.std *
                             std::sqrt(static_cast<float>(effective_nantennas)),
                         2);
            float dof = 2 * _config.cb_tscrunch() * _config.cb_fscrunch() *
                        _config.npol();
            return x + (scale * std::sqrt(2 * dof) / _config.output_level());
        };

        auto get_offset_ib = [&](float x, Statistics const& y) {
            float scale = std::pow(y.std, 2);
            float dof   = 2 * _config.ib_tscrunch() * _config.ib_fscrunch() *
                        effective_nantennas * _config.npol();
            return x + (scale * dof);
        };

        auto get_scale_ib = [&](float x, Statistics const& y) {
            float scale = std::pow(y.std, 2);
            float dof   = 2 * _config.ib_tscrunch() * _config.ib_fscrunch() *
                        effective_nantennas * _config.npol();
            return x + (scale * std::sqrt(2 * dof) / _config.output_level());
        };

        // CB scaling and  offsets
        // Note: we use effective_nantennas below instead of the
        // _config.nantennas() because we want scaling factors that are valid
        // for valid data. Using the _config.nantnnas() would skew the stats.
        const int scale_factor =
            _config.cb_fscrunch() * _config.npol() * effective_nantennas;

        for(std::uint32_t out_chan_idx = 0; out_chan_idx < reduced_nchans_cb;
            ++out_chan_idx) {
            const std::uint32_t oidx =
                beamset_idx * reduced_nchans_cb + out_chan_idx;
            const std::uint32_t f_idx = _config.cb_fscrunch() * out_chan_idx;
            const std::uint32_t pa    = _config.npol() * _config.nantennas();
            const std::uint32_t start_idx = f_idx * pa;
            const std::uint32_t end_idx = (f_idx + _config.cb_fscrunch()) * pa;
            _cb_offsets_h[oidx]         = std::accumulate(&_stats_h[start_idx],
                                                  &_stats_h[end_idx],
                                                  0.0f,
                                                  get_offset_cb) /
                                  scale_factor;
            _cb_scaling_h[oidx] = std::accumulate(&_stats_h[start_idx],
                                                  &_stats_h[end_idx],
                                                  0.0f,
                                                  get_scale_cb) /
                                  scale_factor;
            BOOST_LOG_TRIVIAL(debug)
                << "Coherent beam power offset (beamset " << beamset_idx
                << "): " << _cb_offsets_h[oidx];
            BOOST_LOG_TRIVIAL(debug)
                << "Coherent beam power scaling (beamset " << beamset_idx
                << "): " << _cb_scaling_h[oidx];
        }

        // scaling for incoherent beamformer
        for(std::uint32_t out_chan_idx = 0; out_chan_idx < reduced_nchans_ib;
            ++out_chan_idx) {
            const std::uint32_t oidx =
                beamset_idx * reduced_nchans_ib + out_chan_idx;
            const std::uint32_t f_idx = _config.ib_fscrunch() * out_chan_idx;
            const std::uint32_t pa    = _config.npol() * _config.nantennas();
            const std::uint32_t start_idx = f_idx * pa;
            const std::uint32_t end_idx = (f_idx + _config.ib_fscrunch()) * pa;
            _ib_offsets_h[oidx]         = std::accumulate(&_stats_h[start_idx],
                                                  &_stats_h[end_idx],
                                                  0.0f,
                                                  get_offset_ib) /
                                  scale_factor;
            _ib_scaling_h[oidx] = std::accumulate(&_stats_h[start_idx],
                                                  &_stats_h[end_idx],
                                                  0.0f,
                                                  get_scale_ib) /
                                  scale_factor;
            BOOST_LOG_TRIVIAL(debug)
                << "Incoherent beam power offset (beamset " << beamset_idx
                << "): "
                << _ib_offsets_h[beamset_idx * reduced_nchans_ib +
                                 out_chan_idx];
            BOOST_LOG_TRIVIAL(debug)
                << "Incoherent beam power scaling (beamset " << beamset_idx
                << "): "
                << _ib_scaling_h[beamset_idx * reduced_nchans_ib +
                                 out_chan_idx];
        }
    }
    // At this stage, all scaling vectors are available on the host and
    // could be written to disk.

    // Copying these back to the device
    _cb_offsets_d = _cb_offsets_h;
    _cb_scaling_d = _cb_scaling_h;
    _ib_offsets_d = _ib_offsets_h;
    _ib_scaling_d = _ib_scaling_h;
}

void StatisticsCalculator::dump_all_scalings() const
{
    // Data can be written to the same place as the transient dumps
    // Total volume will typically by 2048 * 4 * 2 * 2 = 32 kB per target
    std::time_t now = std::time(0);
    std::tm* now_tm = std::gmtime(&now);
    char timestamp[42];
    std::strftime(timestamp, 42, "%Y-%m-%dT%X", now_tm);
    dump_scalings(timestamp, "cb_offsets", _config.output_dir(), _cb_offsets_h);
    dump_scalings(timestamp, "cb_scaling", _config.output_dir(), _cb_scaling_h);
    dump_scalings(timestamp, "ib_offsets", _config.output_dir(), _ib_offsets_h);
    dump_scalings(timestamp, "ib_scaling", _config.output_dir(), _ib_scaling_h);
}

void StatisticsCalculator::dump_scalings(
    std::string const& timestamp,
    std::string const& tag,
    std::string const& path,
    thrust::host_vector<float> const& ar) const
{
    std::ofstream writer;
    std::string filename = path + "/" + timestamp + "_" + tag + "_" +
                           std::to_string(_config.centre_frequency()) +
                           "Hz.bin";
    writer.open(filename, std::ofstream::out | std::ofstream::binary);
    if(writer.is_open()) {
        BOOST_LOG_TRIVIAL(info) << "Opened output file " << filename;
    } else {
        std::stringstream error_message;
        error_message << "Could not open file " << filename;
        BOOST_LOG_TRIVIAL(error) << error_message.str();
        return;
    }
    writer.write((char*)ar.data(), ar.size() * sizeof(float));
    writer.close();
}

StatisticsCalculator::StatisticsVectorDType const&
StatisticsCalculator::statistics() const
{
    return _stats_d;
}

StatisticsCalculator::ScalingVectorDType const&
StatisticsCalculator::cb_offsets() const
{
    return _cb_offsets_d;
}

StatisticsCalculator::ScalingVectorDType const&
StatisticsCalculator::cb_scaling() const
{
    return _cb_scaling_d;
}

StatisticsCalculator::ScalingVectorDType const&
StatisticsCalculator::ib_offsets() const
{
    return _ib_offsets_d;
}

StatisticsCalculator::ScalingVectorDType const&
StatisticsCalculator::ib_scaling() const
{
    return _ib_scaling_d;
}

} // namespace skyweaver
