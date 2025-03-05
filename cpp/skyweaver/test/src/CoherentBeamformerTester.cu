#include "psrdada_cpp/common.hpp"
#include "psrdada_cpp/cuda_utils.hpp"
#include "skyweaver/IncoherentBeamformer.cuh"
#include "skyweaver/beamformer_utils.cuh"
#include "skyweaver/skyweaver_constants.hpp"
#include "skyweaver/test/CoherentBeamformerTester.cuh"
#include "skyweaver/test/test_utils.cuh"

#include <cmath>
#include <complex>
#include <fstream>
#include <random>

namespace skyweaver
{
namespace test
{

template <typename BfTraits>
CoherentBeamformerTester<BfTraits>::CoherentBeamformerTester()
    : ::testing::Test(), _stream(0)
{
}

template <typename BfTraits>
CoherentBeamformerTester<BfTraits>::~CoherentBeamformerTester()
{
}

template <typename BfTraits>
void CoherentBeamformerTester<BfTraits>::SetUp()
{
    CUDA_ERROR_CHECK(cudaStreamCreate(&_stream));
}

template <typename BfTraits>
void CoherentBeamformerTester<BfTraits>::TearDown()
{
    CUDA_ERROR_CHECK(cudaStreamDestroy(_stream));
}

template <typename BfTraits>
void CoherentBeamformerTester<BfTraits>::beamformer_c_reference(
    VoltageVectorTypeH const& ftpa_voltages,
    WeightsVectorTypeH const& fbpa_weights,
    HostPowerVectorType& btf_powers,
    int nchannels,
    int tscrunch,
    int fscrunch,
    int nsamples,
    int nbeams,
    int nantennas,
    int npol,
    float const* scales,
    float const* offsets,
    float const* antenna_weights,
    int const* beamset_mapping)
{
    for(int channel_idx = 0; channel_idx < nchannels; channel_idx += fscrunch) {
        BOOST_LOG_TRIVIAL(debug)
            << "Beamformer C reference: "
            << static_cast<int>(100.0f * (channel_idx + 1.0f) / nchannels)
            << "% complete";
        for(int sample_idx = 0; sample_idx < nsamples; sample_idx += tscrunch) {
            for(int beam_idx = 0; beam_idx < nbeams; ++beam_idx) {
                const int beamset_idx = beamset_mapping[beam_idx];
                typename BfTraits::RawPowerType power    = BfTraits::zero_power;
                typename BfTraits::RawPowerType ib_power = BfTraits::zero_power;
                for(int sub_channel_idx = channel_idx;
                    sub_channel_idx < channel_idx + fscrunch;
                    ++sub_channel_idx) {
                    for(int sample_offset = 0; sample_offset < tscrunch;
                        ++sample_offset) {
                        cuFloatComplex p0_accumulator = {0, 0};
                        cuFloatComplex p1_accumulator = {0, 0};
                        // new loop
                        for(int antenna_idx = 0; antenna_idx < nantennas;
                            ++antenna_idx) {
                            const int ftpa_voltages_partial_idx =
                                nantennas * npol * nsamples * sub_channel_idx +
                                nantennas * npol *
                                    (sample_idx + sample_offset) +
                                +antenna_idx;
                            const int p0_idx = ftpa_voltages_partial_idx;
                            const int p1_idx =
                                ftpa_voltages_partial_idx + nantennas;
                            char2 ant_p0_c2       = ftpa_voltages[p0_idx];
                            char2 ant_p1_c2       = ftpa_voltages[p1_idx];
                            cuFloatComplex ant_p0 = make_cuFloatComplex(
                                static_cast<float>(ant_p0_c2.x),
                                static_cast<float>(ant_p0_c2.y));
                            cuFloatComplex ant_p1 = make_cuFloatComplex(
                                static_cast<float>(ant_p1_c2.x),
                                static_cast<float>(ant_p1_c2.y));
                            BfTraits::integrate_weighted_stokes(
                                ant_p0,
                                ant_p1,
                                ib_power,
                                antenna_weights[beamset_idx * nantennas +
                                                antenna_idx]);
                            int fbpa_weights_idx =
                                nantennas * nbeams * sub_channel_idx +
                                nantennas * beam_idx + antenna_idx;
                            char2 weight_c2 = fbpa_weights[fbpa_weights_idx];
                            cuFloatComplex weight = make_cuFloatComplex(
                                static_cast<float>(weight_c2.x),
                                static_cast<float>(weight_c2.y));
                            cuFloatComplex p0 = cuCmulf(weight, ant_p0);
                            cuFloatComplex p1 = cuCmulf(weight, ant_p1);
                            p0_accumulator    = cuCaddf(p0_accumulator, p0);
                            p1_accumulator    = cuCaddf(p1_accumulator, p1);
                        }
                        BfTraits::integrate_stokes(p0_accumulator,
                                                   p1_accumulator,
                                                   power);
                        // end new loop
                    }
                }
                const int output_sample_idx = sample_idx / tscrunch;
                const int output_chan_idx   = channel_idx / fscrunch;
                const int nchans_out        = nchannels / fscrunch;

                /* For BTF outputs
                const int nsamps_out        = nsamples / tscrunch;
                const int tf_size           = nsamps_out * nchans_out;
                const int output_idx        = beam_idx * tf_size +
                                              output_sample_idx * nchans_out +
                                              output_chan_idx;
                */
                // For TFB outputs
                const int output_idx = output_sample_idx * nbeams * nchans_out +
                                       output_chan_idx * nbeams + beam_idx;

                const int scloff_idx =
                    beamset_idx * nchans_out + output_chan_idx;
#if SKYWEAVER_IB_SUBTRACTION
                typename BfTraits::RawPowerType powerf32 =
                    BfTraits::ib_subtract(power,
                                          ib_power,
                                          16129.0f,
                                          scales[scloff_idx]);
#else
                typename BfTraits::RawPowerType powerf32 =
                    BfTraits::rescale(power,
                                      offsets[scloff_idx],
                                      scales[scloff_idx]);
#endif // SKYWEAVER_IB_SUBTRACTION
                btf_powers[output_idx] = BfTraits::clamp(powerf32);
            }
        }
    }
}

template <typename BfTraits>
void CoherentBeamformerTester<BfTraits>::compare_against_host(
    VoltageVectorTypeD const& ftpa_voltages_gpu,
    WeightsVectorTypeD const& fbpa_weights_gpu,
    ScalingVectorTypeD const& scales_gpu,
    ScalingVectorTypeD const& offsets_gpu,
    ScalingVectorTypeD const& antenna_weights,
    MappingVectorTypeD const& beamset_mapping,
    PowerVectorTypeD& btf_powers_gpu,
    int nsamples)
{
    VoltageVectorTypeH ftpa_voltages_host = ftpa_voltages_gpu;
    WeightsVectorTypeH fbpa_weights_host  = fbpa_weights_gpu;
    HostPowerVectorType btf_powers_cuda   = btf_powers_gpu;
    HostPowerVectorType btf_powers_host;
    btf_powers_host.like(btf_powers_gpu);

    HostScalingVectorType scales               = scales_gpu;
    HostScalingVectorType offsets              = offsets_gpu;
    HostScalingVectorType antenna_weights_host = antenna_weights;
    MappingVectorTypeH beamset_mapping_host    = beamset_mapping;

    beamformer_c_reference(
        ftpa_voltages_host,
        fbpa_weights_host,
        btf_powers_host,
        _config.nchans(),
        _config.cb_tscrunch(),
        _config.cb_fscrunch(),
        nsamples,
        _config.nbeams(),
        _config.nantennas(),
        _config.npol(),
        thrust::raw_pointer_cast(scales.data()),
        thrust::raw_pointer_cast(offsets.data()),
        thrust::raw_pointer_cast(antenna_weights_host.data()),
        thrust::raw_pointer_cast(beamset_mapping_host.data()));
    for(int ii = 0; ii < btf_powers_host.size(); ++ii) {
        expect_near(btf_powers_host[ii], btf_powers_cuda[ii], 1);
    }
}

typedef ::testing::Types<SingleStokesBeamformerTraits<StokesParameter::I>,
                         SingleStokesBeamformerTraits<StokesParameter::Q>,
                         SingleStokesBeamformerTraits<StokesParameter::U>,
                         SingleStokesBeamformerTraits<StokesParameter::V>,
                         StokesTraits<StokesParameter::Q, StokesParameter::U>,
                         StokesTraits<StokesParameter::I, StokesParameter::V>,
                         FullStokesBeamformerTraits>
    StokesTypes;
TYPED_TEST_SUITE(CoherentBeamformerTester, StokesTypes);

/**
Something fishy with this test.
When run on a debug build (-DCMAKE_BUILD_TYPE=DEBUG) the powers returned
from the GPU kernel are always zero. This seems to be due to the kernel
not being executed in this mode.
//TODO: Work out what is going on here with nsys or ncu
*/
TYPED_TEST(CoherentBeamformerTester, representative_noise_test)
{
    using BfTraits = typename TestFixture::BfTraitsType;
    using CBT      = CoherentBeamformerTester<BfTraits>;
    auto& config   = this->_config;

#if SKYWEAVER_IB_SUBTRACTION
    BOOST_LOG_TRIVIAL(info) << "Running with IB subtraction";
#else
    BOOST_LOG_TRIVIAL(info) << "Running without IB subtraction";
#endif
    const int nbeamsets = 2;

    // Calculate the CB scaling values.
    const float input_level = 32.0f;
    const double pi         = std::acos(-1);
    config.output_level(input_level);
    float scale =
        std::pow(127.0f * input_level *
                     std::sqrt(static_cast<float>(config.nantennas())),
                 2);
    float dof = 2 * config.cb_tscrunch() * config.cb_fscrunch() * config.npol();
    float offset_val = (scale * dof);
    float scale_val  = (scale * std::sqrt(2 * dof) / config.output_level());

    // Set constant scales and offsets for all channels and beamsets
    typename CBT::ScalingVectorTypeD cb_scales(
        config.nchans() / config.cb_fscrunch() * nbeamsets,
        scale_val);
    typename CBT::ScalingVectorTypeD cb_offsets(
        config.nchans() / config.cb_fscrunch() * nbeamsets,
        offset_val);

    // Map all beams to the first beamset by default
    typename CBT::MappingVectorTypeD beamset_mapping(config.nbeams(), 0);

    // Enable all antennas in all beamsets
    typename CBT::ScalingVectorTypeD beamset_weights(config.nantennas() *
                                                         nbeamsets,
                                                    1.0f);

    // Weight half of the antennas to zero in beamset 1
    for (int ii =  config.nantennas()/2 ; ii < config.nantennas(); ++ii) {
        beamset_weights[config.nantennas() + ii] = 0.0f;
    }

    /**
    This currently causes tests to fail even though the weights are the same for
    all beamsets. Either the reference or the implementation are using the
    mapping incorrectly.

    Issue is that the incoherent beamformer is not populating beamsets > 0
    */
    beamset_mapping[0] = 1;
    beamset_mapping[config.nbeams() - 1] = 1;

    BOOST_LOG_TRIVIAL(info) << "CB scaling: " << scale_val;
    BOOST_LOG_TRIVIAL(info) << "CB offset: " << offset_val;

    std::default_random_engine generator;
    generator.seed(0x2333425234);
    std::normal_distribution<float> normal_dist(0.0, input_level);
    std::uniform_real_distribution<float> uniform_dist(0.0, 2 * pi);

    CoherentBeamformer<BfTraits> coherent_beamformer(config);
    IncoherentBeamformer<BfTraits> incoherent_beamformer(config);

    std::size_t ntimestamps = max(
        1L,
        SKYWEAVER_CB_PACKET_SIZE / (config.nchans() / config.cb_fscrunch()) /
            (config.nsamples_per_heap() / config.cb_tscrunch()));
    ntimestamps =
        max(ntimestamps,
            SKYWEAVER_CB_NSAMPLES_PER_BLOCK / config.nsamples_per_heap());
    std::size_t input_size =
        (ntimestamps * config.nantennas() * config.nchans() *
         config.nsamples_per_heap() * config.npol());
    BOOST_LOG_TRIVIAL(info) << "FTPA input dims: " << config.nchans() << ", "
                            << ntimestamps * config.nsamples_per_heap() << ", "
                            << config.npol() << ", " << config.nantennas();
    int nsamples = config.nsamples_per_heap() * ntimestamps;

    std::size_t weights_size =
        config.nantennas() * config.nchans() * config.nbeams();

    typename CBT::VoltageVectorTypeH ftpa_voltages_host(
        {config.nchans(),
         ntimestamps * config.nsamples_per_heap(),
         config.npol(),
         config.nantennas()});
    for(int ii = 0; ii < ftpa_voltages_host.size(); ++ii) {
        ftpa_voltages_host[ii].x =
            static_cast<int8_t>(std::lround(normal_dist(generator)));
        ftpa_voltages_host[ii].y =
            static_cast<int8_t>(std::lround(normal_dist(generator)));
    }

    typename CBT::WeightsVectorTypeH fbpa_weights_host(weights_size, char2{127, 0});
    
    for (int f_idx = 0; f_idx < config.nchans(); ++f_idx) {
        for (int b_idx = 0; b_idx < config.nbeams(); ++b_idx) {
            for (int a_idx = 0; a_idx < config.nantennas(); ++a_idx) {
                std::size_t beamset_idx = beamset_mapping[b_idx];
                float scale = beamset_weights[beamset_idx * SKYWEAVER_NANTENNAS + a_idx];
                std::complex<double> val = scale * 127.0f * 
                    std::exp(std::complex<float>(0.0f, uniform_dist(generator)));
                std::size_t idx = f_idx * config.nbeams() * config.nantennas() + b_idx * config.nantennas() + a_idx;
                fbpa_weights_host[idx].x = static_cast<int8_t>(std::lround(val.real()));
                fbpa_weights_host[idx].y = static_cast<int8_t>(std::lround(val.imag()));
            }
        }
    }
    
    typename CBT::VoltageVectorTypeD ftpa_voltages_gpu = ftpa_voltages_host;
    typename CBT::WeightsVectorTypeD fbpa_weights_gpu  = fbpa_weights_host;
    typename CBT::PowerVectorTypeD tfb_powers_gpu;

    // Note that below even though this is for the IB we have to use the
    // CB scrunching parameters to make sure we get the right data out.
    float ib_scale = std::pow(input_level, 2);
    float ib_dof   = 2 * config.cb_tscrunch() * config.cb_fscrunch() *
                   config.nantennas() * config.npol();
    float ib_power_offset = ib_scale * ib_dof;
    float ib_power_scaling =
        ib_scale * std::sqrt(2 * ib_dof) / config.output_level();
    typename CBT::ScalingVectorTypeD ib_scales(
        config.nchans() / config.cb_fscrunch() * nbeamsets,
        ib_power_scaling);
    typename CBT::ScalingVectorTypeD ib_offset(
        config.nchans() / config.cb_fscrunch() * nbeamsets,
        ib_power_offset);
    typename CBT::IBPowerVectorTypeD tf_powers_gpu;
    typename CBT::RawIBPowerVectorTypeD tf_powers_raw_gpu;

    // dump_device_vector(ftpa_voltages_gpu, "ftpa_voltages_gpu.bin");
    incoherent_beamformer.beamform(ftpa_voltages_gpu,
                                   tf_powers_raw_gpu,
                                   tf_powers_gpu,
                                   ib_scales,
                                   ib_offset,
                                   beamset_weights,
                                   nbeamsets,
                                   this->_stream);
    // dump_device_vector(tf_powers_raw_gpu, "tf_powers_raw_gpu.bin");
    // dump_device_vector(fbpa_weights_gpu, "fbpa_weights_gpu.bin");
    coherent_beamformer.beamform(ftpa_voltages_gpu,
                                 fbpa_weights_gpu,
                                 cb_scales,
                                 cb_offsets,
                                 beamset_mapping,
                                 tf_powers_raw_gpu,
                                 tfb_powers_gpu,
                                 nbeamsets,
                                 this->_stream);
    // dump_device_vector(tfb_powers_gpu, "tfb_powers_gpu.bin");
    this->compare_against_host(ftpa_voltages_gpu,
                               fbpa_weights_gpu,
                               cb_scales,
                               cb_offsets,
                               beamset_weights,
                               beamset_mapping,
                               tfb_powers_gpu,
                               nsamples);
}

} // namespace test
} // namespace skyweaver
