#ifndef SKYWEAVER_BEAMFORMERPIPELINE_CUH
#define SKYWEAVER_BEAMFORMERPIPELINE_CUH

#include "cuda.h"
#include "psrdada_cpp/raw_bytes.hpp"
#include "skyweaver/BufferedDispenser.cuh"
#include "skyweaver/CoherentBeamformer.cuh"
#include "skyweaver/CoherentDedisperser.cuh"
#include "skyweaver/DelayManager.cuh"
#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/IncoherentBeamformer.cuh"
#include "skyweaver/ObservationHeader.hpp"
#include "skyweaver/PipelineConfig.hpp"
#include "skyweaver/StatisticsCalculator.cuh"
#include "skyweaver/Timer.hpp"
#include "skyweaver/Transposer.cuh"
#include "skyweaver/WeightsManager.cuh"

#include <functional>
#include <thread>

namespace skyweaver
{

template <typename CBHandler,
          typename IBHandler,
          typename StatsHandler,
          typename BeamformerTraits>
class BeamformerPipeline
{
  public:
    using VoltageVectorTypeH   = TAFTPVoltagesH<char2>;
    using VoltageVectorTypeD   = TAFTPVoltagesD<char2>;
    using TimeType             = long double;
    using CoherentBeamformer   = CoherentBeamformer<BeamformerTraits>;
    using IncoherentBeamformer = IncoherentBeamformer<BeamformerTraits>;

  public:
    /**
     * @brief      Constructs the pipeline object.
     *
     * @param      config                  The pipeline configuration
     * @param      cb_handler               DADA write client for output
     * coherent beam data buffer
     * @param      ib_handler               DADA write client for output
     * incoherent beam data buffer
     * @param[in]  input_data_buffer_size  The input DADA buffer block size
     */
    BeamformerPipeline(PipelineConfig const& config,
                       CBHandler& cb_handler,
                       IBHandler& ib_handler,
                       StatsHandler& stats_handler);
    ~BeamformerPipeline();
    BeamformerPipeline(BeamformerPipeline const&) = delete;

    /**
     * @brief      Initialise the pipeline with a DADA header block
     *
     * @param      header  A RawBytes object wrapping the DADA header block
     */
    void init(ObservationHeader const& header, long double utc_offset = 0.0);

    /**
     * @brief      Process the data in a DADA data buffer
     *
     * @param      data  A RawBytes object wrapping the DADA data block
     */
    bool operator()(VoltageVectorTypeH const& data);

  private:
    void process();

  private:
    PipelineConfig const& _config;
    CoherentDedisperserConfig _dedisperser_config;

    // Data info
    ObservationHeader _header;
    int _nbeamsets;

    // Handlers
    CBHandler& _cb_handler;
    std::unique_ptr<std::thread> _cb_handler_thread;
    IBHandler& _ib_handler;
    StatsHandler& _stats_handler;

    // Streams
    cudaStream_t _h2d_copy_stream;
    cudaStream_t _processing_stream;
    cudaStream_t _d2h_copy_stream;

    // Pipeline components
    std::unique_ptr<DelayManager> _delay_manager;
    std::unique_ptr<WeightsManager> _weights_manager;
    std::unique_ptr<StatisticsCalculator> _stats_manager;
    std::unique_ptr<Transposer> _transposer;
    std::unique_ptr<CoherentBeamformer> _coherent_beamformer;
    std::unique_ptr<IncoherentBeamformer> _incoherent_beamformer;
    std::unique_ptr<CoherentDedisperser> _coherent_dedisperser;
    std::unique_ptr<BufferedDispenser> _dispenser;

    // Buffers
    TAFTPVoltagesD<char2> _taftp_from_host;
    FTPAVoltagesD<char2> _ftpa_post_transpose;
    FTPAVoltagesD<char2> _ftpa_dedispersed;
    DoubleDescribedVector<TFBPowersD<typename BeamformerTraits::QuantisedPowerType>> _btf_cbs;
    BTFPowersD<typename BeamformerTraits::RawPowerType> _tf_ib_raw;
    BTFPowersD<typename BeamformerTraits::QuantisedPowerType> _tf_ib;

    // Variable
    long double _unix_timestamp;
    std::size_t _sample_clock_tick_per_block;
    std::size_t _call_count;
    long double _utc_offset;

    // Trackers
    Timer _timer;
};

} // namespace skyweaver

#include "skyweaver/detail/BeamformerPipeline.cu"

#endif // SKYWEAVER_BEAMFORMERPIPELINE_CUH
