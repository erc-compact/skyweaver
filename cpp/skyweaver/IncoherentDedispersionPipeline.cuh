#ifndef SKYWEAVER_INCOHERENTDEDISPERSIONPIPELINE_CUH
#define SKYWEAVER_INCOHERENTDEDISPERSIONPIPELINE_CUH

#include "skyweaver/AggregationBuffer.cuh"
#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/IncoherentDedisperser.cuh"
#include "skyweaver/ObservationHeader.hpp"
#include "skyweaver/PipelineConfig.hpp"
#include "skyweaver/Timer.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace skyweaver
{

template <typename InputType, typename OutputType, typename Handler>
class IncoherentDedispersionPipeline
{
  public:
    typedef AggregationBuffer<InputType> AggBufferType;
    typedef std::vector<std::unique_ptr<AggBufferType>> AggBufferVector;
    typedef TFBPowersD<InputType> InputVectorType;   // CoherentBeam Data
    typedef TFBPowersH<InputType> InputVectorTypeH;  // CoherentBeam Data
    typedef TDBPowersH<OutputType> OutputVectorType; // Dedispersered Data
    typedef IncoherentDedisperser DedisperserType;
    typedef std::vector<std::unique_ptr<DedisperserType>> DedisperserVector;

  public:
    IncoherentDedispersionPipeline(PipelineConfig const& config,
                                   Handler& handler);
    ~IncoherentDedispersionPipeline();
    IncoherentDedispersionPipeline(IncoherentDedispersionPipeline const&) =
        delete;
    IncoherentDedispersionPipeline&
    operator=(IncoherentDedispersionPipeline const&) = delete;
    void init(ObservationHeader const& header);
    void operator()(InputVectorType const& data, std::size_t ref_dm_idx);

  private:
    void
    agg_buffer_callback(typename InputVectorTypeH::VectorType const& buffer,
                        std::size_t ref_dm_idx);

  private:
    PipelineConfig const& _config;
    Handler& _handler;
    AggBufferVector _agg_buffers;
    DedisperserVector _dedispersers;
    std::vector<std::vector<OutputVectorType>> _output_buffers;
    Timer _timer;
};

} // namespace skyweaver

#include "skyweaver/detail/IncoherentDedispersionPipeline.cu"

#endif // SKYWEAVER_INCOHERENTDEDISPERSIONPIPELINE_CUH
