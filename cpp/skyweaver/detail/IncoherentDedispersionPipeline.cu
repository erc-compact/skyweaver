#include "skyweaver/IncoherentDedispersionPipeline.cuh"

namespace skyweaver
{

template <typename InputType, typename OutputType, typename Handler>
IncoherentDedispersionPipeline<InputType, OutputType, Handler>::
    IncoherentDedispersionPipeline(PipelineConfig const& config,
                                   Handler& handler)
    : _config(config), _handler(handler)
{
    BOOST_LOG_NAMED_SCOPE("IncoherentDedispersionPipeline")
    BOOST_LOG_TRIVIAL(debug) << "Preparing incoherent dedispersion pipeline";
    auto& plan         = _config.ddplan();
    auto const& blocks = plan.blocks();
    _dedispersers.resize(blocks.size());
    _agg_buffers.resize(blocks.size());
    _output_buffers.resize(blocks.size());
    for(std::size_t block_idx = 0; block_idx < blocks.size(); ++block_idx) {
        _dedispersers[block_idx].reset(
            new DedisperserType(_config,
                                blocks[block_idx].incoherent_dms,
                                blocks[block_idx].tscrunch));
        int max_delay = _dedispersers[block_idx]->max_delay();
        BOOST_LOG_TRIVIAL(debug)
            << "Created dedisperser for block = " << block_idx
            << " (max_delay = " << max_delay
            << ", tscrunch = " << blocks[block_idx].tscrunch << ")";
        std::size_t dispatch_size = std::max(max_delay * 10, 8192);
        // The dispatch size needs rounded to the next multiple of tscrunch
        dispatch_size -= dispatch_size % blocks[block_idx].tscrunch;

        std::size_t batch_size = _config.nbeams() * _config.nchans();
        _agg_buffers[block_idx].reset(new AggBufferType(
            std::bind(&IncoherentDedispersionPipeline::agg_buffer_callback,
                      this,
                      std::placeholders::_1,
                      block_idx),
            dispatch_size,
            max_delay,
            batch_size));
        BOOST_LOG_TRIVIAL(debug)
            << "Created aggregation buffer for block = " << block_idx
            << " (dispatch_size = " << dispatch_size << ", "
            << "overlap_size = " << max_delay << ", "
            << "batch_size = " << batch_size << ")";
        // TODO: give these sensible numbers
        // TODO: work out how to get the max delay into the handler
    }
}

template <typename InputType, typename OutputType, typename Handler>
IncoherentDedispersionPipeline<InputType, OutputType, Handler>::
    ~IncoherentDedispersionPipeline()
{
    _timer.show_all_timings();
}

template <typename InputType, typename OutputType, typename Handler>
void IncoherentDedispersionPipeline<InputType, OutputType, Handler>::
    agg_buffer_callback(typename InputVectorTypeH::VectorType const& buffer,
                        std::size_t block_idx)
{
    BOOST_LOG_TRIVIAL(debug)
        << "Agg buffer callback called for block_idx = " << block_idx;
    _timer.start("incoherent dedispersion");
    _dedispersers[block_idx]->dedisperse(buffer, _output_buffers[block_idx]);
    _timer.stop("incoherent dedispersion");
    BOOST_LOG_TRIVIAL(debug) << "Dedispersion complete, calling handler";
    BOOST_LOG_TRIVIAL(debug) << _output_buffers[block_idx].vector().size();
    auto const& plan = _config.ddplan();
    // Set the correct tsamp on the block
    _output_buffers[block_idx].tsamp(_output_buffers[block_idx].tsamp() *
                                     plan[block_idx].tscrunch);
    // Set the correct DMs on the block
    _output_buffers[block_idx].dms(plan[block_idx].incoherent_dms);
    _output_buffers[block_idx].reference_dm(plan[block_idx].coherent_dm);
    _output_buffers[block_idx].frequencies({_config.centre_frequency()});

    BOOST_LOG_TRIVIAL(debug) << "Passing output buffer to handler: "
                            << _output_buffers[block_idx].describe();
    _timer.start("file writing");
    _handler(_output_buffers[block_idx], block_idx);
    _timer.stop("file writing");
}

template <typename InputType, typename OutputType, typename Handler>
void IncoherentDedispersionPipeline<InputType, OutputType, Handler>::init(
    ObservationHeader const& header)
{
    std::size_t nchans = _config.channel_frequencies().size();
    long double chbw   = _config.bandwidth() / _config.nchans();
    long double tsamp  = _config.cb_tscrunch() / chbw;
    std::vector<long double> dm_delays(_config.coherent_dms().size());
    for(std::size_t dm_idx = 0; dm_idx < _config.coherent_dms().size();
        ++dm_idx) {
        dm_delays[dm_idx] = _dedispersers[dm_idx]->max_delay() * tsamp;
    }
    _handler.init(header);
}

template <typename InputType, typename OutputType, typename Handler>
void IncoherentDedispersionPipeline<InputType, OutputType, Handler>::operator()(
    InputVectorType const& data,
    std::size_t dm_idx)
{
    _output_buffers[dm_idx].metalike(data);
    _output_buffers[dm_idx].tsamp(data.tsamp() *
                                  _config.ddplan()[dm_idx].tscrunch);
    _agg_buffers[dm_idx]->push_back(data.vector());
}

} // namespace skyweaver