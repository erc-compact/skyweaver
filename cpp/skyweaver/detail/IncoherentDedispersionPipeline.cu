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
        int max_delay = _dedispersers[block_idx]->max_sample_delay();
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

    // ceil division
    _n_tdb_files = (_config.nbeams() + _config.nbeams_per_file() - 1) / _config.nbeams_per_file();
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
                        std::size_t ref_dm_idx)
{
    BOOST_LOG_TRIVIAL(debug)
        << "Agg buffer callback called for ref_dm_idx = " << ref_dm_idx;
    _timer.start("incoherent dedispersion");
    _dedispersers[ref_dm_idx]->dedisperse(buffer, _output_buffers[ref_dm_idx]);
    _timer.stop("incoherent dedispersion");
    BOOST_LOG_TRIVIAL(debug) << "Dedispersion complete, calling handler";
    BOOST_LOG_TRIVIAL(debug) << _output_buffers[ref_dm_idx].vector().size();
    auto const& plan = _config.ddplan();

    // Set the correct DMs on the block
    _output_buffers[ref_dm_idx].dms(plan[ref_dm_idx].incoherent_dms);
    _output_buffers[ref_dm_idx].reference_dm(plan[ref_dm_idx].coherent_dm);
    _output_buffers[ref_dm_idx].frequencies(_config.channel_frequencies().front());

    BOOST_LOG_TRIVIAL(debug) << "setting centre frequency to " << std::setprecision(15) <<  _output_buffers[ref_dm_idx].frequencies()[0];
    BOOST_LOG_TRIVIAL(debug) << "Passing output buffer to handler: "
                             << _output_buffers[ref_dm_idx].describe();
    _timer.start("file writing");

    if (_n_tdb_files > 1)
    {
        const std::size_t ndms     = _output_buffers[ref_dm_idx].ndms();
        const std::size_t nbeams   = _output_buffers[ref_dm_idx].nbeams();
        const std::size_t nsamples = _output_buffers[ref_dm_idx].nsamples();

        std::size_t nbeams_per_file = _config.nbeams_per_file();
        std::size_t td_input_offset = 0;
        std::size_t td_output_offset = 0;
        std::size_t b_offset = 0;
        std::size_t tdb_file_idx = 0;

        while (b_offset < nbeams)
        {
            if (b_offset + nbeams_per_file > nbeams)
            {
                nbeams_per_file = nbeams - b_offset;
            }

            _beamsplit_buffer.metalike(_output_buffers[ref_dm_idx]);
            _beamsplit_buffer.resize({nsamples, ndms, nbeams_per_file});

            for (std::size_t tdidx = 0; tdidx < nsamples * ndms; ++tdidx)
            {
                td_input_offset = tdidx * nbeams;
                td_output_offset = tdidx * nbeams_per_file;

                std::copy(&_output_buffers[ref_dm_idx][td_input_offset + b_offset],
                          &_output_buffers[ref_dm_idx][td_input_offset + b_offset + nbeams_per_file],
                          &_beamsplit_buffer[td_output_offset]);
            }

            _beamsplit_buffer.beam0_idx(b_offset);
            _handler(_beamsplit_buffer, tdb_file_idx * _config.coherent_dms().size() + ref_dm_idx);

            b_offset += nbeams_per_file;
            tdb_file_idx++;
        }
    }
    else
    {
        _handler(_output_buffers[ref_dm_idx], ref_dm_idx);
    }
    _timer.stop("file writing");
}

template <typename InputType, typename OutputType, typename Handler>
void IncoherentDedispersionPipeline<InputType, OutputType, Handler>::init(
    ObservationHeader const& header)
{
    std::size_t nchans = _config.nchans();
    long double chbw   = _config.bandwidth() / _config.nchans();
    long double tsamp  = _config.cb_tscrunch() * 1.0 / chbw;
    std::vector<long double> dm_delays(_config.coherent_dms().size());
    for(std::size_t ref_dm_idx = 0; ref_dm_idx < _config.coherent_dms().size();
        ++ref_dm_idx) {
        dm_delays[ref_dm_idx] =
            _dedispersers[ref_dm_idx]->max_sample_delay() * tsamp;
    }
    _handler.init(header);
}

template <typename InputType, typename OutputType, typename Handler>
void IncoherentDedispersionPipeline<InputType, OutputType, Handler>::operator()(
    InputVectorType const& data,
    std::size_t ref_dm_idx)
{
    _output_buffers[ref_dm_idx].metalike(data);
    _output_buffers[ref_dm_idx].tsamp(data.tsamp() *
                                      _config.ddplan()[ref_dm_idx].tscrunch);
    _output_buffers[ref_dm_idx].utc_offset(
        data.utc_offset() +
        _dedispersers[ref_dm_idx]->max_sample_delay() * data.tsamp());
    _agg_buffers[ref_dm_idx]->push_back(data.vector());
}

} // namespace skyweaver
