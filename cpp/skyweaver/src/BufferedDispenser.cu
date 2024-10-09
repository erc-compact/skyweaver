#include "skyweaver/BufferedDispenser.cuh"
using namespace skyweaver;
BufferedDispenser::BufferedDispenser(PipelineConfig const& config,
                                    CoherentDedisperserConfig const& dedisp_config,
                                     cudaStream_t stream)
    : _config(config), _stream(stream)
{
    this->_block_length_tpa =
        _config.nantennas() * _config.npol() * _config.gulp_length_samps();
    this->_max_delay_tpa =
        dedisp_config.overlap_samps * _config.nantennas() * _config.npol();

    BOOST_LOG_TRIVIAL(debug) << "BD -> Block length TPA: " << _block_length_tpa;
    BOOST_LOG_TRIVIAL(debug) << "BD -> Max delay TPA: " << _max_delay_tpa;

    // this->_d_prev_ftpa_voltages.resize(_nchans * _max_delay_tpa);

    _d_channeled_tpa_voltages.resize(_config.nchans());
    _d_prev_channeled_tpa_voltages.resize(_config.nchans());
    _first_hoard.resize(_config.nchans(), true);

    for(std::size_t i = 0; i < _config.nchans(); i++) {
        _d_channeled_tpa_voltages[i].resize(
            { _config.gulp_length_samps() +  dedisp_config.overlap_samps,
             _config.npol(),
             _config.nantennas()});
        _d_prev_channeled_tpa_voltages[i].resize(
            { dedisp_config.overlap_samps,
             _config.npol(),
             _config.nantennas()});
    }

    BOOST_LOG_TRIVIAL(debug) << "overlap samples: " << dedisp_config.overlap_samps;

   
}

void BufferedDispenser::hoard(FTPAVoltagesTypeD const& new_ftpa_voltages_in)
{
    auto const& freqs = new_ftpa_voltages_in.frequencies();
    for(std::size_t i = 0; i < _config.nchans(); i++) {
        _d_channeled_tpa_voltages[i].frequencies(freqs[i]);
    }

    typename FTPAVoltagesTypeD::value_type zeros{};
    for(std::size_t i = 0; i < _config.nchans(); i++) {
        if(_first_hoard[i]) { // if first time set overlaps as zeros
            BOOST_LOG_TRIVIAL(debug)
                << "BD -> Filling TPA voltages " << i
                << " with zeros up to length " << _max_delay_tpa;


            thrust::fill(_d_channeled_tpa_voltages[i].begin(),
                         _d_channeled_tpa_voltages[i].begin() + _max_delay_tpa,
                         zeros);
            _first_hoard[i] = false;

        } else { // first add corresponding overlap to output
            BOOST_LOG_TRIVIAL(debug) << "BD -> Copying previous voltages of size: " << _d_prev_channeled_tpa_voltages[i].size();
            thrust::copy(_d_prev_channeled_tpa_voltages[i].begin(),
                         _d_prev_channeled_tpa_voltages[i].end(),
                         _d_channeled_tpa_voltages[i].begin());

        }
        // then add the input data
        BOOST_LOG_TRIVIAL(debug) << "BD -> Copying new voltages of size: " << _block_length_tpa << " starting at " << i * _block_length_tpa ;
        thrust::copy(new_ftpa_voltages_in.begin() + i * _block_length_tpa,
                     new_ftpa_voltages_in.begin() + (i + 1) * _block_length_tpa,
                     _d_channeled_tpa_voltages[i].begin() + _max_delay_tpa);

        // update the overlap for the next hoard
        BOOST_LOG_TRIVIAL(debug) << "BD -> Updating overlap to the data between " << (i + 1) * _block_length_tpa - _max_delay_tpa << " and " << (i + 1) * _block_length_tpa;
        thrust::copy(new_ftpa_voltages_in.begin() +
                         (i + 1) * _block_length_tpa - _max_delay_tpa,
                     new_ftpa_voltages_in.begin() + (i + 1) * _block_length_tpa,
                     _d_prev_channeled_tpa_voltages[i].begin());
    }
}

typename BufferedDispenser::TPAVoltagesTypeD const&
BufferedDispenser::dispense(std::size_t chan_idx) const
{ // implements overlapped buffering of data
    return _d_channeled_tpa_voltages[chan_idx];
}
