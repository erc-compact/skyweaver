#include "skyweaver/BufferedDispenser.cuh"
using namespace skyweaver;
BufferedDispenser::BufferedDispenser(PipelineConfig const& config) : _config(config) {
        this-> _block_length_tpa = _config.nantennas() * _config.npol() * _config.gulp_length_samps();
        this->_kernel_length_tpa = _config.dedisp_kernel_length_samps() * _config.nantennas() * _config.npol();

        // this->_d_prev_ftpa_voltages.resize(_nchans * _kernel_length_tpa);

        _d_channeled_tpa_voltages.resize(_config.nchans());
        _d_prev_channeled_tpa_voltages.resize(_config.nchans());

        for (std::size_t i = 0; i < _config.nchans(); i++){
            _d_channeled_tpa_voltages[i].resize(_block_length_tpa + _kernel_length_tpa);
            _d_prev_channeled_tpa_voltages[i].resize(_kernel_length_tpa);
        }
    }

void BufferedDispenser::hoard(DeviceVoltageType const& new_ftpa_voltages_in, cudaStream_t stream){

    _stream = stream;
    for (std::size_t i = 0; i < _config.nchans(); i++){

        if(!_d_prev_channeled_tpa_voltages){ // if first time set overlaps as zeros
            thrust::fill(_d_channeled_tpa_voltages[i].begin(),
                         _d_channeled_tpa_voltages[i].begin() + kernel_length_tpa, 
                         0); 
        }
        else{ // first add corresponding overlap to output 
            thrust::copy(_d_prev_channeled_tpa_voltages[i].begin()  
                        _d_prev_channeled_tpa_voltages[i].end() 
                        _d_channeled_tpa_voltages[i].begin()); 
        }
        // then add the input data
        thrust::copy(ftpa_voltages_in.begin() + i * _block_length_tpa, 
                     ftpa_voltages_in.begin() + (i + 1) * _block_length_tpa, 
                     _d_channeled_tpa_voltages[i].begin() + kernel_length_tpa);

        // update the overlap for the next hoard
        thrust::copy(ftpa_voltages_in.begin() + (i+1) * _block_length_tpa - kernel_length_tpa, 
                     ftpa_voltages_in.begin() + (i+1) * _block_length_tpa, 
                     _d_prev_channeled_tpa_voltages[i].begin());
    }

}

DeviceVoltageType const& BufferedDispenser::dispense(std::size_t chan_idx) const { // implements overlapped buffering of data

   return _d_channeled_tpa_voltages[chan_idx];
    
}





// void BufferedDispenser::dispense(std::size_t chan_idx, DeviceVoltageType& tpa_voltages_out) { // implements overlapped buffering of data

//     std::size_t offset = _block_length_tpa * chan_idx; // offset to the channel index

//     // copy Kernel length size of previous data to the next buffer
//     thrust::copy(_d_prev_ftpa_voltages.begin() + offset, 
//                     _d_prev_ftpa_voltages.begin() + offset + kernel_length_tpa, 
//                     tpa_voltages_out.begin()); // copy last  Kernel Length size of data to the next buffer

//     // copy current input data to the next buffer
//     thrust::copy(ftpa_voltages_in.begin(), 
//                     ftpa_voltages_in.end(), 
//                     tpa_voltages_out.begin() 
//                         + kernel_length_tpa); // from offset -> end 

//     // copy the last Kernel length size of data to the previous buffer for the next iteration
//     thrust::copy(ftpa_voltages_in.end() - kernel_length_tpa,
//                  ftpa_voltages_in.end(), 
//                  _d_prev_ftpa_voltages.begin() + this->_block_length_tpa * chan_idx); // copy the data to the previous data buffer.

// }


