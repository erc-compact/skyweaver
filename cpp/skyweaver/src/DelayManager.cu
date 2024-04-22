#include "skyweaver/DelayManager.cuh"
#include "psrdada_cpp/cuda_utils.hpp"
#include <errno.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cerrno>

namespace skyweaver {

DelayManager::DelayManager(std::string delay_file, cudaStream_t stream)
    : _copy_stream(stream)
{
    BOOST_LOG_TRIVIAL(debug) << "Constructing new DelayManager instance";
    BOOST_LOG_TRIVIAL(debug) << std::format("Opening delay model file: {}", delay_file);
    _input_stream.open(delay_file, ios::in | ios::binary);
    if (!_input_stream)
    {
        throw std::runtime_error(
            std::format("Error: Unable to open file {} due to error: {}",
                delay_file, std::strerror(errno)
            );
        )
    }
    BOOST_LOG_TRIVIAL(debug) << "Delay model file successfully opened";

    // Read the file header
    _input_stream.read(reinterpret_cast<char*>(&_delay_file_header), sizeof(_delay_file_header));

    // Read the scalar antenna weights from file
    antenna_weights.resize(_delay_file_header.nantennas);
    _input_stream.read(reinterpret_cast<char*>(
        antenna_weights.ptr()),
         _delay_file_header.nantennas * sizeof(float));

    // Resize the arrays for the delay model
    _delays_h.resize(_delay_file_header.nantennas * _delay_file_header.nbeams);
    _delays_d.resize(_delay_file_header.nantennas * _delay_file_header.nbeams);
    read_next_model();
}

DelayManager::~DelayManager()
{
    if (_input_stream.is_open()) {
        _input_stream.close();
    }   

    if (_input_stream.fail()) {
            throw std::runtime_error(
                std::format("Error: Unable to close delay file due to error: {}",
                    std::strerror(errno)
            );
        )
    }
}

DelayManager::DelayVectorType const& DelayManager::delays(double epoch)
{
    // This function should return the delays in GPU memory
    
    // Scan through the model file until we reach model that 
    // contains valid delays for the given epoch or until we 
    // hit EOF (which throws an exception).
    while (!validate_model(epoch))
    {
        read_next_model();
    }
    thrust::copy(_delays_h.begin(), _delays_h.end(), _delays_d.begin());
    return _delays_d;
}

bool DelayManager::validate_model(double epoch) const
{
    return ((epoch >= _delay_model_header.start_epoch) 
            && (epoch <= _delay_model_header.end_epoch));
}

void DelayManager::read_next_model()
{  
    std::size_t nbytes = 0;
    if (_input_stream.eof())
    {
        // Reached the end of the delay model file
        // TODO: Decide what the behaviour should be here.
        throw std::runtime_error("Reached end of delay model file")
    }
    nbytes = _input_stream.read(
        reinterpret_cast<char*>(&_delay_model_header), 
        sizeof(_delay_model_header));
    if (nbytes != sizeof(_delay_model_header))
    {
        throw std::runtime_error("Unable to read delay model header from file");
    }
     nbytes = _input_stream.read(
        reinterpret_cast<char*>(_delays_h.ptr()), 
        _delays_h.size() * sizeof(float2));
    if (nbytes != sizeof(_delay_model_header))
    {
        throw std::runtime_error("Unable to read complete delay model from file");
    }
}

} //namespace skyweaver