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

    //Get first model from file
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
    return ((epoch >= _header.start_epoch) 
            && (epoch <= _header.end_epoch));
}

void DelayManager::safe_read(char* buffer, std::size_t nbytes)
{   
    if (_input_stream.eof())
    {
        // Reached the end of the delay model file
        // TODO: Decide what the behaviour should be here.
        throw std::runtime_error("Reached end of delay model file")
    }
    std::size_t nbytes_read = _input_stream.read(buffer, nbytes);
    if (nbytes_read < nbytes)
    {
        throw std::runtime_error("Could not complete read from delay file")
    }
}

void DelayManager::read_next_model()
{  
    std::size_t nbytes = 0;
    
    // Read the model header
    safe_read(reinterpret_cast<char*>(&_header), sizeof(_header));

    // Resize the arrays for the delay model on the host and GPU
    const std::size_t nelements = _header.nantennas * _header.nbeams;
    _delays_h.resize(nelements);
    _delays_d.resize(nelements);

    // Read the weight, offset, rate tuples from the file
     safe_read(reinterpret_cast<char*>(_delays_h.ptr()), 
       nelements * sizeof(DelayVectorHType::data_type));
}

} //namespace skyweaver