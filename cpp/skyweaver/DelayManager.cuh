#ifndef SKYWEAVER_DELAYMANAGER_HPP
#define SKYWEAVER_DELAYMANAGER_HPP

#include "skyweaver/DelayModel.cuh"
#include <thrust/device_vector.h>
#include <iostream>
#include <format>
#include <exception>
#include <semaphore.h>

namespace skyweaver {

class InvalidDelayEpoch : public std::exception {
    private:
        double epoch;

    public:
        InvalidDelayEpoch(double _epoch) : epoch(_epoch) {}
        
        char* what () {
            return std::format("No valid delay solution for unix epoch {:f}\n", epoch);
        }
};

struct DelayFileHeader
{
    uint32_t version;
    uint32_t nantennas;
    uint32_t nbeams;
};

struct DelayModelHeader
{
    double start_epoch;
    double end_epoch;
};

/**
 * @brief      Class for managing loading of delay models from
 *             file and provisioning of them on the GPU.
 */
class DelayManager
{
public:
    
    typedef thrust::device_vector<float2> DelayVectorHType;
    typedef thrust::host_vector<float2> DelayVectorDType;

public:
    /**
     * @brief      Create a new DelayManager object
     *
     * @detail     The passed pipeline configuration contains the names
     *             of the POSIX shm and sem to connect to for the delay
     *             models.
     */
    DelayManager(std::string delay_file, cudaStream_t stream);
    ~DelayManager();
    DelayManager(DelayManager const&) = delete;

    /**
     * @brief      Get the delay model for the given epoch
     *
     * @detail     Implemented for strictly increasing epochs. 
     *
     * @return     A device vector containing the current delays
     */
    DelayVectorType const& delays(double epoch);

private:
    bool validate_model(double epoch) const;
    void read_next_model();
    
    cudaStream_t _copy_stream;
    DelayFileHeader _delay_file_header;
    DelayModelHeader _delay_model_header;
    std::ifstream _input_stream;
    DelayVectorHType _delays_h;
    DelayVectorDType _delays_d;
};

} //namespace skyweaver

#endif // SKYWEAVER_DELAYMANAGER_HPP