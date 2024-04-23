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

/**
 * @brief A struct wrapping the header of each delay model in the delay file
 * 
 */
struct DelayModelHeader
{
    uint32_t nantennas; // Number of antennas in the model set
    uint32_t nbeams; // Number of beams in the model set
    double start_epoch; // The start of the validity of the model as a unix epoch
    double end_epoch; // The end of the validity of the model as a unix epoch
};

/**
 * @brief      Class for managing loading of delay models from
 *             file and provisioning of them on the GPU.
 */
class DelayManager
{
public:
    
    typedef thrust::device_vector<float3> DelayVectorHType;
    typedef thrust::host_vector<float3> DelayVectorDType;

public:
    /**
     * @brief Construct a new Delay Manager object
     * 
     * @param delay_file A file containing delay models in skyweaver format
     * @param stream A cuda stream on which to execute host to device copies
     */
    DelayManager(std::string delay_file, cudaStream_t stream);
    ~DelayManager();
    DelayManager(DelayManager const&) = delete;

    /**
     * @brief Get the delay model for the given epoch
     * 
     * @param epoch A Unix epoch for which to fetch delays.
     *
     * @details Implemented for strictly increasing epochs. 
     *
     * @return A device vector containing the current delays
     */
    DelayVectorType const& delays(double epoch);

private:
    bool validate_model(double epoch) const;
    void read_next_model();
    void safe_read(char* buffer, std::size_t nbytes);

    cudaStream_t _copy_stream;
    DelayModelHeader _header;
    std::ifstream _input_stream;
    DelayVectorHType _delays_h;
    DelayVectorDType _delays_d;
};

} //namespace skyweaver

#endif // SKYWEAVER_DELAYMANAGER_HPP