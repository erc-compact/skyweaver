#ifndef SKYWEAVER_STATISTICSCALCULATOR_CUH
#define SKYWEAVER_STATISTICSCALCULATOR_CUH

#include "skyweaver/DescribedVector.hpp"
#include "skyweaver/PipelineConfig.hpp"

#include <fstream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

namespace skyweaver
{

struct Statistics {
    float mean     = 0.0f;
    float std      = 0.0f;
    float skew     = 0.0f;
    float kurtosis = 0.0f;
};

inline std::ostream& operator<<(std::ostream& stream, Statistics const& val)
{
    stream << "mean: " << val.mean << ", std: " << val.std
           << ", skew: " << val.skew << ", kurtosis: " << val.kurtosis;
    return stream;
}

struct StatisticsFileHeader {
    uint32_t nantennas;
    uint32_t nchannels;
    uint32_t npol;
    uint64_t timestamp;
    uint32_t naccumulate;
};

/**
 * @brief      Class for wrapping statistics calculations on input data.
 */
class StatisticsCalculator
{
  public:
    typedef float ScalingType;
    typedef thrust::device_vector<ScalingType> ScalingVectorDType;
    typedef thrust::host_vector<ScalingType> ScalingVectorHType;
    typedef FPAStatsD<Statistics> StatisticsVectorDType;
    typedef FPAStatsH<Statistics> StatisticsVectorHType;

  public:
    /**
     * @brief      Create a new StatisticsCalculator object
     *
     * @param      config  The pipeline configuration.
     *
     * @detail     The passed pipeline configuration contains the names
     *             of the sem to connect to for the channel statistics
     */
    StatisticsCalculator(PipelineConfig const& config, cudaStream_t stream);
    ~StatisticsCalculator();
    StatisticsCalculator(StatisticsCalculator const&) = delete;

    /**
     * @brief      Calculate all statistics for the given input data
     */
    void calculate_statistics(FTPAVoltagesD<char2> const& ftpa_voltages);

    /**
     * @brief      Return the current channel input levels on GPU memory
     */
    StatisticsVectorDType const& statistics() const;

    /**
     * @brief      Return the current coherent beam offsets on GPU memory
     */
    ScalingVectorDType const& cb_offsets() const;

    /**
     * @brief      Return the current coherent beam scaling on GPU memory
     */
    ScalingVectorDType const& cb_scaling() const;

    /**
     * @brief      Return the current incoherent beam offsets on GPU memory
     */
    ScalingVectorDType const& ib_offsets() const;

    /**
     * @brief      Return the current incoherent beam scaling on GPU memory
     */
    ScalingVectorDType const& ib_scaling() const;

    /**
     * @brief Update the scaling arrays based on the last statistics
     *        calculation.
     */
    void update_scalings(ScalingVectorHType const& beamset_weights,
                         int nbeamsets);

  private:
    void dump_all_scalings() const;

    void dump_scalings(std::string const& timestamp,
                       std::string const& tag,
                       std::string const& path,
                       thrust::host_vector<float> const& ar) const;

  private:
    PipelineConfig const& _config;
    cudaStream_t _stream;
    StatisticsVectorDType _stats_d;
    StatisticsVectorHType _stats_h;
    ScalingVectorDType _cb_offsets_d;
    ScalingVectorHType _cb_offsets_h;
    ScalingVectorDType _cb_scaling_d;
    ScalingVectorHType _cb_scaling_h;
    ScalingVectorDType _ib_offsets_d;
    ScalingVectorHType _ib_offsets_h;
    ScalingVectorDType _ib_scaling_d;
    ScalingVectorHType _ib_scaling_h;
    std::ofstream _stats_file;
};

} // namespace skyweaver

#endif // SKYWEAVER_STATISTICSCALCULATOR_CUH
