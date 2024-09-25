#ifndef SKYWEAVER_PIPELINECONFIG_HPP
#define SKYWEAVER_PIPELINECONFIG_HPP

#include "psrdada_cpp/common.hpp"
#include "skyweaver/DedispersionPlan.hpp"
#include "skyweaver/skyweaver_constants.hpp"

#include <string>

namespace skyweaver
{
  struct WaitConfig
  {
    int iterations;
    int sleep_time;
    std::size_t min_free_space;
  };

  struct PreWriteConfig
  {
    bool is_enabled;
    WaitConfig wait;
  };


/**
 * @brief      Class for wrapping the skyweaver pipeline configuration.
 */
class PipelineConfig
{
  public:
    PipelineConfig();
    ~PipelineConfig();
    PipelineConfig(PipelineConfig const&) = delete;

    /**
     * @brief      Get the path to the file containing delay solutions
     */
    std::string const& delay_file() const;

    /**
     * @brief      Set the path to the file containing delay solutions
     */
    void delay_file(std::string const&);

    /**
     * @brief      Get the list of input files
     */
    std::vector<std::string> const& input_files() const;

    /**
     * @brief      Set list of input files
     */
    void input_files(std::vector<std::string> const&);

    /**
     * @brief      Check if the input files are contiguous
     */
    bool check_input_contiguity() const;

    /**
     * @brief      Set whether the input files are contiguous
     */
    void check_input_contiguity(bool);

    /**
     * @brief      Get the size of the DADA header in the input files
     */
    std::size_t dada_header_size() const;

    /**
     * @brief      Set the size of the DADA header in the input files
     */
    void dada_header_size(std::size_t);

    /**
     * @brief      Set the list of input files from a text file
     *
     * @details    File format is a newline separated list of
     *             absolute or relative filepaths. Lines beginning
     *             with # are considered to be comments
     */
    void read_input_file_list(std::string filename);

    /**
     * @brief      Get the directory path for output files
     */
    std::string const& output_dir() const;

    /**
     * @brief      Set the directory path for output files
     */
    void output_dir(std::string const&);

    /**
     * @brief Get the prefix for the output files
     */
    std::string const& output_file_prefix() const;

    /**
     * @brief Set the prefix for the output files
     */
    void output_file_prefix(std::string const&);

    /**
     * @brief Get the maximum size of the output files
     */
    std::size_t max_output_filesize() const;

    /**
     * @brief Set the maximum output file size
     */
    void max_output_filesize(std::size_t);

    /**
     * @brief      Get the centre frequency for the subband to
     *             be processed by this instance.
     */
    double centre_frequency() const;

    /**
     * @brief      Set the centre frequency for the subband to
     *             be processed by this instance.
     *             Units of Hz.
     */
    void centre_frequency(double cfreq);

    /**
     * @brief      Get the bandwidth of the subband to
     *             be processed by this instance.
     *             Units of Hz.
     *
     *
     */
    double bandwidth() const;

    /**
     * @brief      Set the bandwidth of the subband to
     *             be processed by this instance.
     *             Units of Hz.
     */
    void bandwidth(double bw);

    /**
     * @brief      Return the centre frequency of each channel in the
     *             subband to be processed.
     *             Units of Hz.
     *
     */
    std::vector<double> const& channel_frequencies() const;

    /**
     * @brief      Get the coherent dm trials to be dedispersed to
     *
     * @details    Helper method to return the coherent DMs from the ddplan
     */
    std::vector<float> const& coherent_dms() const;

    /**
     * @brief      Get a reference to the dedispersion plan
     */
    DedispersionPlan& ddplan();

    /**
     * @brief      Get a const reference to the dedispersion plan
     */
    DedispersionPlan const& ddplan() const;

    /**
     * @brief      configures wait for filesystem space
     */
    void configure_wait(std::string argument);

    /**
     * @brief      Enable/disable incoherent dedispersion based fscrunch after
     * beamforming
     */
    void enable_incoherent_dedispersion(bool enable);

    /**
     * @brief      Check if incoherent dedispersion based fscrunch after
     * beamforming is enabled
     */
    bool enable_incoherent_dedispersion() const;

    /**
     * @brief      Return the number of time samples to be integrated
     *             in the coherent beamformer.
     */
    std::size_t cb_tscrunch() const { return SKYWEAVER_CB_TSCRUNCH; }

    /**
     * @brief      Return the number of frequency channels to be integrated
     *             in the coherent beamformer.
     */
    std::size_t cb_fscrunch() const { return SKYWEAVER_CB_FSCRUNCH; }

    /**
     * @brief      Return the number of time samples to be integrated
     *             in the incoherent beamformer.
     */
    std::size_t ib_tscrunch() const { return SKYWEAVER_IB_TSCRUNCH; }

    /**
     * @brief      Return the number of frequency channels to be integrated
     *             in the incoherent beamformer.
     */
    std::size_t ib_fscrunch() const { return SKYWEAVER_IB_FSCRUNCH; }

    /**
     * @brief     Get the setting for CB-IB subtraction
     */
    bool cb_ib_subtract() const { return SKYWEAVER_IB_SUBTRACTION; }

    /**
     * @brief      Return the number of beams to be formed by the coherent
     * beamformer
     */
    std::size_t nbeams() const { return SKYWEAVER_NBEAMS; }

    /**
     * @brief      Return the number of samples that will be processed
     *             in each batch.
     */
    std::size_t nsamples_per_block() const
    {
        return SKYWEAVER_CB_NSAMPLES_PER_BLOCK;
    }

    PreWriteConfig pre_write_config() const
    {
        return _pre_write_config;
    }

    /**
     * @brief Return the total number of samples to read from file in each gulp.
     *
     * @details Must be a multiple of nsamps per heap and greater than\
     *          the dedispersion kernel size.
     */
    std::size_t gulp_length_samps() const;

    /**
     * @brief Set the total number of samples to read from file in each gulp.
     *
     * @details Must be a multiple of nsamps per heap and greater than\
     *          the dedispersion kernel size.
     */
    void gulp_length_samps(std::size_t);

    float start_time() const;
    void start_time(float);

    float duration() const;
    void duration(float);

    /**
     * Below are methods to get and set the power scaling and offset in the
     * beamformer these are tricky parameters that need to be set correctly for
     * the beamformer to function as expected. The values are used when
     * downcasting from floating point to 8-bit integer at the end stage of
     * beamforming. The scaling is the last step in the code and as such the
     * factors can be quite large.
     *
     * The scaling and offset are applied such that:
     *
     *    int8_t scaled_power = static_cast<int8_t>((power - offset) / scaling);
     *
     * In the case above, the variable power is the power after summing all
     * antennas, timesamples and frequency channels requested (tscrunch and
     * fscrunch, respectively). The offset and scaling can be estimated if the
     * power per input F-engine stream is known. The interested reader can
     * examine the PipelineConfig::update_power_offsets_and_scalings method to
     * see how the exact scaling and offsets are calculated for the coherent and
     * incoherent beamformer
     *
     * Note: We do not assume different scaling per channel, if there are
     * significantly different power levels in each channel the scaling should
     * always be set to accommodate the worst cast scenario.
     */

    /**
     * @brief      Set the output standard deviation for data out
     *             of both the coherent and incoherent beamformers
     */
    void output_level(float level);

    /**
     * @brief      Get the output standard deviation for data out
     *             of both the coherent and incoherent beamformers
     */
    float output_level() const;

    /**
     * @brief      Return the total number of antennas that will be beamformed
     */
    std::size_t nantennas() const { return SKYWEAVER_NANTENNAS; }

    /**
     * @brief      Return the total number of frequency channels in the input
     * data
     */
    std::size_t nchans() const { return SKYWEAVER_NCHANS; }

    /**
     * @brief      Return the F-engine channelisation mode
     * data
     */
    std::size_t total_nchans() const;

    /**
     * @brief      Set the F-engine channelisation mode
     * data
     */
    void total_nchans(std::size_t);

    /**
     * @brief      Return the number of polarisations in the observation
     *
     * @note       This better be 2 otherwise who knows what will happen...
     */
    std::size_t npol() const { return SKYWEAVER_NPOL; }

    /**
     * @brief      Return the Stokes mode
     */
    std::string stokes_mode() const { return _stokes_mode; }

    /**
     * @brief      Set the Stokes mode
     */
    void stokes_mode(std::string const& stokes_mode_)
    {
        _stokes_mode = stokes_mode_;
    }

    /**
     * @brief      Return the number of time samples per F-engine SPEAD heap.
     *
     * @note       This corresponds to the inner "T" dimension in the input
     *             TAF[T]P order data.
     */
    std::size_t nsamples_per_heap() const
    {
        return SKYWEAVER_NSAMPLES_PER_HEAP;
    }

  private:
    std::size_t convertMemorySize(const std::string& str) const;
    void calculate_channel_frequencies() const;
    void update_power_offsets_and_scalings();

  private:
    /* input file options*/
    std::vector<std::string> _input_files;
    bool _check_input_contiguity;
    std::size_t _dada_header_size;

    /* Other required files*/
    std::string _delay_file;

    std::string _output_dir;
    std::size_t _max_output_filesize;
    std::string _output_file_prefix;
    bool _enable_incoherent_dedispersion;
    double _cfreq;
    double _bw;
    mutable bool _channel_frequencies_stale;
    std::size_t _gulp_length_samps;
    float _start_time;
    float _duration;
    std::size_t _total_nchans;
    std::string _stokes_mode;
    float _output_level;
    DedispersionPlan _ddplan;
    mutable std::vector<double> _channel_frequencies;
    PreWriteConfig _pre_write_config;
};

} // namespace skyweaver

#endif // SKYWEAVER_PIPELINECONFIG_HPP
