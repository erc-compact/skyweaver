#ifndef SKYWEAVER_INCOHERENTDEDISPERSER_CUH
#define SKYWEAVER_INCOHERENTDEDISPERSER_CUH

#include "skyweaver/PipelineConfig.hpp"

namespace skyweaver
{

//TODO: Document this interface

class IncoherentDedisperser
{
  public:
    IncoherentDedisperser(PipelineConfig const& config,
                          std::vector<float> const& dms,
                          std::size_t tscrunch = 1);
    IncoherentDedisperser(PipelineConfig const& config,
                          float dm,
                          std::size_t tscrunch = 1);
    ~IncoherentDedisperser();
    IncoherentDedisperser(IncoherentDedisperser const&)            = delete;
    IncoherentDedisperser& operator=(IncoherentDedisperser const&) = delete;

    template <typename InputVectorType, typename OutputVectorType>
    void dedisperse(InputVectorType const& tfb_powers,
                    std::vector<OutputVectorType>& tdb_powers);
    std::vector<int> const& delays() const;
    int max_sample_delay() const;

  private:
    void prepare();

    PipelineConfig const& _config;
    std::vector<float> _dms;
    std::size_t _tscrunch;
    int _max_sample_delay;
    float _scale_factor;
    std::vector<int> _delays;
};

} // namespace skyweaver

#endif // SKYWEAVER_INCOHERENTDEDISPERSER_CUH
