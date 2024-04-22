#include "skyweaver/BeamformerPipeline.cuh"

#include "psrdada_cpp/Header.hpp"
#include "psrdada_cpp/cuda_utils.hpp"
#include "ascii_header.h"
#include "cuda.h"
#include <stdexcept>
#include <exception>
#include <cstdlib>

#define FBFUSE_SAMPLE_CLOCK_START_KEY "SAMPLE_CLOCK_START"
#define FBFUSE_SAMPLE_CLOCK_KEY "SAMPLE_CLOCK"
#define FBFUSE_SYNC_TIME_KEY "SYNC_TIME"


/*
According to dmsmear for COMPACT we need something like a 8192 sample kernel for 
filter. This is even smaller for BVRUSE due to the high frequency.

For memory purposes we assume we are working on a 64 channel subband
and that there are 64 antennas in the array. 

---- Pseudocode pipeline ----

data_loader = DataLoader(bunch of files)

input_buffer = OverlapSaveBuffer(overlap)

output_buffer = DedispersingOutputBuffer(???)

dedisperser = Dedisperser(???)

data_loader.seek(to time that I care about);
while data_loader.tell() < time that I want to end at:
    
    block = data_loader.load_block(nbytes)

    input_buffer.write(block) # includes copy to GPU

    block_to_process = input_buffer.read() # this is memory on the GPU

    # Transpose TAFTP to ATFP or AFPT whichever dspsr prefers

    >>>> Minimum Memory: 64 * 64 * 2 * 2 * 8192 * (1+4) = ~0.6-0.7 GB 
    >>>> assumes that we need a char and float copy in memory at the same time

    # get delays for timestamp if new set needed

    # make weights from delays always

    # recalculate scalings is new set needed (see delay model flags)

    for dm in dms:
        dedisperser.set_dm(dm)

        for ii, antenna in enumerate(antennas):
            dedispersed_antenna = dedisperser.dedisperse(block_to_process[ii])

        >>>> Minimum Memory: 64 * 64 * 2 * 2 * 8192 * 4 = 0.5 GB
        >>>> This is in addition to the input buffers so we are at 1 GB
        >>>> now 

        # Transpose AFPT (or whatever) to beamformer order FTPA
        >>>> Minimum Memory: 64 * 64 * 2 * 2 * 8192 * 1 = 0.1 GB
        >>>> We need to go back to 8 bit here otherwise the beamformer
        >>>> will not work.
        >>>> At 1.1 GB total now.

        # Form beams
        >>>>> Minimum Memory: 800 * 64 * 1 * 8192/4 = 100 MB
        >>>>> At this point we are using < 2 GB of memory
        
        output_buffer.write(beams) 

        # in the output buffer we dedisperse somehow
        # for a worst case DM of 1000 we have a delay of almost a second
        # this would require 4 times the length of data for the lowest frequency
        # band. Also if we want to do this on the GPU then all DM trials need
        # to be on GPU memory (is this actually true, with good overlap this 
        # doesn't matter and we can use dedisp). With all trials and the extra 
        # time it is 160 times the memory requirement here. Something like 16 GB.

        [NOTE] Dedisp doesn't work with 8 bit data

        # copy beams back to host 

        # Alternative option incoherent dedisperse on host with FDMT or whatever.


-- To be implemented
 - Data loader -> vvk
 - Input buffer + GPU copy -> vvk
 - Dedispersion code (waiting for input from Willem)
 - Delay and weight management -> eb
 - Updated coherent beamformer -> eb
 - Updated 32-bit incoherent beamformer -> eb



    
    
    











*/



namespace skyweaver {

template<typename CBHandler, typename IBHandler>
OfflinePipeline<CBHandler, IBHandler>::OfflinePipeline(PipelineConfig const& config,
    CBHandler& cb_handler,
    IBHandler& ib_handler,
    std::size_t input_data_buffer_size)
    : _config(config)
    , _sample_clock_start(0)
    , _sample_clock(0.0)
    , _sync_time(0.0)
    , _unix_timestamp(0.0)
    , _sample_clock_tick_per_block(0)
    , _call_count(0)
    , _cb_handler(cb_handler)
    , _ib_handler(ib_handler)
{
    BOOST_LOG_TRIVIAL(debug) << "Verifying all DADA buffer capacities";
    // Here we check the sizes of the block size intended to be passed to the 
    // pipeline. It is necessary to know this in advance to allow for allocation
    // of all intermediate buffers required for the beamformer.
    //
    // Input buffer checks:
    //
    std::size_t heap_group_size = (_config.total_nantennas() * _config.nchans()
        * _config.nsamples_per_heap() * _config.npol()) * sizeof(char2);
    
    //
    if (input_data_buffer_size % heap_group_size != 0)
    {
        throw std::runtime_error("Input DADA buffer is not a multiple "
            "of the expected heap group size");
    }
    _nheap_groups_per_block = input_data_buffer_size / heap_group_size;
    _nsamples_per_dada_block = _nheap_groups_per_block * _config.nsamples_per_heap();
    BOOST_LOG_TRIVIAL(debug) << "Number of heap groups per block: " << _nheap_groups_per_block;
    BOOST_LOG_TRIVIAL(debug) << "Number of samples/spectra per block: " << _nsamples_per_dada_block;
    if (_nsamples_per_dada_block % _config.cb_nsamples_per_block() != 0)
    {
        throw std::runtime_error("Input DADA buffer does not contain an integer "
            "multiple of the required number of samples per device block");
    }
    _taftp_db.resize(input_data_buffer_size / sizeof(char2), {0,0});

    /*
    In the original FBFUSE code, a series of output buffer checks were
    done here. These are now obsolete with the change to using handlers
    */
    std::size_t expected_cb_size = (_config.cb_nbeams() * _nsamples_per_dada_block
        / _config.cb_tscrunch() * _config.nchans() / _config.cb_fscrunch()) * sizeof(int8_t);
    _tbtf_db.resize(expected_cb_size, 0);

    std::size_t expected_ib_size = (_config.ib_nbeams() * _nsamples_per_dada_block
        / _config.ib_tscrunch() * _config.nchans() / _config.ib_fscrunch()) * sizeof(int8_t);
    _tf_db.resize(expected_ib_size, 0);

    // Calculate the timestamp step per block
    _sample_clock_tick_per_block = 2 * _config.total_nchans() * _nsamples_per_dada_block;
    BOOST_LOG_TRIVIAL(debug) << "Sample clock tick per block: " << _sample_clock_tick_per_block;

    //Default channel scales
    /*
    Note: these channel scalings are currently unused. The value is set to 1.0f here
    and in the voltage scaler this is used to multiply the input data. This could be
    left for future use as a channel mask. If there is no use case, it should be removed.
    */
    _channel_scalings.resize(_config.nchans(), 1.0f);

    BOOST_LOG_TRIVIAL(debug) << "Allocating CUDA streams";
    CUDA_ERROR_CHECK(cudaStreamCreate(&_h2d_copy_stream));
    CUDA_ERROR_CHECK(cudaStreamCreate(&_processing_stream));
    CUDA_ERROR_CHECK(cudaStreamCreate(&_d2h_copy_stream));

    BOOST_LOG_TRIVIAL(debug) << "Constructing delay and weights managers";
    _delay_manager.reset(new DelayManager(_config, _h2d_copy_stream));
    _gain_manager.reset(new GainManager(_config, _h2d_copy_stream));
    _weights_manager.reset(new WeightsManager(_config, _processing_stream));
    _stats_manager.reset(new ChannelScalingManager(_config, _processing_stream));
    _split_transpose.reset(new SplitTranspose(_config));
    _coherent_beamformer.reset(new CoherentBeamformer(_config));
    _incoherent_beamformer.reset(new IncoherentBeamformer(_config));
}

template<typename CBHandler, typename IBHandler
OfflinePipeline<CBHandler, IBHandler>::~OfflinePipeline()
{
    /*
    cb_handler.await();
    ib_hadler.await();
    */
    CUDA_ERROR_CHECK(cudaStreamDestroy(_h2d_copy_stream));
    CUDA_ERROR_CHECK(cudaStreamDestroy(_processing_stream));
    CUDA_ERROR_CHECK(cudaStreamDestroy(_d2h_copy_stream));
}

template<typename CBHandler, typename IBHandler
void OfflinePipeline<CBHandler, IBHandler>::set_header(RawBytes& header)
{
    Header parser(header);
    parser.purge();
    // There is a bug in DADA that results in keys made of subkeys not being writen if
    // the superkey is writen first. To get around this the order of key writes needs to
    // be carefully considered.
    parser.set<long double>(FBFUSE_SAMPLE_CLOCK_KEY, _sample_clock);
    parser.set<long double>(FBFUSE_SYNC_TIME_KEY, _sync_time);
    parser.set<std::size_t>(FBFUSE_SAMPLE_CLOCK_START_KEY, _sample_clock_start);
    header.used_bytes(header.total_bytes());
}

template<typename CBHandler, typename IBHandler
void OfflinePipeline<CBHandler, IBHandler>::init(RawBytes& header)
{
    BOOST_LOG_TRIVIAL(debug) << "Parsing DADA header";
    // Extract the time from the header and convert it to a double epoch
    Header parser(header);
    _sample_clock_start = parser.get<std::size_t>(FBFUSE_SAMPLE_CLOCK_START_KEY);
    _sample_clock = parser.get<long double>(FBFUSE_SAMPLE_CLOCK_KEY);
    _sync_time = parser.get<long double>(FBFUSE_SYNC_TIME_KEY);

    // Need to set the header information on the coherent beam output block
    std::vector<char> cb_header(4096);
    RawBytes cb_header_rb = RawBytes(cb_header.ptr(), cb_header.size());
    set_header(cb_header_rb);
    cb_handler.init(cb_header_rb);

    // Need to set the header information on the coherent beam output block
    std::vector<char> ib_header(4096);
    RawBytes ib_header_rb = RawBytes(ib_header.ptr(), ib_header.size());
    set_header(ib_header_rb);
    ib_handler.init(ib_header_rb);

}

template<typename CBHandler, typename IBHandler
void OfflinePipeline<CBHandler, IBHandler>::process(
    VoltageVectorType& taftp_vec,
    PowerVectorType& tbtf_vec, 
    PowerVectorType& tf_vec,
    )
{
    BOOST_LOG_TRIVIAL(debug) << "Executing coherent beamforming pipeline";

    /*
     * For the COMPACT project the gains are all 1.0f, as the data is 
     * recorded from PTUSE. As such, this code could be skipped or hidden
     * behind an if condition.
     */
    BOOST_LOG_TRIVIAL(debug) << "Checking for complex gain updates";
    auto const& gains = _gain_manager->gains();

    // Need to add the unix timestmap to the delay manager here
    // to fetch valid delays for this epoch.
    BOOST_LOG_TRIVIAL(debug) << "Checking for delay updates";
    auto const& delays = _delay_manager->delays(_unix_timestamp);

    // Stays the same
    BOOST_LOG_TRIVIAL(debug) << "Calculating weights at unix time: " << _unix_timestamp;
    auto const& weights = _weights_manager->weights(delays, _unix_timestamp, _delay_manager->epoch());

    // Stays the same
    BOOST_LOG_TRIVIAL(debug) << "Checking if channel statistics update request";
    _stats_manager->channel_statistics(taftp_vec);

    // This is no longer needed if no gains are being applied
    BOOST_LOG_TRIVIAL(debug) << "Applying complex gain corrections";
    voltage_scaling(taftp_vec, taftp_vec, gains, _channel_scalings, _processing_stream);
    
    // Stays the same    
    BOOST_LOG_TRIVIAL(debug) << "Transposing input data from TAFTP to FTPA order";
    _split_transpose->transpose(taftp_vec, _split_transpose_output, _processing_stream);

    // As the beamformer will now output the CB-IB product
    // it necessary to calculate new scaling factors.
    BOOST_LOG_TRIVIAL(debug) << "Forming coherent beams";

    // Not yet implemented
    auto const& cb_minus_ib_scaling = _stats_manager->cb_minus_ib_scaling();
    auto const& cb_minus_ib_offsets = _stats_manager->cb_minus_ib_offsets();

    // Need a new beamformer call to do cb-ib with Stokes selection
    _coherent_beamformer->beamform(
        _split_transpose_output, weights,
        cb_minus_ib_scaling, cb_minus_ib_offsets,
        tbtf_vec, _processing_stream);

    BOOST_LOG_TRIVIAL(debug) << "Forming incoherent beam";
    auto const& ib_scaling = _stats_manager->ib_scaling();
    auto const& ib_offsets = _stats_manager->ib_offsets();
    _incoherent_beamformer->beamform(
            taftp_vec, tf_vec, ib_scaling, ib_offsets,
             _processing_stream);
}

template<typename CBHandler, typename IBHandler
bool OfflinePipeline<CBHandler, IBHandler>::operator()(RawBytes& data)
{
    // This update goes at the top of the method to ensure
    // that it is always incremented on each pass
    ++_call_count;

    BOOST_LOG_TRIVIAL(debug) << "Pipeline::operator() called (count = " << _call_count << ")";
    // We first need to synchronize the h2d copy stream to ensure that
    // last host to device copy has completed successfully. When this is
    // done we are free to call swap on the double buffer without affecting
    // any previous copy.
    CUDA_ERROR_CHECK(cudaStreamSynchronize(_h2d_copy_stream));
    _taftp_db.swap();

    if (data.used_bytes() != _taftp_db.a().size()*sizeof(char2))
    {
        throw std::runtime_error(std::string("Unexpected buffer size, expected ")
            + std::to_string(_taftp_db.a().size()*sizeof(char2))
            + " but got "
            + std::to_string(data.used_bytes()));
    }
    CUDA_ERROR_CHECK(cudaMemcpyAsync(static_cast<void*>(_taftp_db.a_ptr()),
        static_cast<void*>(data.ptr()), data.used_bytes(),
        cudaMemcpyHostToDevice, _h2d_copy_stream));

    // If we are on the first call we can exit here as there is no
    // data on the GPU yet to process.
    if (_call_count == 1)
    {
        return false;
    }
    // Here we block on the processing stream before swapping
    // the processing buffers
    CUDA_ERROR_CHECK(cudaStreamSynchronize(_processing_stream));
    _tbtf_db.swap();
    _tf_db.swap();
    // Calculate the unix timestamp for the block that is about to be processed
    // (which is the block passed the last time that operator() was called)
    _unix_timestamp = (_sync_time + (_sample_clock_start +
        ((_call_count - 2) * _sample_clock_tick_per_block))
    / _sample_clock);
    process(_taftp_db.b(), _tbtf_db.a(), _tf_db.a());

    // If we are on the second call we can exit here as there is not data
    // that has completed processing at this stage.
    if (_call_count == 2)
    {
        return false;
    }

    CUDA_ERROR_CHECK(cudaStreamSynchronize(_d2h_copy_stream));

    // Need to check that the sizes here are all in bytes
    RawBytes tbtf_db_rb = RawBytes(_tbtf_db.b_ptr(), _tbtf_db.size())
    tbtf_db_rb.used_bytes(_tbtf_db.size());
    RawBytes tf_db_rb = RawBytes(_tf_db.b_ptr(), _tf_db.size())
    tbtf_db_rb.used_bytes(_tf_db.size());

    // Currently these handler calls are going to result in synchronous
    // Down stream calls which will pause execution on the GPU
    // Asnync handlers are required to avoid this (or semi-async where
    // the synchronous part is limited to a H2H memcopy).
    cb_handler(tbtf_db_rb);
    ib_handler(tf_db_rb);
    return false;
}

} //namespace skyweaver