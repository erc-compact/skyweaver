#ifndef SKYWEAVER_DESCRIBEDVECTORS_HPP
#define SKYWEAVER_DESCRIBEDVECTORS_HPP

#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <thrust/mr/allocator.h>
#include <thrust/mr/universal_memory_resource.h>
#include <type_traits>
#include <vector>

namespace skyweaver
{

template <typename T, typename A>
inline std::ostream& operator<<(std::ostream& stream,
                                std::vector<T, A> const& vec)
{
    bool first = true;
    stream << "(";
    for(T const& val: vec) {
        if(!first) {
            stream << ", ";
        } else {
            first = false;
        }
        stream << val;
    }
    stream << ")";
    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, char2 const& val)
{
    stream << "(" << static_cast<int>(val.x) << "," << static_cast<int>(val.y)
           << ")";
    return stream;
}

// Define the Dimension enum
enum Dimension {
    TimeDim,
    FreqDim,
    BeamDim,
    AntennaDim,
    PolnDim,
    DispersionDim
};

std::string dimensions_to_string(std::vector<Dimension> const& dims);

// Helper to check if a Dimension is in a list of Dimensions
template <Dimension D, Dimension... Ds>
struct contains;

template <Dimension D>
struct contains<D>: std::false_type {
};

template <Dimension D, Dimension Head, Dimension... Tail>
struct contains<D, Head, Tail...>
    : std::conditional_t<D == Head, std::true_type, contains<D, Tail...>> {
};

template <typename VectorType_, Dimension... dims>
struct DescribedVector;

template <typename T, typename U>
struct is_dv_copyable: std::false_type {
};

template <template <typename, typename> class Container1,
          template <typename, typename>
          class Container2,
          typename T,
          typename A1,
          typename A2,
          Dimension... dims>
struct is_dv_copyable<DescribedVector<Container1<T, A1>, dims...>,
                      DescribedVector<Container2<T, A2>, dims...>>
    : std::true_type {
};

/**
 * @brief A class for holding vectors of specific dimension with associated
 * metadata
 *
 * @tparam VectorType_ The container type used for storing the data
 * @tparam dims        The dimension types of the data
 */
template <typename VectorType_, Dimension... dims>
struct DescribedVector {
    using VectorType             = VectorType_;
    using value_type             = typename VectorType::value_type;
    using FrequenciesType        = std::vector<double>;
    using DispersionMeasuresType = std::vector<float>;

    /**
     * @brief Construct a new Described Vector object
     *
     */
    DescribedVector()
        : _dms_stale(true), _frequencies_stale(true), _dims{dims...},
          _tsamp(0.0)
    {
    }

    /**
     * @brief Construct a new Described Vector object of specific size
     *
     * @param sizes The sizes of the dimensions (must match the number of
     * dimensions)
     */
    DescribedVector(std::initializer_list<std::size_t> sizes)
        : _dms_stale(true), _frequencies_stale(true), _sizes(sizes),
          _dims{dims...}, _tsamp(0.0)
    {
        if(_sizes.size() != sizeof...(dims)) {
            throw std::invalid_argument(
                "Number of sizes must match number of dimensions");
        }
        _vector.resize(calculate_nelements());
    }

    /**
     * @brief Construct a new Described Vector object of specific size
     *
     * @param sizes The sizes of the dimensions (must match the number of
     * dimensions)
     */
    DescribedVector(std::initializer_list<std::size_t> sizes, value_type default_value)
        : _dms_stale(true), _frequencies_stale(true), _sizes(sizes),
          _dims{dims...}, _tsamp(0.0)
    {
        if(_sizes.size() != sizeof...(dims)) {
            throw std::invalid_argument(
                "Number of sizes must match number of dimensions");
        }
        _vector.resize(calculate_nelements(), default_value);
    }

    /**
     * @brief Destroy the Described Vector object
     *
     */
    ~DescribedVector() {};

    /**
     * @brief Construct a new Described Vector object
     *
     * @tparam OtherDescribedVector  The type of the vector that is being copied
     * from
     * @param other                  The instance of the vector that is being
     * copied from
     *
     * @details The type of the other vector is checked at compile time to
     * assure it is copyable to the current vector. This allows for reuse of the
     * implicit H2D and D2H copy constructors of thrust vectors.
     */
    template <typename OtherDescribedVector,
              typename = std::enable_if_t<
                  is_dv_copyable<DescribedVector, OtherDescribedVector>::value>>
    DescribedVector(const OtherDescribedVector& other)
        : _sizes(other._sizes), _dims{other._dims},
          _vector(other._vector.begin(), other._vector.end()), _dms(other._dms),
          _frequencies(other._frequencies), _dms_stale(other._dms_stale),
          _tsamp(other._tsamp), _utc_offset(other._utc_offset),
          _reference_dm(other._reference_dm),
          _frequencies_stale(other._frequencies_stale)
    {
    }

    /**
     * @brief Resize the dimensions and set the metadata to match another vector
     *
     * @tparam OtherDescribedVector The type of the other vector
     * @param other                 The instance of the other vector
     */
    template <typename OtherDescribedVector,
              typename = std::enable_if_t<
                  is_dv_copyable<DescribedVector, OtherDescribedVector>::value>>
    void like(const OtherDescribedVector& other)
    {
        resize(other._sizes);
        metalike(other);
    }

    /**
     * @brief Set the metadata to match another vector
     *
     * @tparam OtherDescribedVector The type of the other vector
     * @param other                 The instance of the other vector
     */
    template <typename OtherDescribedVector>
    void metalike(const OtherDescribedVector& other)
    {
        _frequencies       = other._frequencies;
        _frequencies_stale = other._frequencies_stale;
        _dms               = other._dms;
        _dms_stale         = other._dms_stale;
        _tsamp             = other._tsamp;
        _utc_offset        = other._utc_offset;
        _reference_dm      = other._reference_dm;
    }

    /**
     * @brief Get the data at an index
     *
     * @param idx     The index to fetch
     * @return auto&
     *
     * @details Despite the dimensions being tracked, the indexing to the
     * underlying data is linear (flat indexing).
     */
    auto& operator[](std::size_t idx) { return _vector[idx]; }

    /**
     * @brief Get the data at an index
     *
     * @param idx     The index to fetch
     * @return auto&
     *
     * @details Despite the dimensions being tracked, the indexing to the
     * underlying data is linear (flat indexing).
     */
    auto const& operator[](std::size_t idx) const { return _vector[idx]; }


    // also the at() method
    auto& at(std::size_t idx) { return _vector[idx]; }
    auto const& at(std::size_t idx) const { return _vector[idx]; }
    

    /**
     * @brief Resize the dimensions of the vector
     *
     * @param sizes The extents of each dimension
     */
    void resize(std::initializer_list<std::size_t> sizes)
    {
        _sizes = sizes;
        _vector.resize(calculate_nelements());
    }

    /**
     * @brief Resize the dimensions of the vector
     *
     * @param sizes The extents of each dimension
     */
    void resize(std::vector<std::size_t> const& sizes)
    {
        _sizes = sizes;
        _vector.resize(calculate_nelements());
    }

    /**
     * @brief Get the size of the underlying vector
     *
     * @return std::size_t
     */
    std::size_t size() const { return _vector.size(); }

    /**
     * @brief Get a pointer-like object to the underlying data
     *
     * @return auto
     */
    auto data() noexcept { return _vector.data(); }

    /**
     * @brief Get a pointer-like object to the underlying data
     *
     * @return auto
     */
    const auto data() const noexcept { return _vector.data(); }

    /**
     * @brief Get an iterator that points to the start of the underlying data
     *
     * @return auto
     */
    auto begin() noexcept { return _vector.begin(); }

    /**
     * @brief Get a const iterator that points to the start of the underlying
     * data
     *
     * @return auto
     */
    const auto begin() const noexcept { return _vector.begin(); }

    /**
     * @brief Get an iterator that points to the end of the underlying data
     *
     * @return auto
     */
    auto end() noexcept { return _vector.end(); }

    /**
     * @brief Get a const iterator that points to the end of the underlying data
     *
     * @return auto
     */
    const auto end() const noexcept { return _vector.end(); }

    /**
     * @brief Get a const reference to the underlying container
     *
     * @return auto
     */
    VectorType const& vector() const { return _vector; }

    /**
     * @brief Get a mutable reference to the underlying container
     *
     * @return auto
     */
    VectorType& vector() { return _vector; }

    std::size_t calculate_nelements() const
    {
        return std::accumulate(_sizes.begin(),
                               _sizes.end(),
                               1,
                               std::multiplies<>());
    }

    template <Dimension dim>
    std::size_t get_dim_extent() const
    {
        if constexpr(!contains<dim, dims...>::value) {
            return 1;
        } else {
            std::size_t total_size = 1;
            for(std::size_t i = 0; i < _dims.size(); ++i) {
                if(_dims[i] == dim) {
                    total_size *= _sizes[i];
                }
            }
            return total_size;
        }
    }

    /**
     * @brief Get the frequencies associated with any channels
     *
     * @return FrequenciesType const&
     */
    FrequenciesType const& frequencies() const
    {
        if(_frequencies_stale) {
            throw std::runtime_error(
                "Frequencies must be explicitly set before they can be used");
        }
        return _frequencies;
    }

    /**
     * @brief Set the frequencies of the channels
     *
     * @param freqs
     */
    void frequencies(FrequenciesType const& freqs)
    {
        if(freqs.size() != get_dim_extent<FreqDim>()) {
            throw std::runtime_error("Invalid number of frequecies passed.");
        }
        _frequencies_stale = false;
        _frequencies       = freqs;
    }

    /**
     * @brief Set the frequencies of the channels
     *
     * @param freq
     */
    void frequencies(typename FrequenciesType::value_type const& freq)
    {
        if(get_dim_extent<FreqDim>() != 1) {
            throw std::runtime_error("Invalid number of frequecies passed.");
        }
        _frequencies_stale = false;
        _frequencies.resize(1, freq);
    }

    /**
     * @brief Return the number of frequency channels
     *
     * @return std::size_t
     *
     * @details Repeated dimensions of the same type will be aggregated
     */
    std::size_t nchannels() const { return get_dim_extent<FreqDim>(); }

    /**
     * @brief Return the number of antennas
     *
     * @return std::size_t
     *
     * @details Repeated dimensions of the same type will be aggregated
     */
    std::size_t nantennas() const { return get_dim_extent<AntennaDim>(); }

    /**
     * @brief Return the number of time samples
     *
     * @return std::size_t
     *
     * @details Repeated dimensions of the same type will be aggregated
     */
    std::size_t nsamples() const { return get_dim_extent<TimeDim>(); }

    /**
     * @brief Return the number of beams
     *
     * @return std::size_t
     *
     * @details Repeated dimensions of the same type will be aggregated
     */
    std::size_t nbeams() const { return get_dim_extent<BeamDim>(); }

    /**
     * @brief Return the number of polarisation dimensions
     *
     * @return std::size_t
     *
     * @details Repeated dimensions of the same type will be aggregated.
     *
     * @note The underlying data type may also encode the polarisation
     *       e.g. char4 for Stokes data. In this case the number of
     * polarisations reported by this method will be 1.
     */
    std::size_t npol() const { return get_dim_extent<PolnDim>(); }

    /**
     * @brief Return the number of DMs
     *
     * @return std::size_t
     */
    std::size_t ndms() const { return get_dim_extent<DispersionDim>(); }

    /**
     * @brief Set the time resolution of the data
     *
     * @param tsamp_ The time resolution in seconds
     */
    void tsamp(double tsamp_) { _tsamp = tsamp_; }

    /**
     * @brief Get the time resolution of the data
     *
     * @return double
     */
    double tsamp() const { return _tsamp; }

    /**
     * @brief Get the latency of this data w.r.t. the stream
     *
     * @return double
     *
     * @details This is intended to allow filter delays to be
     *          propagated through the code. e.g. the DM delay
     *          from incoherent dedispersion may be reflected
     *          in this parameter.
     */
    double utc_offset() const { return _utc_offset; }


    /**
     * @brief Set the latency of this data w.r.t. the stream
     *
     * @param offset The latency in seconds
     */
    void utc_offset(double offset) { _utc_offset = offset; }

    /**
     * @brief Get the sizes of each of the dimensions
     *
     * @return std::vector<std::size_t> const&
     */
    std::vector<std::size_t> const& extents() const { return _sizes; }

    /**
     * @brief Get the list of DMs
     *
     * @return DispersionMeasuresType const&
     */
    DispersionMeasuresType const& dms() const
    {
        if(_dms_stale) {
            throw std::runtime_error(
                "Dms must be explicitly set before they can be used");
        }
        return _dms;
    }

    /**
     * @brief Set the list of DMs
     *
     * @param dms
     */
    void dms(DispersionMeasuresType const& dms)
    {
        if(dms.size() != get_dim_extent<DispersionDim>()) {
            throw std::runtime_error(
                "Invalid number of dispersion measures passed.");
        }
        _dms_stale = false;
        _dms       = dms;
    }

    /**
     * @brief Set the DM of the data
     *
     * @param dm
     */
    void dms(typename DispersionMeasuresType::value_type const& dm)
    {
        if(get_dim_extent<DispersionDim>() != 1) {
            throw std::runtime_error(
                "Invalid number of dispersion measures passed.");
        }
        _dms_stale = false;
        _dms.resize(1);
        _dms[0] = dm;
    }

    /**
     * @brief Set the reference coherent DM of the channels
     *
     * @param dm The coherent DM to which the channels have been referenced
     */
    void reference_dm(float dm) { _reference_dm = dm; }

    /**
     * @brief Get the reference coherent DM of the channels
     *
     * @return float
     */
    float reference_dm() const { return _reference_dm; }

    /**
     * @brief Get the dimensions of the data as a string
     *
     * @return std::string
     *
     * @details The dimensions are ordered slowest to fastest (e.g. TAFTP)
     *          where T [time] is the slowest dimension and P [poln.] is
     *          fastest dimension.
     */
    std::string dims_as_string() const { return dimensions_to_string(_dims); }

    /**
     * @brief Return a string describing the vector in detail
     *
     * @return std::string
     */
    std::string describe() const
    {
        std::stringstream stream;
        stream << "DescribedVector: " << dimensions_to_string(_dims) << " ("
               << size() << " elements)\n";
        stream << "  dimensions: " << _sizes << "\n";
        stream << "  nsamples: " << nsamples() << "\n";
        stream << "  nchans: " << nchannels() << "\n";
        stream << "  nantennas: " << nantennas() << "\n";
        stream << "  nbeams: " << nbeams() << "\n";
        stream << "  ndms: " << ndms() << "\n";
        stream << "  npol: " << npol() << "\n";
        stream << "  ndms: " << ndms() << "\n";
        stream << std::setprecision(15);
        stream << "  frequencies (Hz): " << _frequencies << "\n";
        stream << "  DMs (pc cm^-3): " << _dms << "\n";
        stream << "  Coherent DM (pc cm^-3): " << reference_dm() << "\n";
        stream << "  Time resolution (s): " << _tsamp << "\n";
        stream << "  Time offset (s): " << _utc_offset << "\n";
        stream << std::setprecision(6);
        return stream.str();
    }

    DispersionMeasuresType _dms;
    bool _dms_stale;
    FrequenciesType _frequencies;
    bool _frequencies_stale;
    std::vector<std::size_t> _sizes;
    std::vector<Dimension> _dims;
    double _tsamp       = 0.0;
    double _utc_offset  = 0.0;
    float _reference_dm = 0.0;
    VectorType _vector;
};

template <typename T>
struct is_device_vector: std::false_type {
};

template <typename T, typename Alloc>
struct is_device_vector<thrust::device_vector<T, Alloc>>: std::true_type {
};

template <typename T, typename Alloc, Dimension... Dims>
struct is_device_vector<
    DescribedVector<thrust::device_vector<T, Alloc>, Dims...>>: std::true_type {
};

using MemoryResource = thrust::universal_host_pinned_memory_resource;

template <typename T>
using PinnedAllocator =
    thrust::mr::stateless_resource_allocator<T, MemoryResource>;

// Pipeline inputs
template <typename T>
using TAFTPVoltagesH =
    DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                    TimeDim,
                    AntennaDim,
                    FreqDim,
                    TimeDim,
                    PolnDim>;
template <typename T>
using TAFTPVoltagesD = DescribedVector<thrust::device_vector<T>,
                                       TimeDim,
                                       AntennaDim,
                                       FreqDim,
                                       TimeDim,
                                       PolnDim>;
// Beamformer inputs
template <typename T>
using FTPAVoltagesH =
    DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                    FreqDim,
                    TimeDim,
                    PolnDim,
                    AntennaDim>;
template <typename T>
using FTPAVoltagesD = DescribedVector<thrust::device_vector<T>,
                                      FreqDim,
                                      TimeDim,
                                      PolnDim,
                                      AntennaDim>;
// Coherent dedisperser inputs
template <typename T>
using TPAVoltagesH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                     TimeDim,
                                     PolnDim,
                                     AntennaDim>;
template <typename T>
using TPAVoltagesD =
    DescribedVector<thrust::device_vector<T>, TimeDim, PolnDim, AntennaDim>;
// Coherent beamformer outputs
template <typename T>
using TFBPowersH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                   TimeDim,
                                   FreqDim,
                                   BeamDim>;
template <typename T>
using TFBPowersD =
    DescribedVector<thrust::device_vector<T>, TimeDim, FreqDim, BeamDim>;
// Incoherent beamformer outputs
template <typename T>
using BTFPowersH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                   BeamDim,
                                   TimeDim,
                                   FreqDim>;
template <typename T>
using BTFPowersD =
    DescribedVector<thrust::device_vector<T>, BeamDim, TimeDim, FreqDim>;
// Incoherent dedisperser outputs, also used for skycleaver
template <typename T>
using TDBPowersH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                   TimeDim,
                                   DispersionDim,
                                   BeamDim>;
template <typename T>
using TDBPowersD =
    DescribedVector<thrust::device_vector<T>, TimeDim, DispersionDim, BeamDim>;
// Statistics outputs
template <typename T>
using FPAStatsH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                  FreqDim,
                                  PolnDim,
                                  AntennaDim>;
template <typename T>
using FPAStatsD =
    DescribedVector<thrust::device_vector<T>, FreqDim, PolnDim, AntennaDim>;

//skycleaver output vectors

template <typename T>
using TFPowersH = DescribedVector<thrust::host_vector<T, PinnedAllocator<T>>,
                                  TimeDim,
                                  FreqDim>;

} // namespace skyweaver

#endif // SKYWEAVER_DESCRIBEDVECTORS_HPP