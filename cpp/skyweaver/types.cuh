#ifndef SKYWEAVER_TYPES_CUH
#define SKYWEAVER_TYPES_CUH

#include <iostream>

// Operator overloading for CUDA vector types
template <typename T>
struct is_vec4: std::false_type {
};

template <>
struct is_vec4<float4>: std::true_type {
};

template <>
struct is_vec4<char4>: std::true_type {
};

template <typename T>
inline constexpr bool is_vec4_v = is_vec4<T>::value;

template <typename T>
struct value_traits {};

template <>
struct value_traits<char>
{
    typedef char type;
    static constexpr char zero = 0;
    static constexpr char one = 1;
};

template <>
struct value_traits<float>
{
    typedef float type;
    static constexpr float zero = 0.0f;
    static constexpr float one = 1.0f;
};

template <>
struct value_traits<char4>
{
    typedef char type;
    static constexpr char4 zero = {0, 0, 0, 0};
    static constexpr char4 one = {1, 1, 1, 1};
};

template <>
struct value_traits<float4>
{
    typedef float type;
    static constexpr float4 zero = {0.0f, 0.0f, 0.0f, 0.0f};
    static constexpr float4 one = {1.0f, 1.0f, 1.0f, 1.0f};
};


inline std::ostream& operator<<(std::ostream& stream, char4 const& val) {
    stream << "(" << static_cast<int>(val.x) 
    << "," << static_cast<int>(val.y) 
    << "," << static_cast<int>(val.z) 
    << "," << static_cast<int>(val.w) 
    << ")";
    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, float4 const& val) {
    stream << "(" << val.x << "," << val.y << "," << val.z << "," << val.w << ")";
    return stream;
}

/**
 * vector - vector operations
 * explicit static_casts used to avoid Wnarrowing errors for char types due to integral promotion 
 * (over/underflow is the expected behaviour here).
 */
template <typename T>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T>, T>::type operator+(const T& lhs, const T& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x + rhs.x), 
            static_cast<typename value_traits<T>::type>(lhs.y + rhs.y),
            static_cast<typename value_traits<T>::type>(lhs.z + rhs.z),
            static_cast<typename value_traits<T>::type>(lhs.w + rhs.w)};
}

template <typename T>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T>, T>::type operator-(const T& lhs, const T& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x - rhs.x), 
            static_cast<typename value_traits<T>::type>(lhs.y - rhs.y),
            static_cast<typename value_traits<T>::type>(lhs.z - rhs.z),
            static_cast<typename value_traits<T>::type>(lhs.w - rhs.w)};
}

template <typename T>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T>, T>::type operator*(const T& lhs, const T& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x * rhs.x), 
            static_cast<typename value_traits<T>::type>(lhs.y * rhs.y),
            static_cast<typename value_traits<T>::type>(lhs.z * rhs.z),
            static_cast<typename value_traits<T>::type>(lhs.w * rhs.w)};
}

template <typename T>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T>, T>::type operator/(const T& lhs, const T& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x / rhs.x), 
            static_cast<typename value_traits<T>::type>(lhs.y / rhs.y),
            static_cast<typename value_traits<T>::type>(lhs.z / rhs.z),
            static_cast<typename value_traits<T>::type>(lhs.w / rhs.w)};
}

template <typename T>
__host__ __device__ inline  typename std::enable_if<is_vec4_v<T>, bool>::type operator==(const T& lhs, const T& rhs) {
    return  (lhs.x == rhs.x) && 
            (lhs.y == rhs.y) &&
            (lhs.z == rhs.z) &&
            (lhs.w == rhs.w);
}


/**
 * vector - scalar operations
 */
template <typename T, typename X>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T> && std::is_arithmetic_v<X>, T>::type operator*(const T& lhs, const X& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x * rhs), 
            static_cast<typename value_traits<T>::type>(lhs.y * rhs),
            static_cast<typename value_traits<T>::type>(lhs.z * rhs),
            static_cast<typename value_traits<T>::type>(lhs.w * rhs)};
}

template <typename T, typename X>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T> && std::is_arithmetic_v<X>, T>::type operator/(const T& lhs, const X& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x / rhs), 
            static_cast<typename value_traits<T>::type>(lhs.y / rhs),
            static_cast<typename value_traits<T>::type>(lhs.z / rhs),
            static_cast<typename value_traits<T>::type>(lhs.w / rhs)};
}

template <typename T, typename X>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T> && std::is_arithmetic_v<X>, T>::type operator+(const T& lhs, const X& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x + rhs), 
            static_cast<typename value_traits<T>::type>(lhs.y + rhs),
            static_cast<typename value_traits<T>::type>(lhs.z + rhs),
            static_cast<typename value_traits<T>::type>(lhs.w + rhs)};
}

template <typename T, typename X>
__host__ __device__ inline typename std::enable_if<is_vec4_v<T> && std::is_arithmetic_v<X>, T>::type operator-(const T& lhs, const X& rhs) {
    return {static_cast<typename value_traits<T>::type>(lhs.x - rhs), 
            static_cast<typename value_traits<T>::type>(lhs.y - rhs),
            static_cast<typename value_traits<T>::type>(lhs.z - rhs),
            static_cast<typename value_traits<T>::type>(lhs.w - rhs)};

}

#endif //SKYWEAVER_TYPES_CUH