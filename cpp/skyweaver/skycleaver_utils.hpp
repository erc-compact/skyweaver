#ifndef SKYWEAVER_SKYCLEAVER_UTILS_HPP
#define SKYWEAVER_SKYCLEAVER_UTILS_HPP
namespace skyweaver
{
    
template <typename T>
std::vector<T>
get_list_from_string(const std::string& value,
                     T epsilon = std::numeric_limits<T>::epsilon())
{
    std::vector<T> output;
    std::vector<std::string> comma_chunks;

    // Split the input string by commas
    std::stringstream ss(value);
    std::string token;
    while(std::getline(ss, token, ',')) { comma_chunks.push_back(token); }

    for(const auto& comma_chunk: comma_chunks) {
        // Check if the chunk contains a colon (indicating a range)
        if(comma_chunk.find(':') == std::string::npos) {
            output.push_back(static_cast<T>(std::atof(comma_chunk.c_str())));
            continue;
        }

        // Split the range chunk by colons
        std::stringstream ss_chunk(comma_chunk);
        std::vector<T> colon_chunks;
        std::string colon_token;
        while(std::getline(ss_chunk, colon_token, ':')) {
            colon_chunks.push_back(
                static_cast<T>(std::atof(colon_token.c_str())));
        }

        // Determine the step size
        T step = colon_chunks.size() == 3 ? colon_chunks[2] : static_cast<T>(1);
        T start = colon_chunks[0];
        T stop  = colon_chunks[1];

        // Loop and add values to the output vector
        if constexpr(std::is_floating_point<T>::value) {
            for(T k = start; k <= stop + epsilon; k += step) {
                output.push_back(k);
            }
        } else {
            for(T k = start; k <= stop; k += step) { output.push_back(k); }
        }
    }
    return output;
}

}
#endif // SKYWEAVER_SKYCLEAVER_UTILS_HPP