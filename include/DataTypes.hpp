#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cfloat>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

// type definitions
typedef long long int64;
typedef unsigned long long uint64;

// Constants
constexpr double DEBYE_MAX_Rd = 1'000'000'000.0; // (1e9) max value to scale up the positions to
constexpr size_t N_BINS = 1'732'100;
constexpr size_t N_SMALL_BINS = 173'210;
constexpr size_t N_STRIP_MINED = 512;
constexpr int64 MAX_BIN_ERROR = 1'732'109'250'000; // see below for calculation - overestimation for normal case by a bit
constexpr int MAX_DELTA_UPDATE_COUNT = 5'324'935; // 964 ordinaryly, but underestimated a bit

/**
 * @brief Calculation Configuration
 * 
 * This struct holds the configuration for the calculation. It is used to control the flow of the calculation.
 * 
 */
struct CalculationConfig
{
    bool useCellList{true}; //< Use cell list for the calculation (CPU and position sorting for GPU)
    bool useMPI{false};     //< Use MPI for the calculation (CPU+Multi-GPU)
    bool useGPU{false};     //< Use GPU for the calculation (GPU)
    // bool hybrid{false};     //< Use hybrid calculation (CPU+GPU)
    bool useLocalHistogram{true}; //< Use local histogram for the calculation (CPU)
    bool smallBins{false}; //< Use small bins for the calculation (CPU+GPU)
    bool pseudoCoal{true}; //< Use pseudo coalescing for the calculation (GPU)
    bool fillGPU{false}; //< Fill the GPU with threads (GPU)
    bool useGPUCellList{false}; //< Whether to use GPU cell list or not (GPU)

    // overload the ofstream operator to print the config
    friend std::ostream &operator<<(std::ostream &os, const CalculationConfig &config)
    {
        os << ">> Config: useCellList: " << config.useCellList << ", useMPI: " << config.useMPI << ", useGPU: " << config.useGPU
           << ", useLocalHistogram: " << config.useLocalHistogram << ", smallBins: " << config.smallBins << "\n";
        os << ">> GPU specific: pseudoCoal: " << config.pseudoCoal << " fillGPU: " << config.fillGPU << " useGPUCellList: " << config.useGPUCellList << "\n";
        return os;
    }

};



namespace helpers{

// Global helpers, need to be static since they are used in multiple files

/**
 * @brief Get the current wall time
 *
 * @return double
 */
inline double get_wall_time() {
    using namespace std::chrono;
    return duration_cast<duration<double>>(system_clock::now().time_since_epoch()).count();
}

inline std::string bool_to_string(bool b) {
    return b ? "\033[32mtrue\033[0m" : "\033[31mfalse\033[0m";
}

/**
* @brief Split a string by a delimiter
*
* @param inputString The string to split
* @param delimiter The delimiter to split by
* @return std::vector<std::string> splitted segments
*/
inline std::vector<std::string> stringSplit(std::string inputString, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = inputString.find(delimiter, pos_start)) != std::string::npos) {
        token = inputString.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(inputString.substr(pos_start));
    return res;
}


/**
* @brief Creates a linspace vector

*/
inline std::vector<double> linspace(double start, double end, int nSteps) {
    std::vector<double> result;
    result.reserve(nSteps);  // Reserve space for the vector
    if (nSteps < 2) {
        result.push_back(start);  // Return single value if nSteps < 2
        return result;
    }

    double step = (end - start) / (nSteps - 1);  // Calculate step size

    for (int i = 0; i < nSteps; ++i) {
        result.push_back(start + i * step);  // Add the values to the vector
    }

    return result;
}
}



/* python code to get the bin error and delta count
scaling
scaling_factor = 1_000_000_0 # remove last zero
rh = 173210 # add one zero
limit = rh * rh + rh + 1
limit *= scaling_factor
limit -= 750000
center = rh * rh * scaling_factor
int64_max = 9223372036854775807
print(f"max bin error {limit - center}")

// Explanation from Rose-X
@@@ (A+B)^2 = A^2 + B^2 + 2AB taken B = 0.50 it will be A^2 + 0.25 + A finally by ceiling the value it is A^2 + A + 1 - 0.75
Bin limits and values for histogram -
limits - lower limit of center^2 - but calculated using center of the bin
limits = ciel((center+0.50)^2) = center^2 - center + 1 - still in 1e12 = * 1000000 to get 1e18
shift the limit by 0.75 to get the lower limit of the bin ??
value = center^2 - direct scaling to 1e18
*/