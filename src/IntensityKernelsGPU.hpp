#pragma once
#include <DataTypes.hpp>

#ifdef USE_GPU
std::vector<double>
calculateIntensityGPU(std::vector<double> const& qVector,
                      std::vector<double> const& centers,
                      const std::vector<double> &counts);
// #endif
#else

inline std::vector<double>
calculateIntensityGPU(std::vector<double> const& qVector,
                      std::vector<double> const& centers,
                      const std::vector<double> &counts)
{
    std::cerr << "GPU not enabled\n";
    exit(1);
}

#endif