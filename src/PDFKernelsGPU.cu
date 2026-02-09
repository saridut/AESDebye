#include <functional>
#include <iomanip>
#include <cuda/std/limits>

#include <DataTypesGPU.hpp>
#include <PDFKernelsGPU.hpp>

void setupStreamAttributes(const PositionsGPU &positions1_dev);
void printDeviceInfo(int64 start, int64 stop, ParallelHelper &parallelHelper, int &current_device);

std::tuple<dim3, dim3> getLaunchConfiguration(int64 start, int64 stop, int64 sizeJ, const CalculationConfig &config,
                                              int current_device);


#pragma region inline_functions


template <bool PseudoCoal>
__global__ void calculateHistogramKernel(int nPairs, int *cellPairsList,
                                         int64 *X1, int64 *Y1, int64 *Z1,
                                         int64 *X2, int64 *Y2, int64 *Z2,
                                         int *cellHeads1, int *cellHeads2, int *atomCounts1, int *atomCounts2,
                                         SBin *histogram, uint *counts, int64 *deltas,
                                         uint *countsOverflowCounts, int *deltaOverflowCounts,
                                         bool bothSame)
{
    const auto t = threadIdx.x;
    const auto iStart = blockIdx.y * blockDim.y + threadIdx.y;
    const auto iStride = gridDim.y * blockDim.y;
    for (auto pairIdx = iStart; pairIdx < nPairs; pairIdx += iStride)
    {
        const int pair1 = cellPairsList[2 * pairIdx];
        const int pair2 = cellPairsList[2 * pairIdx + 1];
        // setup the start and end indices for the cells
        const int i_start = cellHeads1[pair1];
        const int i_stop = i_start + atomCounts1[pair1]; // go till however many atoms are there in the cell
        const int j_start = cellHeads2[pair2];
        int j_stop = j_start + atomCounts2[pair2];

        const bool bothCellEqual = ((pair1 == pair2) && bothSame);

        for (int i = i_start; i < i_stop; i++)
        {
            const auto xi = X1[i];
            const auto yi = Y1[i];
            const auto zi = Z1[i];
            j_stop = bothCellEqual ? i : j_stop;

            for (int j = j_start; j < j_stop; j++)
            {
                if constexpr (PseudoCoal)
                {
                    __syncthreads();
                }

                const int64 dx = xi - X2[j];
                const int64 dy = yi - Y2[j];
                const int64 dz = zi - Z2[j];
                const int64 r2 = dx * dx + dy * dy + dz * dz;
                const auto Id = (uint)((sqrt((float)r2) + 500) / 1000);
                const int64 id_ = (int64)Id;
                const int64 delta = r2 - (id_ * id_ * 1000000LL);
                if constexpr (PseudoCoal)
                {
                    SBin &currentHistBin = histogram[Id];
                    const auto ret = atomicAdd(
                        (0 == t) ? &currentHistBin.count : reinterpret_cast<uint64 *>(&currentHistBin.delta),
                        (0 == t) ? 1 : delta);

                    if (1 == t)
                    {
                        const auto oldDelta = (int64)ret;
                        if (delta > 0 && oldDelta > cuda::std::numeric_limits<int64>::max() - delta)
                        {
                            atomicAdd(&deltaOverflowCounts[Id], 1);
                        }
                        else if (delta < 0 && oldDelta < cuda::std::numeric_limits<int64>::min() - delta)
                        {
                            atomicAdd(&deltaOverflowCounts[Id], -1);
                        }
                    }
                    // no need to check for counts overflow, the counter is already uint64
                    // for SBin.
                }
                else
                {
                    const auto cnt = atomicAdd(&counts[Id], 1);
                    const auto oldDelta = (int64)atomicAdd(reinterpret_cast<uint64 *>(&deltas[Id]), delta);

                    if (cnt == UINT_MAX)
                        atomicAdd(&countsOverflowCounts[Id], 1);
                    if (delta > 0 && oldDelta > cuda::std::numeric_limits<int64>::max() - delta)
                    {
                        atomicAdd(&deltaOverflowCounts[Id], 1);
                    }
                    else if (delta < 0 && oldDelta < cuda::std::numeric_limits<int64>::min() - delta)
                    {
                        atomicAdd(&deltaOverflowCounts[Id], -1);
                    }
                }
            }
        }
    }
}

template <bool PseudoCoal>
__global__ void calculateHistogramKernel(int64 start, int64 N1, int64 N2, int64 *X1, int64 *Y1, int64 *Z1,
                                         int64 *X2, int64 *Y2, int64 *Z2,
                                         SBin *histogram, uint *counts, int64 *deltas,
                                         uint *countsOverflowCounts, int *deltaOverflowCounts,
                                         bool bothSame)
{
    const auto t = threadIdx.x; // this will be ignored if PseudoCoal is false
    const auto iStart = blockIdx.y * blockDim.y + threadIdx.y;
    const auto iStride = gridDim.y * blockDim.y;

    const auto jStart = blockIdx.z * blockDim.z + threadIdx.z;
    const auto jStride = gridDim.z * blockDim.z;

    for (auto i = iStart + start; i < N1; i += iStride)
    {
        const auto xi = X1[i];
        const auto yi = Y1[i];
        const auto zi = Z1[i];

        if (bothSame)
            N2 = i;

        for (auto j = jStart; j < N2; j += jStride)
        {
            if constexpr (PseudoCoal)
            {
                __syncthreads();
            }

            const int64 dx = xi - X2[j];
            const int64 dy = yi - Y2[j];
            const int64 dz = zi - Z2[j];
            const int64 r2 = dx * dx + dy * dy + dz * dz;
            const auto Id = (uint)((sqrt((float)r2) + 500) / 1000);
            const int64 id_ = (int64)Id;
            const int64 delta = r2 - (id_ * id_ * 1000000LL);

            if constexpr (PseudoCoal)
            {
                SBin &currentHistBin = histogram[Id];
                const auto ret = atomicAdd(
                    (0 == t) ? &currentHistBin.count : reinterpret_cast<uint64 *>(&currentHistBin.delta),
                    (0 == t) ? 1 : delta);

                if (1 == t)
                {
                    const auto oldDelta = (int64)ret;
                    if (delta > 0 && oldDelta > cuda::std::numeric_limits<int64>::max() - delta)
                    {
                        atomicAdd(&deltaOverflowCounts[Id], 1);
                    }
                    else if (delta < 0 && oldDelta < cuda::std::numeric_limits<int64>::min() - delta)
                    {
                        atomicAdd(&deltaOverflowCounts[Id], -1);
                    }
                }

                // no need to check for counts overflow, the counter is already uint64
                // for SBin.
            }
            else
            {
                const auto cnt = atomicAdd(&counts[Id], 1);
                const auto oldDelta = (int64)atomicAdd(reinterpret_cast<uint64 *>(&deltas[Id]), delta);

                if (cnt == UINT_MAX)
                    atomicAdd(&countsOverflowCounts[Id], 1);
                if (delta > 0 && oldDelta > cuda::std::numeric_limits<int64>::max() - delta)
                {
                    atomicAdd(&deltaOverflowCounts[Id], 1);
                }
                else if (delta < 0 && oldDelta < cuda::std::numeric_limits<int64>::min() - delta)
                {
                    atomicAdd(&deltaOverflowCounts[Id], -1);
                }
            }
        }
    }
}



void printLaunchConfig(dim3 blockSize, dim3 gridSize, ParallelHelper& ParallelHelper)
{
    ParallelHelper << "Block size: " << blockSize.x << ", " << blockSize.y << ", " << blockSize.z << "\n";
    ParallelHelper << "Grid size: " << gridSize.x << ", " << gridSize.y << ", " << gridSize.z << "\n";
}

dim3 maxGrid(dim3 blockSize, int64 sizeI, int64 sizeJ)
{
    return dim3(1, std::min((int)std::ceil((sizeI + blockSize.y - 1) / blockSize.y), 65535),
                (std::min((int)std::ceil((sizeJ + blockSize.z - 1) / blockSize.z), 65535)));
}

dim3 maxGridCells(dim3 blockSize, int64 sortedCellListSize)
{
    return dim3(1, std::min((int)std::ceil((sortedCellListSize + blockSize.y - 1) / blockSize.y), 65535));
}

void calculatePDFGPU(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                     const PairsList &sortedCellPairs, int64 start, int64 stop, bool samePositions,
                     ParallelHelper &parallelHelper, CalculationConfig &config)
{
    Positions gpuPositionsI(positionsI);

    // get number of devices
    int deviceCount; cudaGetDeviceCount(&deviceCount);
    int current_device{parallelHelper.worldRank % deviceCount};
    cudaSetDevice(current_device);
    printDeviceInfo(start, stop, parallelHelper, current_device);

    auto positions1_dev = PositionsGPU(gpuPositionsI);
    auto positions2_dev = PositionsGPU(positionsJ);
    auto histogram_dev = HistogramGPU();
    auto sortedPairsList_dev = CellPairsGPU(sortedCellPairs);
    auto [blockSize, gridSize] = getLaunchConfiguration(start, stop, positionsJ.size(), config, current_device);

    setupStreamAttributes(positions1_dev);

    if (config.pseudoCoal)
    {
        blockSize.x = 2; // for updating the delta
        if (!config.useGPUCellList)
        {
            calculateHistogramKernel<true><<<gridSize, blockSize>>>(start, stop, positionsJ.size(),
                                                                    positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
                                                                    positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
                                                                    histogram_dev.histogram_SBin,
                                                                    histogram_dev.counts, histogram_dev.deltas,
                                                                    histogram_dev.countsOverFlowCounts,
                                                                    histogram_dev.deltaOverflowCounts,
                                                                    samePositions);
        }
        else
        {
            parallelHelper << "Using GPU cell list" << "\n";
            blockSize = dim3(2, 512);
            gridSize = config.fillGPU ? maxGridCells(blockSize, sortedCellPairs.size()) : dim3(1, (int) 108 * 32);
            printLaunchConfig(blockSize, gridSize, parallelHelper);
            calculateHistogramKernel<true><<<gridSize, blockSize>>>(sortedCellPairs.size(),
                                                                    sortedPairsList_dev.sortedPairsList,
                                                                    positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
                                                                    positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
                                                                    positions1_dev.cellHeads, positions2_dev.cellHeads,
                                                                    positions1_dev.atomCounts, positions2_dev.atomCounts,
                                                                    histogram_dev.histogram_SBin,
                                                                    histogram_dev.counts, histogram_dev.deltas,
                                                                    histogram_dev.countsOverFlowCounts,
                                                                    histogram_dev.deltaOverflowCounts,
                                                                    samePositions);
        }
    }
    else
    {
        if (!config.useGPUCellList)
        {
            calculateHistogramKernel<false><<<gridSize, blockSize>>>(start, stop, positionsJ.size(),
                                                                     positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
                                                                     positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
                                                                     histogram_dev.histogram_SBin,
                                                                     histogram_dev.counts, histogram_dev.deltas,
                                                                     histogram_dev.countsOverFlowCounts,
                                                                     histogram_dev.deltaOverflowCounts,
                                                                     samePositions);
        }
        else
        {
            parallelHelper << "Using GPU cell list" << "\n";
            blockSize = dim3(1, 1024);
            gridSize = config.fillGPU ? maxGridCells(blockSize, sortedCellPairs.size()) : dim3(1, (int) 108 * 32);
            printLaunchConfig(blockSize, gridSize, parallelHelper);
            calculateHistogramKernel<false><<<gridSize, blockSize>>>(sortedCellPairs.size(),
                                                                    sortedPairsList_dev.sortedPairsList,
                                                                    positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
                                                                    positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
                                                                    positions1_dev.cellHeads, positions2_dev.cellHeads,
                                                                    positions1_dev.atomCounts, positions2_dev.atomCounts,
                                                                    histogram_dev.histogram_SBin,
                                                                    histogram_dev.counts, histogram_dev.deltas,
                                                                    histogram_dev.countsOverFlowCounts,
                                                                    histogram_dev.deltaOverflowCounts,
                                                                    samePositions);
        }
    }

    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        // manually call destroy to avoid memory leak on the gpu
        positions1_dev.destroy();
        positions2_dev.destroy();
        sortedPairsList_dev.destroy();
        histogram_dev.destroy();
        throw std::runtime_error(cudaGetErrorString(error));
    }

    // Hybrid parallelism for calculating the PDF
//    if (config.hybrid)
//    {
//        calculatePDFCPU(positionsI, positionsJ, pdfVector,
//                        sortedCellPairs, start, stop, samePositions,
//                        parallelHelper, config);
//    }

    cudaDeviceSynchronize();
    histogram_dev.copyToHost(pdfVector, config.pseudoCoal);
}

std::tuple<dim3, dim3> getLaunchConfiguration(int64 start, int64 stop, int64 sizeJ,  const CalculationConfig &config,
                                              int current_device) {
//    int threads= 1024;
    dim3 blockSize= ((!config.useLocalHistogram) || (config.pseudoCoal)) ? dim3(1, 16, 16) : dim3(1, 32, 32);
    dim3 gridSize = maxGrid(blockSize, stop - start, sizeJ);
    if (!config.fillGPU)
    {
        blockSize = dim3(1, 512, 1);
        gridSize = dim3(1, 108, 32);
        cudaDeviceProp deviceProp{};
        cudaGetDeviceProperties(&deviceProp, current_device); // 0-th device
        gridSize.y = deviceProp.multiProcessorCount;
    }

    return {blockSize, gridSize};
}

void printDeviceInfo(int64 start, int64 stop, ParallelHelper &parallelHelper, int &current_device) {
    cudaDeviceProp device_prop{};
    cudaGetDevice(&current_device);
    cudaGetDeviceProperties(&device_prop, current_device);
    parallelHelper << "GPU: " << device_prop.name << "\n";
    parallelHelper << "L2 Cache Size: " << device_prop.l2CacheSize / 1024 / 1024
              << " MB" << "\n";
    parallelHelper << "Max Persistent L2 Cache Size: "
              << device_prop.persistingL2CacheMaxSize / 1024 / 1024 << " MB"
              << "\n";

    std::cout  << "Device: " << current_device << ", start:" << start << ", stop: " << stop << "\n";
}

void setupStreamAttributes(const PositionsGPU &positions1_dev) {
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    cudaStreamAttrValue stream_attribute;                                                      // Stream level attributes data structure
    stream_attribute.accessPolicyWindow.base_ptr = reinterpret_cast<void *>(positions1_dev.X); // Global Memory data pointer
    stream_attribute.accessPolicyWindow.num_bytes = positions1_dev.size * 3 * sizeof(int64);     // Number of bytes for persistence access.

// (Must be less than cudaDeviceProp::accessPolicyMaxWindowSize)
    stream_attribute.accessPolicyWindow.hitRatio = 1.;                                         // Hint for cache hit ratio
    stream_attribute.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;                // Type of access property on cache hit
    stream_attribute.accessPolicyWindow.missProp = cudaAccessPropertyStreaming;                // Type of access property on cache miss.

    // Set the attributes to a CUDA stream of type cudaStream_t
    cudaStreamSetAttribute(stream, cudaStreamAttributeAccessPolicyWindow, &stream_attribute);
//    checkLastError(false);

    // check if there are any errors in stream setup
    cudaStreamQuery(stream);

    // check for errors
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        std::cerr << "CUDA error: " << cudaGetErrorString(error) << std::endl;
    }

    cudaDeviceSynchronize();
}

// __global__ void calculateHistogramKernel(int nPairs, int *cellPairsList,
//                                         int64 *X1, int64 *Y1, int64 *Z1,
//                                          int64 *X2, int64 *Y2, int64 *Z2,
//                                         int *cellHeads1, int *cellHeads2, int *atomCounts1, int *atomCounts2,
//                                          uint *counts, uint *countsOverflowCounts,
//                                          int *deltaOverflowCounts, int64 *deltas,
//                                          bool bothSame)
// {
//     const auto iStart = blockIdx.x * blockDim.x + threadIdx.x;
//     const auto iStride = gridDim.x * blockDim.x;
//     for (auto pairIdx = iStart; pairIdx < nPairs; pairIdx += iStride)
//     {
//         const int pair1 = cellPairsList[2*pairIdx];
//         const int pair2 = cellPairsList[2*pairIdx + 1];
//         // setup the start and end indices for the cells
//         const int i_start = cellHeads1[pair1];
//         const int i_stop = i_start + atomCounts1[pair1]; // go till however many atoms are there in the cell
//         const int j_start = cellHeads2[pair2];
//         int j_stop = j_start + atomCounts2[pair2];

//         const bool bothCellEqual = ((pair1 == pair2) && bothSame);

//         for (int i = i_start; i < i_stop; i++)
//         {
//             const auto xi = X1[i];
//             const auto yi = Y1[i];
//             const auto zi = Z1[i];

//             j_stop = bothCellEqual ? i : j_stop;
//             for (int j = j_start; j < j_stop; j++)
//             {
//                 const int64 dx = xi - X2[j];
//                 const int64 dy = yi - Y2[j];
//                 const int64 dz = zi - Z2[j];
//                 const int64 r2 = dx * dx + dy * dy + dz * dz;
//                 const auto Id = (uint)((sqrt((float)r2) + 500) / 1000);
//                 const int64 id_ = (int64)Id;
//                 const auto cnt = atomicAdd(&counts[Id], 1);
//                 if (cnt == UINT_MAX)
//                     atomicAdd(&countsOverflowCounts[Id], 1);

//                 const int64 delta = r2 - (id_ * id_ * 1000000LL);

//                 const auto oldDelta = (int64)atomicAdd(reinterpret_cast<uint64 *>(&deltas[Id]), delta);
//                 if (delta > 0 && oldDelta > cuda::std::numeric_limits<int64>::max() - delta)
//                 {
//                     atomicAdd(&deltaOverflowCounts[Id], 1);
//                 }
//                 if (delta < 0 && oldDelta < cuda::std::numeric_limits<int64>::min() - delta)
//                 {
//                     atomicAdd(&deltaOverflowCounts[Id], -1);
//                 }
//             }
//         }
//     }
// }

//   // celllist based histogram calculation
//     blockSize = dim3(64);
//     gridSize = dim3(std::min((int) std::ceil((cellPairsList.size() + blockSize.x - 1) / blockSize.x), 65535));

//     startGPU = get_wall_time();
//     calculateHistogramKernel<<<gridSize, blockSize>>>(cellPairsList.size(), cellPairsList_dev.sortedCellList,
//                                                         positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
//                                                         positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
//                                                         positions1_dev.cellHeads, positions2_dev.cellHeads,
//                                                         positions1_dev.atomCounts, positions2_dev.atomCounts,
//                                                         histogram_dev.counts, histogram_dev.countsOverFlowCounts,
//                                                         histogram_dev.deltaOverflowCounts, histogram_dev.deltas,
//                                                         (&positions1 == &positions2));

//     numAtom = positions1.size();
//     checkLastError();
//     cudaDeviceSynchronize();
//     endGPU = get_wall_time();
//     timeGPU = endGPU - startGPU;

//     std::cout << "CellList" << ":" << std::endl;
//     std::cout << "   Time GPU:  " << timeGPU << std::endl;
//     std::cout << "   GPD/s GPU: " << 1e-9 * numAtom * (numAtom - 1) / 2 / timeGPU << std::endl;

//     histogram_dev.copyToHost(histogram, false);

// if (pseudoCoal)
// {
//     dim3 blockSize(2, 512, 1);
//     dim3 gridSize(1, 108, 32);
//     // dim3 gridSize(1, std::min((int)std::ceil((positions1.size() + blockSize.y - 1) / blockSize.y), 65535),
//     //             (std::min((int)std::ceil((positions1.size() + blockSize.z - 1) / blockSize.z), 65535)));
//     printLaunchConfig(blockSize, gridSize);
//     calculateHistogramKernel<<<gridSize, blockSize>>>(positions1.size(), positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
//                                                       histogram_dev.histogram_SBin, histogram_dev.deltaOverflowCounts);
// }
// else
// {
//     dim3 blockSize(32, 32);
//     if (crystalline)
//     {
//         blockSize = dim3(16, 16); // (16, 16) works best in case of ideal crystal
//     }

//     dim3 gridSize(std::min((int)std::ceil((positions1.size() + blockSize.x - 1) / blockSize.x), 65535),
//                   (std::min((int)std::ceil((positions1.size() + blockSize.y - 1) / blockSize.y), 65535)));

//     printLaunchConfig(blockSize, gridSize);
//     calculateHistogramKernel<<<gridSize, blockSize>>>(positions1.size(), positions2.size(),
//                                                       positions1_dev.X, positions1_dev.Y, positions1_dev.Z,
//                                                       positions2_dev.X, positions2_dev.Y, positions2_dev.Z,
//                                                       histogram_dev.counts, histogram_dev.countsOverFlowCounts,
//                                                       histogram_dev.deltaOverflowCounts, histogram_dev.deltas);
// }
