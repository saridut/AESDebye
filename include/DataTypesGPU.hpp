#pragma once

#include <DataTypes.hpp>
#include <PDF.hpp>
#include <Positions.hpp>
#include <CellList.hpp>


struct SBin {
    uint64 count{0};
    int64 delta{0};

    SBin() = default;
};


struct PositionsGPU {
    int64 *X, *Y, *Z;
    int *cellHeads;
    int *atomCounts;
    size_t size, nAtoms;

    PositionsGPU(const Positions &positions) {
        nAtoms = positions.X.size();
        size = nAtoms * sizeof(int64);

        cudaMalloc((void **) &X, size);
        cudaMalloc((void **) &Y, size);
        cudaMalloc((void **) &Z, size);

        cudaMemcpy(X, positions.X.data(), size, cudaMemcpyHostToDevice);
        cudaMemcpy(Y, positions.Y.data(), size, cudaMemcpyHostToDevice);
        cudaMemcpy(Z, positions.Z.data(), size, cudaMemcpyHostToDevice);

        cudaMalloc((void **) &cellHeads, positions.cellHeads.size() * sizeof(int));
        cudaMalloc((void **) &atomCounts, positions.atomCounts.size() * sizeof(int));
        cudaMemcpy(cellHeads, positions.cellHeads.data(), positions.cellHeads.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(atomCounts, positions.atomCounts.data(), positions.atomCounts.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
    }

    void destroy() {
        cudaFree(X);
        cudaFree(Y);
        cudaFree(Z);
        cudaFree(cellHeads);
        cudaFree(atomCounts);
    }

    ~PositionsGPU() {
        destroy();
    }
};

struct HistogramGPU {
    size_t numBins;
    size_t size;
    SBin *histogram_SBin;
    size_t size_SBin;

    // AOS
    uint *counts;
    uint *countsOverFlowCounts;
    int *deltaOverflowCounts;
    int64 *deltas;

    HistogramGPU() {
        // numBins = histogram.size();
        numBins = N_BINS;
        size_SBin = numBins * sizeof(SBin);
        cudaMalloc((void **) &this->histogram_SBin, size_SBin);
        cudaMemset(this->histogram_SBin, 0, size_SBin);

        // SOA
        cudaMalloc((void **) &counts, numBins * sizeof(uint));
        cudaMalloc((void **) &countsOverFlowCounts, numBins * sizeof(uint));
        cudaMalloc((void **) &deltaOverflowCounts, numBins * sizeof(int));
        cudaMalloc((void **) &deltas, numBins * sizeof(uint64));

        cudaMemset(counts, 0, numBins * sizeof(uint));
        cudaMemset(countsOverFlowCounts, 0, numBins * sizeof(uint));
        cudaMemset(deltaOverflowCounts, 0, numBins * sizeof(int));
        cudaMemset(deltas, 0, numBins * sizeof(int64));

    }

    void reset() {
        cudaMemset(this->histogram_SBin, 0, size_SBin);
        cudaMemset(counts, 0, numBins * sizeof(uint));
        cudaMemset(countsOverFlowCounts, 0, numBins * sizeof(uint));
        cudaMemset(deltaOverflowCounts, 0, numBins * sizeof(int));
        cudaMemset(deltas, 0, numBins * sizeof(int64));
    }


    void copyToHost(std::vector<HistBin> &histogram, bool useSBin = false) { // only works with SOA
        // empty counts, countsOverflowCounts, deltaOverflowCounts, deltas on cpu
        std::vector<uint> counts(numBins, 0);
        std::vector<uint> countsOverFlowCounts(numBins, 0);
        std::vector<int> deltaOverflowCounts(numBins, 0);
        std::vector<int64> deltas(numBins, 0);
        std::vector<SBin> histogram_SBin(numBins);

        // copy counts, countsOverflowCounts, deltaOverflowCounts, deltas from gpu
        cudaMemcpy(counts.data(), this->counts, numBins * sizeof(uint), cudaMemcpyDeviceToHost);
        cudaMemcpy(countsOverFlowCounts.data(), this->countsOverFlowCounts, numBins * sizeof(uint),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(deltaOverflowCounts.data(), this->deltaOverflowCounts, numBins * sizeof(int),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(deltas.data(), this->deltas, numBins * sizeof(int64), cudaMemcpyDeviceToHost);
        cudaMemcpy(histogram_SBin.data(), this->histogram_SBin, numBins * sizeof(SBin), cudaMemcpyDeviceToHost);

        // copy to histogram datastructure 

        uint64 total_count = 0;
        if (useSBin) {
#pragma omp parallel for schedule(static) reduction(+:total_count)
            for (size_t i = 0; i < numBins; i++) {
                histogram[i].count = histogram_SBin[i].count;
                total_count += histogram_SBin[i].count;
                histogram[i].delta = histogram_SBin[i].delta;
                histogram[i].deltaOverflowCount = 2 * deltaOverflowCounts[i];
            }
        } else {
#pragma omp parallel for schedule(static) reduction(+:total_count)
            for (size_t i = 0; i < numBins; i++) {
                histogram[i].count = (uint64) counts[i] + (uint64) countsOverFlowCounts[i] * UINT32_MAX;
                total_count += histogram[i].count;
                // these deltas won't match the CPU deltas due to not offsetting by INT64_MAX
                histogram[i].delta = deltas[i];

                // 2 here is important, because we don't offset the delta by INT64_MAX in the GPU
                // But we do it in the CPU
                // instead its gets offset by 2 * INT64_MAX because of the cyclic overflow/underflow
                // important to check if not offsetting is a better approach than offsetting, because,
                // it seems more likely that the overflow counter will get accesed more often if we don't offset
                // for example delta close ~ max, we add delta > 0, it goes close to min, now if we add delta < 0
                // it will again underflow and go to max - leading to update to the overflow/underflow counter
                histogram[i].deltaOverflowCount = 2 * deltaOverflowCounts[i];
            }
        }
    }

    void destroy() {
        cudaFree(histogram_SBin);
        cudaFree(counts);
        cudaFree(countsOverFlowCounts);
        cudaFree(deltaOverflowCounts);
        cudaFree(deltas);
    }

    ~HistogramGPU() {
        destroy();
    }
};


struct CellPairsGPU {
    int *sortedPairsList{};
    size_t size;


    CellPairsGPU(const std::vector<std::array<int, 2>> &sortedPairsList) {
        size = sortedPairsList.size() * 2 * sizeof(int);

        // flatten the sortedPairsList to a 1D array

        std::vector<int> flattenedPairsList;
        for (auto pair: sortedPairsList) {
            flattenedPairsList.push_back(pair[0]);
            flattenedPairsList.push_back(pair[1]);
        }

        cudaMalloc((void **) &this->sortedPairsList, flattenedPairsList.size() * sizeof(int));
        cudaMemcpy(this->sortedPairsList, flattenedPairsList.data(), flattenedPairsList.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
    }

    void destroy() {
        cudaFree(sortedPairsList);
    }

    ~CellPairsGPU() {
        destroy();
    }
};
