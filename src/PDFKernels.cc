/**
 * @file PDFCPUKernels.cc
 * This file contains the main CPU kernels for the PDF calculation
 * 
 * There are two implementations, one with and one without cell list
 * the one without cell list only meant for very tiny systems
 * 
 */
#include <PDFKernels.hpp>


// pdfVector reduction declaration
#pragma omp declare                                                                                                                         \
    reduction(                                                                                                                              \
            + : std::vector<HistBin> : std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<HistBin>())) \
    initializer(                                                                                                                            \
            omp_priv = decltype(omp_orig)(omp_orig.size()))


template
void calculatePDFKernel<true>(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                              int64 start, int64 stop, bool samePositions);

template
void
calculatePDFKernel<false>(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                          int64 start, int64 stop, bool samePositions);


template
void calculatePDFKernel<true>(Positions const &positionsI, Positions const &positionsJ, PairsList const &cellPairsList,
                              std::vector<HistBin> &pdfVector, int64 start, int64 stop, bool samePositions);

template
void calculatePDFKernel<false>(Positions const &positionsI, Positions const &positionsJ, PairsList const &cellPairsList,
                               std::vector<HistBin> &pdfVector, int64 start, int64 stop, bool samePositions);

/**
 * @brief LocalBin struct to store the values before actually offsetting them to the thread PDF
 * 
 * This requires less memory than the HistBin struct.
 */
struct LocalBin // 128
{
    uint count{0}; /*< the cout of the bin*/
    uint tillOverFlowReset{MAX_DELTA_UPDATE_COUNT}; /*< the number of updates until the overflow is reset*/
    int64 delta{
            0}; /*< the delta value of the bin - this will never overflow, since the counter above will run out first*/
};

#pragma region inline_functions

/**
 * @brief Find the bin index and the delta value for the given pair of positions ii and jj
 * 
 * @param positions1 Positions container 1
 * @param positions2 Positions container 2
 * @param ii index of the first position
 * @param jj index of the second position
 * @param binId the bin index calculated - notice that this is a reference
 * @param delta delta value calculated - notice that this is a reference
 */
inline void findBinIndex(Positions const &positions1, Positions const &positions2,
                         uint ii, uint jj,
                         int64 &binId, int64 &delta, float scaleDown, int64 scaleUp) {
    int64 dx = positions1.X[ii] - positions2.X[jj];
    int64 dy = positions1.Y[ii] - positions2.Y[jj];
    int64 dz = positions1.Z[ii] - positions2.Z[jj];
    int64 r2 = dx * dx + dy * dy + dz * dz;
    binId = (int64) ((std::sqrt((float) r2) + 500) * scaleDown);
    delta = r2 - (binId * binId * scaleUp);
}


/**
 * @brief Update the delta value of the global bin and check for overflow
 * 
 * This function is only applicable in cases where LocalBin is not used. 
 * 
 * @param globalBin the global bin to update
 * @param delta the delta value to add
 */
inline void updateDeltaOverflow(HistBin &globalBin, int64 delta) {
    // first check if delta is positive or negative, since the calculations are required to be
    // done in uint64, 
    // finally update the tillOverFlowReset value
    if (globalBin.delta > 0) {
        uint64 diff = LLONG_MAX - globalBin.delta;
        diff /= MAX_BIN_ERROR; // max bin error - if you subtract the center of the last bin from the edge.
        diff = (diff > (int64) INT_MAX) ? INT_MAX : diff;
        globalBin.tillOverFlowReset = (int) diff - 1; // -1 isnt tested yet
    } else if (globalBin.delta < 0) {
        uint64 diff = (LLONG_MAX) + globalBin.delta;
        diff /= MAX_BIN_ERROR; // max bin error
        diff = (diff > (int64) INT_MAX) ? INT_MAX : diff;
        globalBin.tillOverFlowReset = (int) diff - 1; // -1 isnt tested yet
    }

    // now update the overflow values, only if both the current and incoming values have same sign
    // otherwise we know that there will be no overflow/underflow
    if ((globalBin.delta > 0) && (delta > 0)) {
        // check if adding the delta will cause an overflow
        if (LLONG_MAX - globalBin.delta < delta) {
            // if so, then offset the delta
            globalBin.delta -= LLONG_MAX;

            // update the overflow count
            globalBin.deltaOverflowCount++;

            // reset the overflow counter
            globalBin.tillOverFlowReset = MAX_DELTA_UPDATE_COUNT;
        }
    } else if ((globalBin.delta < 0) && (delta < 0)) // same as above, but for negative values
    {
        if (-LLONG_MAX - globalBin.delta > delta) {
            // offset the delta
            globalBin.delta += LLONG_MAX;

            // reduce the overflow count (since underflow occurs)
            globalBin.deltaOverflowCount--;

            // reset the overflow counter
            globalBin.tillOverFlowReset = MAX_DELTA_UPDATE_COUNT;
        }
    }
}

/**
 * @brief Function to offload the LocalBin values to HistBin
 * 
 * @param currentHistBin Local bin to offload
 * @param globalBin Global bin to offload to
 */
inline void addToGlobalBin(LocalBin &currentHistBin, HistBin &globalBin) {

    // only need to care if both have the same sign
    if ((globalBin.delta > 0) && (currentHistBin.delta > 0)) {
        // check if adding the delta will cause an overflow
        if (LLONG_MAX - globalBin.delta < currentHistBin.delta) {
            // offset the delta
            globalBin.delta -= LLONG_MAX;

            // update the overflow count
            globalBin.deltaOverflowCount++;
        }
    } else if ((globalBin.delta < 0) && (currentHistBin.delta < 0)) // again, same as above, but for negative values
    {
        // check if subtracting the delta will cause an underflow
        if (-LLONG_MAX - globalBin.delta > currentHistBin.delta) {
            // offset the delta
            globalBin.delta += LLONG_MAX;

            // reduce the overflow count (since underflow occurs)
            globalBin.deltaOverflowCount--;
        }
    }

    // offload the count and delta values - no need to worry about overflow as we already took care of it
    globalBin.count += currentHistBin.count;
    globalBin.delta += currentHistBin.delta;

    // reset the local bin
    currentHistBin.count = 0;
    currentHistBin.delta = 0;
    currentHistBin.tillOverFlowReset = MAX_DELTA_UPDATE_COUNT;
}

/**
 * @brief Update the bin with the given delta value
 * 
 * @tparam USE_LOCAL_HIST This will make sure to call appropriate overflow logic
 * @param localHist LocalBin vector - notice how we provide this regardless of the value of USE_LOCAL_HIST, if its false, we ignore it
 * @param globalHist Thread HistBin vector 
 * @param binId bin id as calculated by findBinIndex
 * @param delta delta value as calculated by findBinIndex
 */
template<bool USE_LOCAL_HIST>
inline void updateBin(std::vector<LocalBin> &localHist, std::vector<HistBin> &globalHist, int64 binId, int64 delta) {
    // handle the USE_LOCAL_HIST case - this is compile time check with constexpr
    if constexpr (USE_LOCAL_HIST) {
        // get the local bin
        LocalBin &currentHistBin = localHist[binId];

        // check if we need to offload the local bin to the global bin
        if (currentHistBin.tillOverFlowReset < 0) { // sign change here from == to <
            // get the global bin
            HistBin &globalBin = globalHist[binId];

            // offload the local bin to the global bin
            addToGlobalBin(currentHistBin, globalBin);
        }

        // update the local bin
        currentHistBin.count++;
        currentHistBin.delta += delta; // this will never overflow, since the counter will run out first

        // decrement the overflow reset counter
        currentHistBin.tillOverFlowReset--;
    } else // non-local hist case
    {
        // get the global bin
        HistBin &globalBin = globalHist[binId];

        // check if we need to update the overflow value
        if (globalBin.tillOverFlowReset < 0) { // sign change here from <= to <
            // update the overflow value
            updateDeltaOverflow(globalBin, delta);
        }

        // update the global bin
        globalBin.count++;
        globalBin.delta += delta;

        // decrement the overflow reset counter
        globalBin.tillOverFlowReset--;
    }
}


/**
 * @brief Offload the complete local PDF to the thread PDF
 * 
 * @param localPDFVector LocalBin vector to offload
 * @param histogram Thread HistBin vector 
 */
void offloadToGlobal(std::vector<LocalBin> &localPDFVector, std::vector<HistBin> &histogram) {
    // go through all the bins and offload them to the thread PDF 
    // notice how we don't have a parallel loop here, since we are already in a parallel region
    // we need to do this operation for all the threadPDFs individually
    for (int64 binId = 0; binId < localPDFVector.size(); binId++) {
        // get the local and global bins
        LocalBin &currentHistBin = localPDFVector[binId];
        HistBin &globalBin = histogram[binId];

        // offload the local bin to the global bin
        addToGlobalBin(currentHistBin, globalBin);
    }
}


/**
 * @brief The main sub-kernel to calculate the PDF over the j index
 * 
 * This is separated out to allow for better readability and maintainability. Same function can be used for both
 * the cell list and non-cell list cases with appropriate parameters.
 * 
 * @tparam USE_LOCAL_HIST To use LocalBin or not
 * @param positionsI Positions container 1
 * @param positionsJ Positions container 2
 * @param i Which index of the first positions to calculate the PDF for
 * @param j Start index of the second positions 
 * @param j_stop Stop index of the second positions
 * @param j_step Step size for the second positions
 * @param samePositions If both positions/cells are the same - does the half distance matrix calculation
 * @param idxStripMined strip mined array for bin index - only used here, but defined outside to avoid re-allocation
 * @param deltaStripMined strip mined array for delta - only used here, but defined outside to avoid re-allocation
 * @param localPDFVector local histogram to store the values before offloading to the global histogram - again ignored if USE_LOCAL_HIST is false
 * @param histogram Thread PDF to update
 */
template<bool USE_LOCAL_HIST>
inline void histOverJ(Positions const &positionsI, Positions const &positionsJ,
                      uint i, uint j, uint j_stop, uint j_step, bool samePositions,
                      std::array<int, N_STRIP_MINED> &idxStripMined,
                      std::array<int64, N_STRIP_MINED> &deltaStripMined,
                      std::vector<LocalBin> &localPDFVector,
                      std::vector<HistBin> &histogram) {
    // check if the positions are the same, if so, we only need to calculate the half distance matrix
    // set the stop index to the current i index 
    // same thing works in case of cell list as well.

    // update the scaling factors according to nBins
    int64 scaleUP = 1'000'000;
    float scaleDown = 0.001;
    // std::cout << "Size: " << histogram.size() << std::endl;
    if (histogram.size() == 173210) {
        // std::cout << "updating scaling" << std::endl;
        scaleUP = 1'000'000'00; // two extra zeros because of squared center
        scaleDown = 0.0001; // single extra zero because already scaled down
    }

    if (samePositions)
        j_stop = i;

    // define the strip mined variables
    uint j_strip_end, j_strip;

    // define the binId and delta - these will passed by reference to the findBinIndex function and updated there
    int64 binId, delta;

    // go through the j index
    for (j = j; j < j_stop; j += j_step) {
        // make sure that we don't go over the stop index
        j_strip_end = (j + N_STRIP_MINED < j_stop) ? j + N_STRIP_MINED : j_stop;

        // strip mine the loop
#pragma unroll
        for (j_strip = j; j_strip < j_strip_end; j_strip++) {
            // find the bin index and delta value
            findBinIndex(positionsI, positionsJ, i, j_strip, binId, delta, scaleDown, scaleUP);

            // store the values in the strip mined arrays
            idxStripMined[j_strip - j] = (int) binId;
            deltaStripMined[j_strip - j] = delta;
        }

        // and over the strip mined values, update the bins
#pragma unroll
        for (j_strip = j; j_strip < j_strip_end; j_strip++) {
            // get the binId and delta value
            binId = idxStripMined[j_strip - j];
            delta = deltaStripMined[j_strip - j];

#ifndef ROOFLINE
            // update the bin
            updateBin<USE_LOCAL_HIST>(localPDFVector, histogram, binId, delta);
#else
            histogram[i].count += 1;
            histogram[i].delta += delta;
#endif
        }
    }
}


#pragma endregion inline_functions


// main kernels
#pragma region kernels

/**
 * @brief PDF calculation kernel for non cell list case
 * 
 * @tparam USE_LOCAL_HIST whether to use LocalBin or not
 * @param positionsI Positions container 1
 * @param positionsJ Positions container 2
 * @param pdfVector thread PDF to update
 * @param start start index for i
 * @param stop stop index for i
 */
template<bool USE_LOCAL_HIST>
void calculatePDFKernel(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                        int64 start, int64 stop, bool samePositions) {

#ifdef ROOFLINE
    std::cout << "Roofline calculation\n";
    if (positionsI.size() > 1.7e6 || positionsJ.size() > 1.7e6){
        throw std::invalid_argument("n particles cannot be larger than NBins for roofline");
    }
#endif

// open a parallel region with reduction on the pdfVector
#pragma omp parallel reduction(+ : pdfVector) default(none) shared(positionsI, positionsJ, start, stop, samePositions)
    {
        // define the vars, j_stop is till the end of the second positions
        uint i, j_step = N_STRIP_MINED, j_stop = positionsJ.size();

        // define the strip mined arrays
        std::array<int, N_STRIP_MINED> idxStripMined{};
        std::array<int64, N_STRIP_MINED> deltaStripMined{};

        // define the local pdfVector
        std::vector<LocalBin> localPDFVector(pdfVector.size());

// parallel loop over the i index, dynamic scheduling works best
#pragma omp for schedule(dynamic)
        for (i = start; i < stop; i++) {
            // for each i, calculate the PDF over the j index
            histOverJ<USE_LOCAL_HIST>(positionsI, positionsJ, i, 0, j_stop, j_step, samePositions,
                                      idxStripMined, deltaStripMined, localPDFVector, pdfVector);
        }

        if constexpr (USE_LOCAL_HIST) {
            offloadToGlobal(localPDFVector, pdfVector);
        }
    }
}

/**
 * @brief PDF calculation kernel for cell list case
 * 
 * @tparam USE_LOCAL_HIST Whether to use LocalBin or not
 * @param positionsI Positions container 1
 * @param positionsJ Positions container 2
 * @param cellPairsList Cell pairs list (std::vector<std::array<int, 2>>), contains the cell pairs to calculate the PDF for
 * @param pdfVector Current thread PDF to update
 * @param start Start index for the cell pairs list
 * @param stop Stop index for the cell pairs list (step is always 1)
 */
template<bool USE_LOCAL_HIST>
void calculatePDFKernel(Positions const &positionsI, Positions const &positionsJ, PairsList const &cellPairsList,
                        std::vector<HistBin> &pdfVector, int64 start, int64 stop, bool samePositions) {

#ifdef ROOFLINE
    std::cout << "---------------- Roofline calculation -------------------\n";
    if (positionsI.size() > 1.7e6 || positionsJ.size() > 1.7e6){
        throw std::invalid_argument("n particles cannot be larger than NBins for roofline");
    }
#endif

// open a parallel region with reduction on the pdfVector
#pragma omp parallel reduction(+ : pdfVector) default(none) shared(positionsI, positionsJ, cellPairsList, start, stop, samePositions)
    {
        // define the vars
        uint i, pair_idx, j_step = N_STRIP_MINED;
        uint i_start, i_stop, j_start, j_stop;

        // define the strip mined arrays
        std::array<int, N_STRIP_MINED> idxStripMined{};
        std::array<int64, N_STRIP_MINED> deltaStripMined{};

        // define the local pdfVector
        std::vector<LocalBin> localPDFVector(pdfVector.size());

        // define flag to check if both cells are equal
        bool sameCell;


// parallel loop over the cell pairs list, dynamic scheduling works best again
#pragma omp for schedule(dynamic)
        for (pair_idx = start; pair_idx < stop; pair_idx += 1) {
            // get the pair of cells
            auto const &pair = cellPairsList[pair_idx];

            // check if both cells are equal
            // this is only required for cases where both positions are also the same
            // otherwise this will always be false
            sameCell = (samePositions && (pair[0] == pair[1]));

            // setup the start and end indices for the cells
            i_start = positionsI.cellHeads[pair[0]];
            i_stop = i_start + positionsI.atomCounts[pair[0]]; // go till however many atoms are there in the cell
            j_start = positionsJ.cellHeads[pair[1]];
            j_stop = j_start + positionsJ.atomCounts[pair[1]];

            // go through the i index for the first cell
            for (i = i_start; i < i_stop; i++) {
                // for each i, calculate the PDF over the j index
                histOverJ<USE_LOCAL_HIST>(positionsI, positionsJ, i, j_start, j_stop, j_step, sameCell,
                                          idxStripMined, deltaStripMined, localPDFVector, pdfVector);
            }
        }

        // finally offload the local pdfVector to the global pdfVector
        if constexpr (USE_LOCAL_HIST) {
            offloadToGlobal(localPDFVector, pdfVector);
        }
    }
}


void calculatePDFCPU(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                     const PairsList &sortedCellPairs, int64 start, int64 stop, bool samePositions,
                     ParallelHelper &parallelHelper, CalculationConfig &config) {
    if (config.useCellList) {
        start = 0;
        stop = sortedCellPairs.size(); // still need to deal with cellPairsList for multiple GPUs
        std::cout << ">> Thread(" << parallelHelper.worldRank << ") CellPairsList size: " << stop - start << "\n";
        // calculate the PDF - notice how we only pass the vector datastructure
        if (config.useLocalHistogram) {
            calculatePDFKernel<true>(positionsI, positionsJ,
                                     sortedCellPairs, pdfVector, start, stop, samePositions);
        } else {
            calculatePDFKernel<false>(positionsJ, positionsJ,
                                      sortedCellPairs, pdfVector, start, stop, samePositions);
        }
    } else {
        parallelHelper << ">> Non-cell list calculation\n";
        // non cell list based PDF calculation
        if (config.useLocalHistogram) {
            calculatePDFKernel<true>(positionsI, positionsJ,
                                     pdfVector, start, stop, samePositions);
        } else {
            calculatePDFKernel<false>(positionsI, positionsJ,
                                      pdfVector, start, stop, samePositions);
        }
    }
}


// not to be included 
// template <bool USE_LOCAL_HIST>
// void calculatePDFKernelSubSampled(std::vector<std::pair<Positions, Positions>> const &positionsPairs, 
//                                         std::vector<HistBin> &histogram,
//                                         int64 start, int64 stop, int64 step)
// {

//     LIKWID_MARKER_REGISTER("SUBSAMPLED");

// #pragma omp parallel reduction(+ : histogram)
//     {
//         LIKWID_MARKER_START("SUBSAMPLED");

//         std::array<int, N_STRIP_MINED> idxStripMined;
//         std::array<int64, N_STRIP_MINED> deltaStripMined;
//         std::vector<LocalBin> localPDFVector(N_BINS);

// #pragma omp for schedule(dynamic)
//         for (uint positionPairIndex=0; positionPairIndex < positionsPairs.size(); positionPairIndex++)
//         {
//             Positions const &positions1 = positionsPairs[positionPairIndex].first;
//             Positions const &positions2 = positionsPairs[positionPairIndex].second;

//             uint i, j_step = N_STRIP_MINED, j_stop = positions2.X.size();
//             bool bothEqual = &positions1 == &positions2;

//             for (i = start; i < stop; i += step)
//             {
//                 histOverJ<USE_LOCAL_HIST>(positions1, positions2, i, 0, j_stop, j_step, bothEqual, 
//                                                     idxStripMined, deltaStripMined, localPDFVector, histogram);
//             }
//         }
//     }
// }


#pragma endregion kernels