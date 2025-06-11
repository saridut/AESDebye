#pragma once
#include <DataTypes.hpp>
#include <PDF.hpp>
#include <ParallelHelper.hpp>
#include <CellList.hpp>
#include <PDFKernels.hpp>


#ifdef USE_GPU

void calculatePDFGPU(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                     const PairsList &sortedCellPairs, int64 start, int64 stop, bool samePositions,
                     ParallelHelper &parallelHelper, CalculationConfig &config);

#else

static void calculatePDFGPU(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                           const PairsList &sortedCellPairs, int64 start, int64 stop, bool samePositions,
                           ParallelHelper &parallelHelper, CalculationConfig &config)
{

    throw std::runtime_error("GPU support is not enabled in the current build. Please recompile with GPU support.");
}

#endif