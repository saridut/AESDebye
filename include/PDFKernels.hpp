#pragma once

#include <CellList.hpp>
#include <DataTypes.hpp>
#include <PDF.hpp>
#include <Positions.hpp>
#include <ParallelHelper.hpp>

template <bool USE_LOCAL_HIST>
void calculatePDFKernel(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                        int64 start, int64 stop, bool samePositions);

template <bool USE_LOCAL_HIST>
void calculatePDFKernel(Positions const &positionsI, Positions const &positionsJ, PairsList const &cellPairsList,
                        std::vector<HistBin> &pdfVector, int64 start, int64 stop, bool samePositions);



void calculatePDFCPU(Positions const &positionsI, Positions const &positionsJ, std::vector<HistBin> &pdfVector,
                            const PairsList &sortedCellPairs, int64 start, int64 stop, bool samePositions,
                            ParallelHelper &parallelHelper, CalculationConfig &config);
