
#pragma once

#include <CellList.hpp>

#pragma region pairs_reduction

// custom reduction for pairIndexHistogram per bin first
// only used the in the user defined reduction for the pairsHistogram
[[maybe_unused]] std::vector<int> addPairIndexBin(std::vector<int> a, std::vector<int> b)
{
    std::vector<int> result(a.size() + b.size());
    for (int i = 0; i < a.size(); i++)
    {
        result[i] = a[i];
    }

    for (int i = 0; i < b.size(); i++)
    {
        result[i + a.size()] = b[i];
    }

    return result;
}

// user defined reduction for std::vector<std::vector<int>> for the pairsHistogram
#pragma omp declare reduction(accumulate:std::vector<std::vector<int>> : \
std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), addPairIndexBin)) \
initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

#pragma endregion pairs_reduction

void CellList::createCellList(Positions &positions) const
{
    if (positions.cellHeads.size() == nTotalCells)
    {
        std::cout << "Cell list already created\n";
        return;
    } // already created the cell list

    int64 cellSize = maxBoxSize / nCells;
    cellSize += 1;

    // Initialize a cell list and sorted positions
    Positions sortedPositions(positions.boxSize);
    std::vector<std::vector<int>> cellList(nTotalCells, std::vector<int>());

    // put atoms into the cells
    for (int i = 0; i < positions.X.size(); i++)
    {
        auto cellX = (int) (positions.X[i] / cellSize);
        auto cellY = (int) (positions.Y[i] / cellSize);
        auto cellZ = (int) (positions.Z[i] / cellSize);
        auto cellIdx = (int) (cellZ + cellY * nCells + cellX * nCells * nCells);
        cellList[cellIdx].push_back(i);
        assert(cellX < nCells);
        assert(cellY < nCells);
        assert(cellZ < nCells);
    }

    // initialize the cell heads and atom counts for flattening
    auto& cellHeads = sortedPositions.cellHeads; // use the sortedPositions to store the cellHeads - NOT the positions
    auto& atomCounts = sortedPositions.atomCounts; // use the sortedPositions to store the atomCounts 
    cellHeads.resize(nTotalCells, -1);
    atomCounts.resize(nTotalCells, 0);

    // put the atoms in the sorted positions
    int currentHeadPosition = 0;
    for (int cellIdx = 0; cellIdx < nTotalCells; cellIdx++)
    {
        atomCounts[cellIdx] = (int) cellList[cellIdx].size();
        if (atomCounts[cellIdx] == 0)
            continue;
        cellHeads[cellIdx] = currentHeadPosition;

        for (int atomicId : cellList[cellIdx])
        {
            sortedPositions.push_back(positions, atomicId);
            currentHeadPosition++;
        }
    }
    sortedPositions.element = positions.element;
    positions = sortedPositions; // will update the positions, add the cellHeads and atomCounts to the positions as well
}

void CellList::createCellPairs(bool sorted)
{
    // SETUP: create cell positions and cell pairs histogram
    // make sure the number of cells is less than 30, otherwise the number of pairs will be too large and
    // create memory issues!. Also, 30 cell pairs are more than sufficient almost always!
    if (nCells > 30)
    {
        throw std::range_error("Cannot create pairs list for more than 30 cells, reduce the number of cells!");
    }

    for (int i = 0; i < nCells; i++)
    {
        for (int j = 0; j < nCells; j++)
        {
            for (int k = 0; k < nCells; k++)
            {
                cellPositions.push_back({i, j, k});
            }
        }
    }
    auto maxDistance2 = (float) (cellPositions[0].distance2(cellPositions[nTotalCells - 1]) + 0.1);

    // number of bins for the histogram of cell pairs: arbitrarily chosen, should be enough for most cases
    size_t nPairsHistBins = 10000;
    auto scalingFactor = (float) nPairsHistBins / maxDistance2;
    std::vector<std::array<int, 2>> cellPairs(nTotalCells * nTotalCells, {-1, -1});
    std::vector<std::vector<int>> pairsHistogram(nPairsHistBins);  
    for (auto v : pairsHistogram)
    {
        v.reserve(10000); // reserve some space for the pairs to avoid reallocation during push_back
    }

    // Main: Calculate the distances between all cells
#pragma omp parallel for reduction(accumulate : pairsHistogram) if (nTotalCells > 50000)
    for (int i = 0; i < cellPositions.size(); i++)
    {
        for (int j = 0; j <= i; j++)         // only need one side of the matrix with diagonal
        {
            cellPairs[i * nTotalCells + j] = {i, j};
            if (!sorted)
                continue;
            float distance = cellPositions[i].distance2(cellPositions[j]);
            int binIndex = distance * scalingFactor;
            pairsHistogram[binIndex].push_back(i * nTotalCells + j);
        }
    }

    // create sorted pairs list
    if (sorted) // use the pairsHistogram to create the sortedPairsList
    {
        for (int i = 0; i < pairsHistogram.size(); i++)
        {
            for (int j = 0; j < pairsHistogram[i].size(); j++)
            {
                int pairIndex = pairsHistogram[i][j];
                auto pair = cellPairs[pairIndex];
                sortedPairsList.push_back(pair);
            }
        }
    }
    else // use the cellPairs to create the sortedPairsList
    {
        for (int i = 0; i < cellPairs.size(); i++)
        {
            if (cellPairs[i][0] == -1) // skip empty pairs
                continue;
            sortedPairsList.push_back(cellPairs[i]);
        }
    }
};

PairsList CellList::getCellPairsList(Positions const &positionsI, Positions const &positionsJ, bool samePositions)
{
    PairsList pairsList;
    if (!samePositions){
        pairsList.reserve(nTotalCells * nTotalCells); // reserve space for the pairs
        for (auto pair : sortedPairsList)
        {
            pairsList.push_back(pair);
            if (pair[0] != pair[1]) // self-pairs don't need to be added twice
                pairsList.push_back({pair[1], pair[0]});   // add reverse pairs as well (i, j) and (j, i)
        }
    }
    else{
        pairsList = sortedPairsList; // otherwise just copy the sortedPairsList
    }

    // filter out pairs with empty cells
    PairsList filteredPairsList;
    for (auto pair: pairsList)
    {
        int pairIndexI = pair[0];
        int pairIndexJ = pair[1];
        if (positionsI.atomCounts[pairIndexI] == 0 || positionsJ.atomCounts[pairIndexJ] == 0)
            continue; 
        filteredPairsList.push_back(pair);
    }
    return filteredPairsList;
}


std::tuple<std::vector<int>, std::vector<int>> 
    CellList::setupRankBasedDivisions(Positions const &positions1, Positions const &positions2,
                            PairsList const &pairsList, int nRanks, bool loadBalancing)
{
    std::vector<int> rankStarts;
    std::vector<int> rankEnds;
    if (!loadBalancing || nRanks == 1 ) // if no load balancing or only one rank
    {    
        rankStarts.resize(nRanks, 0);
        rankEnds.resize(nRanks, 0);
        for (int i = 0; i < nRanks; i++)
        {
            rankStarts[i] = i * pairsList.size() / nRanks;
            rankEnds[i] = (i + 1) * pairsList.size() / nRanks;
        }

        return {rankStarts, rankEnds};
    }

    std::vector<int64> atomicPairCounts(pairsList.size(), 0);
    bool same = (&positions1 == &positions2);

    #pragma omp parallel for schedule(static)
    for (uint i = 0; i < pairsList.size(); i++)
    {
        auto pair = pairsList[i];
        int pairIndexI = pair[0];
        int pairIndexJ = pair[1];
        if (positions1.atomCounts[pairIndexI] == 0 || positions2.atomCounts[pairIndexJ] == 0)
        {
            atomicPairCounts[i] = 0;
            continue;
        }

        if ((pairIndexI == pairIndexJ) && same)
        {
            atomicPairCounts[i] = positions1.atomCounts[pairIndexI] * (positions1.atomCounts[pairIndexI] - 1) / 2;
        }
        else
        {
            atomicPairCounts[i] = positions1.atomCounts[pairIndexI] * positions2.atomCounts[pairIndexJ];
        }
    }

    // total pairCounts
    int64 totalAtomicPairs;
    if (same)
    {
        totalAtomicPairs = (int64) positions1.size() * ( (int64) positions1.size() - 1) / 2;
    }
    else
    {
        totalAtomicPairs = (int64) positions1.size() * (int64) positions2.size();
    }

    // efficient partitioning - based on pair counts
    int64 nPairsTillNow = 0;
    int64 nPairsPerRank = totalAtomicPairs / nRanks;

    nPairsPerRank += (totalAtomicPairs % nRanks) ? 1 : 0;
    rankStarts.push_back(0); // first rank starts at 0
    int64 pairsIndex = 0;
    for (int i = 1; i < nRanks; i++)
    {
        int64 nPairsTillNow = 0;
        while (nPairsTillNow < nPairsPerRank && pairsIndex < atomicPairCounts.size())
        {
            nPairsTillNow += atomicPairCounts[pairsIndex];
            pairsIndex++;
        }
        rankEnds.push_back(pairsIndex);
        rankStarts.push_back(pairsIndex);
    }
    rankEnds.push_back(pairsList.size());

    assert(rankStarts.size() == nRanks);
    assert(rankEnds.size() == nRanks);

    return {rankStarts, rankEnds};
}

void CellList::createFlattenedPairsList()
{
    flattendPairsList = std::vector<int>(sortedPairsList.size() * 2);
    for (int i = 0; i < sortedPairsList.size(); i++)
    {
        auto pair = sortedPairsList[i];
        flattendPairsList[2 * i] = pair[0];
        flattendPairsList[2 * i + 1] = pair[1];
    }
}

std::tuple<uint, uint> CellList::getMinMaxBin(Cell cellI, Cell cellJ){
    auto [minDistance, maxDistance] = getMinMaxDistance(cellI, cellJ);

    float longestDiagonal = std::sqrt(3) * nCells;
    float eps = std::numeric_limits<float>::epsilon();
    uint minBin = uint((minDistance - eps) * (float) N_BINS / longestDiagonal);
    uint maxBin = uint((maxDistance + eps) * (float) N_BINS / longestDiagonal);
    minBin = std::min(minBin, uint(N_BINS - 1));
    maxBin = std::min(maxBin, uint(N_BINS - 1));
    return {minBin, maxBin};
}


std::tuple<float, float> CellList::getMinMaxDistance(Cell cellI, Cell cellJ)
{
    assert (cellI < nCells);
    assert (cellJ < nCells);

    float dx = std::abs(cellI.x - cellJ.x);
    float dy = std::abs(cellI.y - cellJ.y);
    float dz = std::abs(cellI.z - cellJ.z);

    float dx_active = (dx > 0);
    float dy_active = (dy > 0);
    float dz_active = (dz > 0);

    // only use the active ones for the min distance, all for the max distance
    float minDistance = std::sqrt((dx-1)*(dx-1)*dx_active + (dy-1)*(dy-1)*dy_active + (dz-1)*(dz-1)*dz_active); 
    float maxDistance = std::sqrt((dx+1)*(dx+1) + (dy+1)*(dy+1) + (dz+1)*(dz+1)); // always works

    return {minDistance, maxDistance};
}

std::tuple<float, float> CellList::getMinMaxDistance(int i, int j){
    return getMinMaxDistance(get3DIndex(i), get3DIndex(j));
}

std::tuple<uint, uint> CellList::getMinMaxBin(int i, int j){
    return getMinMaxBin(get3DIndex(i), get3DIndex(j));
}

int CellList::get1DIndex(int i, int j, int k)
{
    return k + j * nCells + i * nCells * nCells;
}

std::tuple<int, int, int> CellList::get3DIndex(int idx)
{
    int k = idx % nCells;
    int j = (idx / nCells) % nCells;
    int i = idx / (nCells * nCells);
    return {i, j, k};
}




        // // same cell
        // if (dx == 0 && dy == 0 && dz == 0)
        // {
        //     return {0, maxDistance};
        // }

        // // if two of dx, dy, dz are zero, then they are in the same row or column or height in 3d grid
        // // min distance then becomes the non-zero distance, max distance would by sqrt((non_zero_distance + 1)^2 + 1)
        // // same rows/columns
        // if (dx == 0 && dy == 0)
        // {
        //     return {dz - 1, maxDistance};
        // }
        // if (dx == 0 && dz == 0)
        // {
        //     return {dy - 1, maxDistance};
        // }
        // if (dy == 0 && dz == 0)
        // {
        //     return {dx - 1, maxDistance};
        // }


        // if (dx == 0)
        // {
        //     return {std::sqrt((dy-1)*(dy-1) + (dz-1)*(dz-1)), maxDistance};
        // }
        // if (dy == 0)
        // {
        //     return {std::sqrt((dx-1)*(dx-1) + (dz-1)*(dz-1)), maxDistance};
        // }
        // if (dz == 0)
        // {
        //     return {std::sqrt((dx-1)*(dx-1) + (dy-1)*(dy-1)), maxDistance};
        // }

// #pragma once

// struct Pair
// {
//     int i, j;
//     double distance;
//     int64 numPairs;

//     // overload the << operator to print the pair
//     friend std::ostream &operator<<(std::ostream &os, const Pair &pair)
//     {
//         os << "Pair: " << pair.i << " " << pair.j << " " << pair.distance;
//         return os;
//     }

//     static bool comparePair(Pair a, Pair b)
//     {
//         if (a.distance == b.distance)
//         {
//             if (a.i == b.i)
//             {
//                 return a.j < b.j;
//             }
//             return a.i < b.i;
//         }
//         return a.distance < b.distance;
//     }
// };

// bool comparePair(Pair a, Pair b);
// #include "dataTypes.h"
// // #include "utils.h"

// // Define the cell structure
// struct Pair {
//     int i, j;
//     double distance;
// };

// bool comparePair(Pair a, Pair b) {
//     if (a.distance == b.distance) {
//         if (a.i == b.i) {
//             return a.j < b.j;
//         }
//         return a.i < b.i;
//     }
//     return a.distance < b.distance;
// }

// // Create the cell list class
// class CellList {
// public:
//     std::vector<std::vector<int>> cellList;
//     std::vector<int> cellHeads;
//     std::vector<int> atomCounts;

//     Positions sortedPositions;
//     std::vector<std::array<int, 2>> sortedPairsList;

//     int64_t maxBoxSize = 1'000'000'000;

// public:
//     CellList(Positions positions, int nCells) {
//         int64_t cellSize = maxBoxSize / nCells;
//         cellSize += 1;

//         int nCells = nCells, nCells = nCells, nCells = nCells;
//         int totalCells = nCells * nCells * nCells;

//         // Initialize the cell list
//         cellList.resize(totalCells, std::vector<int>());

//         // put atoms into the cells
//         for (int i = 0; i < positions.X.size(); i++) {
//             int cellX = positions.X[i] / cellSize;
//             int cellY = positions.Y[i] / cellSize;
//             int cellZ = positions.Z[i] / cellSize;

//             int cellIdx = cellZ + cellY * nCells + cellX * nCells * nCells;
//             cellList[cellIdx].push_back(i);
//             assert (cellX < nCells);
//             assert (cellY < nCells);
//             assert (cellZ < nCells);
//         }

//         // initialize the cell heads and atom counts for flattening
//         cellHeads.resize(totalCells, -1);
//         atomCounts.resize(totalCells, 0);
//         int currentHeadPosition = 0;
//         for (int cellIdx=0; cellIdx < totalCells; cellIdx++) {
//             atomCounts[cellIdx] = cellList[cellIdx].size();
//             if (atomCounts[cellIdx] == 0) continue;
//             cellHeads[cellIdx] = currentHeadPosition;

//             for (int atomicId : cellList[cellIdx]) {
//                 sortedPositions.X.push_back(positions.X[atomicId]);
//                 sortedPositions.Y.push_back(positions.Y[atomicId]);
//                 sortedPositions.Z.push_back(positions.Z[atomicId]);
//                 currentHeadPosition++;
//             }
//         }

//         createSortedPairsList(nCells, nCells, nCells);
//     }

//     void createSortedPairsList(int nCells, int nCells, int nCells) {

//         Positions cellPositions;

//         // Create meshgrid of cell positions
//         for (int i = 0; i < nCells; i++) {
//             for (int j = 0; j < nCells; j++) {
//                 for (int k = 0; k < nCells; k++) {
//                     cellPositions.X.push_back(i);
//                     cellPositions.Y.push_back(j);
//                     cellPositions.Z.push_back(k);
//                 }
//             }
//         }

//         std::vector<Pair> pairs;

//         // Calculate the distances between all cells
//         for (int i = 0; i < cellPositions.X.size(); i++) {
//             for (int j = 0; j < cellPositions.X.size(); j++) {
//                 if (i > j) continue; // only need one side of the matrix with diagonal

//                 if (atomCounts[i] == 0 || atomCounts[j] == 0) continue; // if either cell is empty, skip (no pairs to be made

//                 // int cellIdxI = cellPositions.Z[i] + cellPositions.Y[i] * nCells + cellPositions.X[i] * nCells * nCells;
//                 // int cellIdxJ = cellPositions.Z[j] + cellPositions.Y[j] * nCells + cellPositions.X[j] * nCells * nCells;
//                 // if (atomCounts[cellIdxI] == 0 || atomCounts[cellIdxJ] == 0) continue; // if either cell is empty, skip (no pairs to be made

//                 double distance = std::sqrt(std::pow(cellPositions.X[i] - cellPositions.X[j], 2) + std::pow(cellPositions.Y[i] - cellPositions.Y[j], 2) + std::pow(cellPositions.Z[i] - cellPositions.Z[j], 2));
//                 pairs.push_back({i, j, distance});
//             }
//         }

//         // Sort the pairs - logic in the Pair struct
//         std::sort(pairs.begin(), pairs.end(), comparePair);
//         for (int i = 0; i < pairs.size(); i++) {
//             sortedPairsList.push_back({pairs[i].i, pairs[i].j});
//         }
//         std::cout << "Size of the pair list: " << sortedPairsList.size() << "\n";
//     }

//     std::vector<int> createFlattenedPairsList(){
//         std::vector<int> flattendPairsList(sortedPairsList.size()*2);
//         for (int i = 0; i < sortedPairsList.size(); i++)
//         {
//             auto pair = sortedPairsList[i];
//             flattendPairsList[2*i] = pair[0];
//             flattendPairsList[2*i + 1] = pair[1];
//         }
//         return flattendPairsList;
//     }

// };

// // int main() {
// //     // Example usage
// //     int nRepeats = 18;
// //     double lattice = 3.89070;
// //     double noise = 0.1;
// //     std::string filename = "../data/test_dir/test.xyz";
// //     generateData(nRepeats, lattice, filename, noise);

// //     auto [positions, boxSize] = readXYZ(filename);
// //     std::cout << positions.X.size() << std::endl;

// //     int64 minX, minY, minZ, maxX, maxY, maxZ;

// //     for (int i = 0; i < positions.X.size(); i++) {
// //         if (i == 0) {
// //             minX = positions.X[i];
// //             minY = positions.Y[i];
// //             minZ = positions.Z[i];
// //             maxX = positions.X[i];
// //             maxY = positions.Y[i];
// //             maxZ = positions.Z[i];
// //         } else {
// //             if (positions.X[i] < minX) {
// //                 minX = positions.X[i];
// //             }
// //             if (positions.Y[i] < minY) {
// //                 minY = positions.Y[i];
// //             }
// //             if (positions.Z[i] < minZ) {
// //                 minZ = positions.Z[i];
// //             }
// //             if (positions.X[i] > maxX) {
// //                 maxX = positions.X[i];
// //             }
// //             if (positions.Y[i] > maxY) {
// //                 maxY = positions.Y[i];
// //             }
// //             if (positions.Z[i] > maxZ) {
// //                 maxZ = positions.Z[i];
// //             }
// //         }
// //     }

// //     int64 boxSize = std::max({maxX - minX, maxY - minY, maxZ - minZ});
// //     std::cout << "Max X: " << maxX << " Min X: " << minX << std::endl;

// //     std::cout << "Box size: " << boxSize / 10000000 << std::endl;
// //     int64 cellSize = 10000000;

// //     std::cout << "positions: " << positions.X.size() << std::endl;
// //     CellList cellList(positions, 10);
// //     auto sortedPositions = cellList.sortedPositions;

// //     auto nAtoms = cellList.atomCounts;

// //     // int sum = 0;
// //     // for (int i = 0; i < nAtoms.size(); i++) {
// //     //     if (nAtoms[i] > 1)
// //     //         std::cout << nAtoms[i] << "\n";
// //     //     sum += nAtoms[i];
// //     // }

// //     // std::cout << "size: " << nAtoms.size() << std::endl;
// //     // std::cout << "Sum: " << sum << std::endl;

// //     // cellList.findAllEquivalentCells();

// //     std::cout << "Size of sorted pairs: " << cellList.sortedPairs.size() << std::endl;
// //     return 0;
// // }

// // Create the cell list class
// // class CellList {
// // public:
// //     std::vector<std::vector<int>> cellList;
// //     std::vector<int> cellHeads;
// //     std::vector<int> atomCounts;

// //     Positions sortedPositions;
// //     std::vector<std::array<int, 2>> sortedPairsList;

// //     int64_t maxBoxSize = 1'000'000'000;

// // public:

// // };

// // int main() {
// //     // Example usage
// //     int nRepeats = 18;
// //     double lattice = 3.89070;
// //     double noise = 0.1;
// //     std::string filename = "../data/test_dir/test.xyz";
// //     generateData(nRepeats, lattice, filename, noise);

// //     auto [positions, boxSize] = readXYZ(filename);
// //     std::cout << positions.X.size() << std::endl;

// //     int64 minX, minY, minZ, maxX, maxY, maxZ;

// //     for (int i = 0; i < positions.X.size(); i++) {
// //         if (i == 0) {
// //             minX = positions.X[i];
// //             minY = positions.Y[i];
// //             minZ = positions.Z[i];
// //             maxX = positions.X[i];
// //             maxY = positions.Y[i];
// //             maxZ = positions.Z[i];
// //         } else {
// //             if (positions.X[i] < minX) {
// //                 minX = positions.X[i];
// //             }
// //             if (positions.Y[i] < minY) {
// //                 minY = positions.Y[i];
// //             }
// //             if (positions.Z[i] < minZ) {
// //                 minZ = positions.Z[i];
// //             }
// //             if (positions.X[i] > maxX) {
// //                 maxX = positions.X[i];
// //             }
// //             if (positions.Y[i] > maxY) {
// //                 maxY = positions.Y[i];
// //             }
// //             if (positions.Z[i] > maxZ) {
// //                 maxZ = positions.Z[i];
// //             }
// //         }
// //     }

// //     int64 boxSize = std::max({maxX - minX, maxY - minY, maxZ - minZ});
// //     std::cout << "Max X: " << maxX << " Min X: " << minX << std::endl;

// //     std::cout << "Box size: " << boxSize / 10000000 << std::endl;
// //     int64 cellSize = 10000000;

// //     std::cout << "positions: " << positions.X.size() << std::endl;
// //     CellList cellList(positions, 10);
// //     auto sortedPositions = cellList.sortedPositions;

// //     auto nAtoms = cellList.atomCounts;

// //     // int sum = 0;
// //     // for (int i = 0; i < nAtoms.size(); i++) {
// //     //     if (nAtoms[i] > 1)
// //     //         std::cout << nAtoms[i] << "\n";
// //     //     sum += nAtoms[i];
// //     // }

// //     // std::cout << "size: " << nAtoms.size() << std::endl;
// //     // std::cout << "Sum: " << sum << std::endl;

// //     // cellList.findAllEquivalentCells();

// //     std::cout << "Size of sorted pairs: " << cellList.sortedPairs.size() << std::endl;
// //     return 0;
// // }
