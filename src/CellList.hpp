/**
 * @file CellList.hpp
 * @brief Contains the declaration of the Cell and CellList classes for spatial partitioning.
 */

#pragma once

#include <DataTypes.hpp>
#include <Positions.hpp>

typedef std::vector<std::array<int, 2>> PairsList;

/**
 * @class Cell
 * @brief Represents a single 3D cell in the spatial partitioning grid.
 */
class Cell
{
public:
    float x{0}, y{0}, z{0}; //< Coordinate of the cell in 3D grid space.
    
    /**
     * @brief Default constructor. Initializes coordinates to 0.
     */
    Cell() {}

    /**
     * @brief Construct a Cell with explicit coordinates.
     */
    Cell(int x, int y, int z) : x((float) x), y((float) y), z((float) z) {}

    /**
     * @brief Construct a Cell from a 3D index tuple.
     */
    Cell(std::tuple<int, int, int> idx) : x((float) std::get<0>(idx)), y((float) std::get<1>(idx)), z((float) std::get<2>(idx)) {}

    /**
     * @brief Calculate the Euclidean distance to another cell.
     * @param cell The target cell.
     * @return float The Euclidean distance.
     */
    float distance(Cell cell)
    {
        return std::sqrt(distance2(cell));
    }

    /**
     * @brief Calculate the squared Euclidean distance to another cell.
     * @param cell The target cell.
     * @return float The squared Euclidean distance.
     */
    float distance2(Cell cell)
    {
        float dx = x - cell.x;
        float dy = y - cell.y;
        float dz = z - cell.z;
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * @brief Less-than operator comparison against a scalar threshold.
     */
    bool operator<(const float &rhs) const
    {
        return x < rhs && y < rhs && z < rhs;
    }
};

/**
 * @class CellList
 * @brief Manages cell lists for spatial decomposition of particles to speed up distance calculations.
 */
class CellList
{
private:
    PairsList sortedPairsList;
    std::vector<std::vector<int>> cellList;
    std::vector<Cell> cellPositions;
    uint nCells;
    uint nTotalCells;
    std::vector<int> flattendPairsList;
    int64 maxBoxSize = 1'000'000'000;

public:
    /**
     * @brief Default constructor.
     */
    CellList() {};

    /**
     * @brief Construct a CellList with a specified grid size.
     * @param nCells Number of cells along one dimension.
     */
    CellList(uint nCells): nCells(nCells), nTotalCells(nCells * nCells * nCells) {};

    /**
     * @brief Populates the cell list with the given positions.
     * @param positions Reference to the positions of the atoms.
     */
    void createCellList(Positions &positions) const;

    /**
     * @brief Generates all cell pairs.
     * @param sorted Flag indicating if the cell pairs should be sorted by distance.
     */
    void createCellPairs(bool sorted=true);

    /**
     * @brief Generates a flattened representation of the cell pairs list for optimized access.
     */
    void createFlattenedPairsList();

    /**
     * @brief Retrieve cell pairs list for two sets of positions.
     */
    PairsList getCellPairsList(Positions const &positionsI, Positions const &positionsJ, bool samePositions);

    /**
     * @brief Sets up divisions of work based on rank for MPI execution.
     */
    std::tuple<std::vector<int>, std::vector<int>> 
    setupRankBasedDivisions(Positions const &positions1, Positions const &positions2,
                            PairsList const &pairsList, int nRanks, bool loadBalancing=true);

    /**
     * @brief Get the minimum and maximum bin indices corresponding to distance between two cells.
     */
    std::tuple<uint, uint> getMinMaxBin(Cell cellI, Cell cellJ);

    /**
     * @brief Get the minimum and maximum distance range between two cells.
     */
    std::tuple<float, float> getMinMaxDistance(Cell cellI, Cell cellJ);

    /**
     * @brief Get the minimum and maximum distance range between two cells identified by their indices.
     */
    std::tuple<float, float> getMinMaxDistance(int i, int j);

    /**
     * @brief Get the minimum and maximum bin range between two cells identified by their indices.
     */
    std::tuple<uint, uint> getMinMaxBin(int i, int j);

    /**
     * @brief Map 3D grid indices to a 1D flat index.
     */
    int get1DIndex(int i, int j, int k);

    /**
     * @brief Map a 1D flat index back to 3D grid indices.
     */
    std::tuple<int, int, int> get3DIndex(int idx);
};
