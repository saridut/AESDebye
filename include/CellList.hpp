
#pragma once

#include <DataTypes.hpp>
#include <Positions.hpp>

typedef std::vector<std::array<int, 2>> PairsList;

class Cell
{
public:
    float x{0}, y{0}, z{0};
    Cell() {}
    Cell(int x, int y, int z) : x((float) x), y((float) y), z((float) z) {}
    Cell(std::tuple<int, int, int> idx) : x((float) std::get<0>(idx)), y((float) std::get<1>(idx)), z((float) std::get<2>(idx)) {}

    float distance(Cell cell)
    {
        return std::sqrt(distance2(cell));
    }

    float distance2(Cell cell)
    {
        float dx = x - cell.x;
        float dy = y - cell.y;
        float dz = z - cell.z;
        return dx * dx + dy * dy + dz * dz;
    }

    // less than a scaler overload
    bool operator<(const float &rhs) const
    {
        return x < rhs && y < rhs && z < rhs;
    }
};



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
    CellList() {};
    CellList(uint nCells): nCells(nCells), nTotalCells(nCells * nCells * nCells) {};
    void createCellList(Positions &positions) const;
    void createCellPairs(bool sorted=true);
    void createFlattenedPairsList();
    PairsList getCellPairsList(Positions const &positionsI, Positions const &positionsJ, bool samePositions);
    std::tuple<std::vector<int>, std::vector<int>> 
    setupRankBasedDivisions(Positions const &positions1, Positions const &positions2,
                            PairsList const &pairsList, int nRanks, bool loadBalancing=true);

    std::tuple<uint, uint> getMinMaxBin(Cell cellI, Cell cellJ);
    std::tuple<float, float> getMinMaxDistance(Cell cellI, Cell cellJ);
    std::tuple<float, float> getMinMaxDistance(int i, int j);
    std::tuple<uint, uint> getMinMaxBin(int i, int j);
    int get1DIndex(int i, int j, int k);
    std::tuple<int, int, int> get3DIndex(int idx);
};

