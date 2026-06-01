#include <Positions.hpp>

Positions::Positions(std::vector<std::string> chemicalSymbols,
                     std::vector<std::array<double, 3>> positions,
                     std::vector<std::string> selectionIds,
                     double boxMin,
                     double boxMax)
{
    if (boxMin == 1e100 && boxMax == -1e100)
    {
        for (auto & position : positions)
        {
            double x = position[0];
            double y = position[1];
            double z = position[2];
            boxMin = std::min({boxMin, x, y, z});
            boxMax = std::max({boxMax, x, y, z});
        }
    }

    if (selectionIds.empty())
    {
        selectionIds = std::vector<std::string>(positions.size(), "0");
    }

    if (chemicalSymbols.size() != positions.size())
    {
        throw std::runtime_error("Number of chemical symbols does not match the number of positions.");
    }

    if (boxMin > boxMax)
    {
        throw std::runtime_error("Box minimum is greater than box maximum.");
    }

    this->chemicalSymbols = std::move(chemicalSymbols);
    this->selectionIds = selectionIds;
    boxSize = boxMax - boxMin + 1e-6;
    resize(positions.size());
    for (int i = 0; i < positions.size(); ++i)
    {
        // check if the positions are within the box
        double x = positions[i][0];
        double y = positions[i][1];
        double z = positions[i][2];
        if (x < boxMin || x > boxMax)
        {
            throw std::runtime_error("Position x-coordinate (" + std::to_string(x) + ") is outside the box boundaries (Min: " + std::to_string(boxMin) + " - Max:" + std::to_string(boxMax) + ").");
        }
        if (y < boxMin || y > boxMax)
        {
            throw std::runtime_error("Position y-coordinate (" + std::to_string(y) + ") is outside the box boundaries (Min: " + std::to_string(boxMin) + " - Max:" + std::to_string(boxMax) + ").");
        }
        if (z < boxMin || z > boxMax)
        {
            throw std::runtime_error("Position z-coordinate (" + std::to_string(z) + ") is outside the box boundaries (Min: " + std::to_string(boxMin) + " - Max:" + std::to_string(boxMax) + ").");
        }
        Xd[i] = positions[i][0] - boxMin;
        Yd[i] = positions[i][1] - boxMin;
        Zd[i] = positions[i][2] - boxMin;
    }

    scaleUp();
}

Positions Positions::subsample(double subsampleRatio, int seed)
{
    uint nTotalAtomsNew = std::floor(size() * subsampleRatio);
    std::vector<int> indices(size());
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 g(seed == -1 ? std::random_device()() : seed);
    indices = std::vector<int>(indices.begin(), indices.begin() + nTotalAtomsNew);
    return fromIndexList(*this, indices);
}


void Positions::push_back(const Positions &positions, uint index)
{
    X.push_back(positions.X[index]);
    Y.push_back(positions.Y[index]);
    Z.push_back(positions.Z[index]);

    Xd.push_back(positions.Xd[index]);
    Yd.push_back(positions.Yd[index]);
    Zd.push_back(positions.Zd[index]);

    chemicalSymbols.push_back(positions.chemicalSymbols[index]);
    selectionIds.push_back(positions.selectionIds[index]);
}


Positions Positions::fromIndexList(const Positions &positions, const std::vector<int> &indexList)
{
    Positions newPositions(positions.boxSize);
    for (int i = 0; i < indexList.size(); i++)
    {
        newPositions.push_back(positions, indexList[i]);
    }
    newPositions.element = positions.element;
    return newPositions;
}

Positions Positions::filterBySelectionId(std::string atomId)
{
    std::vector<int> indices;
    for (int i = 0; i < size(); i++)
    {
        if (selectionIds[i] == atomId)
        {
            indices.push_back(i);
        }
    }
    return fromIndexList(*this, indices);
}

Positions Positions::filterByElement(std::string atomType)
{
    std::vector<int> indices;
    for (int i = 0; i < size(); i++)
    {
        if (chemicalSymbols[i] == atomType)
        {
            indices.push_back(i);
        }
    }
    auto positions = fromIndexList(*this, indices);
    positions.element = atomType;
    return positions;
}

std::vector<std::string> Positions::getUniqueElements() const
{
    std::vector<std::string> uniqueAtomTypes;
    for (int i = 0; i < size(); i++)
    {
        if (std::find(uniqueAtomTypes.begin(), uniqueAtomTypes.end(), chemicalSymbols[i]) == uniqueAtomTypes.end())
        {
            uniqueAtomTypes.push_back(chemicalSymbols[i]);
        }
    }
    return uniqueAtomTypes;
}

std::vector<std::string> Positions::getUniqueSelectionIds() const
{
    std::vector<std::string> uniqueAtomIds;
    for (int i = 0; i < size(); i++)
    {
        if (std::find(uniqueAtomIds.begin(), uniqueAtomIds.end(), selectionIds[i]) == uniqueAtomIds.end())
        {
            uniqueAtomIds.push_back(selectionIds[i]);
        }
    }
    return uniqueAtomIds;
}

void Positions::toXYZ(std::string filename)
{
    std::ofstream file(filename);
    file << size() << "\n";
    file << "Generated by Debye\n";
    for (int i = 0; i < size(); i++)
    {
        file << std::setprecision(10) << std::fixed << chemicalSymbols[i] << " " << Xd[i] << " " << Yd[i] << " " << Zd[i] << " " << selectionIds[i] << "\n";
    }
}

std::vector<std::array<double, 3>> Positions::toStdVector()
{
    std::vector<std::array<double, 3>> positions(size());
    for (int i = 0; i < size(); i++)
    {
        positions[i] = {Xd[i], Yd[i], Zd[i]};
    }
    return positions;
}


void Positions::scaleUp()
{
    double scale = (1.0 / boxSize) * DEBYE_MAX_Rd;
    for (int i = 0; i < size(); i++)
    {
        X[i] = static_cast<int64>(std::round(Xd[i] * scale));
        Y[i] = static_cast<int64>(std::round(Yd[i] * scale));
        Z[i] = static_cast<int64>(std::round(Zd[i] * scale));
    }
}


void Positions::resize(int nAtoms)
{
    X.resize(nAtoms);
    Y.resize(nAtoms);
    Z.resize(nAtoms);
    Xd.resize(nAtoms);
    Yd.resize(nAtoms);
    Zd.resize(nAtoms);
    selectionIds.resize(nAtoms);
    chemicalSymbols.resize(nAtoms);
}

void Positions::sliceCellList(int start, int stop)
{
    // if the sliceCellList is the same as the original, return a copy of the original
    if (start == 0 and stop == size())
    {
        return;
    }

    // deal with edge cases
    {
        if (start < 0 || start >= size())
        {
            throw std::runtime_error("Start index is out of bounds.");
        }

        if (stop < 0 || stop > size())
        {
            throw std::runtime_error("Stop index is out of bounds.");
        }

        if (start >= stop)
        {
            throw std::runtime_error("Start index is greater than or equal to stop index.");
        }
    }

    if (cellHeads.empty()) // no cell list, no need to sliceCellList it
    {
        return;
    }


    // find the starting and end indices of the cell heads
    size_t startingCell = getCellIndex(start);
    size_t endCell = getCellIndex(stop-1); // stop is exclusive

    // fill the irrelevant cell heads with -1, atomCounts with 0
    for (size_t i = 0; i < startingCell; i++)
    {
        cellHeads[i] = -1;
        atomCounts[i] = 0;
    }

    for (size_t i = endCell+1; i < cellHeads.size(); i++)
    {
        cellHeads[i] = -1;
        atomCounts[i] = 0;
    }

    // set the starting cell head to start and atom count
    cellHeads[startingCell] = start;
    atomCounts[startingCell] = cellHeads[startingCell+1] - start;

    // update the last cell head to stop by updating the atom count of the last cell
    atomCounts[endCell] = stop - cellHeads[endCell];

}

size_t Positions::getCellIndex(size_t particleIdx) const {
    // find the first cell where the starting index is strictly greater than the particle index
    auto it = std::upper_bound(cellHeads.begin(), cellHeads.end(), particleIdx);

    if (it == cellHeads.begin()) // the particle is in the first cell
    {
        return 0;
    }

    // the cell before the one found is the cell that contains the particle
    return std::distance(cellHeads.begin(), it) - 1;
}