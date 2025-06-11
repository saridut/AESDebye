#include <utils.hpp>
#include <map>


unsigned int generateRandomSeed()
{
    std::random_device rd;
    return rd();
}


auto setupTypeMap(std::string typeMapping)
{   if (typeMapping.empty())
    {
        return std::map<std::string, std::string>();
    }
    std::map<std::string, std::string> typeMap;
    auto types = helpers::stringSplit(typeMapping, ",");
    for (const auto& type : types)
    {
        auto pair = helpers::stringSplit(type, ":");
        typeMap[pair[0]] = pair[1];
    }
    return typeMap;
}


bool isDelimiter(char ch) {
    return std::isspace(static_cast<unsigned char>(ch));
}

// Function to find and return the first delimiter substring in the string
std::string findFirstDelimiterString(const std::string& str) {
    std::string delimiter;
    bool foundDelimiter = false;
    
    for (size_t i = 0; i < str.length(); ++i) {
        if (isDelimiter(str[i])) {
            foundDelimiter = true;
            delimiter += str[i];
        } else if (foundDelimiter) {
            break;  // Stop when the first delimiter sequence ends
        }
    }
    
    return delimiter;
}



Positions readXYZ(const std::string& filename, std::string delimiter,
                  uint skipLines, const std::string& typeMapping)
{
    auto typeMap = setupTypeMap(typeMapping);
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open file");
    }

    std::string buffer;
    std::getline(file, buffer); // natoms
    long nAtoms = std::stol(buffer);
    for (uint i = 0; i < skipLines; i++) // skip 1 line (or possibly more/less)
    {
        std::getline(file, buffer);
        // std::cout << "Buffer: " << buffer << std::endl;
    }

    std::vector<std::array<double, 3>> positionsArray(nAtoms);
    std::vector<std::string> selectionIds(nAtoms);
    std::vector<std::string> chemicalSymbols(nAtoms);

    // get line and check the number of columns
    std::string line;
    std::getline(file, line);
    auto tokens = helpers::stringSplit(line, delimiter);
    
    // try to change the delimiter
    if (tokens.size() == 1) {
        std::cerr << "Warning: The delimiter is not correct. Trying to find the correct delimiter." << std::endl;
        delimiter = findFirstDelimiterString(line);
        if (delimiter.empty())
        {
            throw std::runtime_error("Could not find the delimiter in the file. check the file format or give proper delimiter.");
        }
        std::cout << "Delimiter: " << delimiter << std::endl;
        tokens = helpers::stringSplit(line, delimiter);
    }

    if (tokens.size() > 5)
    {
        std::cerr << "Line: " << line << std::endl;
        std::cerr << "nColumns: " << tokens.size() << std::endl;
        std::cerr << "Warning: The XYZ file does not have the expected number of columns. "
                << " Reading the file anyway."
                <<  " The first column will be considered the atom type, the next three columns will be considered the positions and the last column will be considered the atom id."
                  << std::endl;
    }
    else if (tokens.size() < 4)
    {
        throw std::runtime_error("The XYZ file does not have the expected number of columns. "
                "The file should have at least 4 columns with atom type, x, y, z and atom id.");
    }

    chemicalSymbols[0] = typeMapping.empty() ? tokens[0] : typeMap[tokens[0]];
    positionsArray[0][0] = std::stod(tokens[1]);
    positionsArray[0][1] = std::stod(tokens[2]);
    positionsArray[0][2] = std::stod(tokens[3]);
    selectionIds[0] = tokens[tokens.size() - 1];

    double boxMin = std::min({positionsArray[0][0], positionsArray[0][1], positionsArray[0][2]});
    double boxMax = std::max({positionsArray[0][0], positionsArray[0][1], positionsArray[0][2]});

    for (long i = 1; i < nAtoms; ++i)
    {
        std::getline(file, line);
        auto tokens = helpers::stringSplit(line, delimiter);
        chemicalSymbols[i] = typeMapping.empty() ? tokens[0] : typeMap[tokens[0]];
        positionsArray[i][0] = std::stod(tokens[1]);
        positionsArray[i][1] = std::stod(tokens[2]);
        positionsArray[i][2] = std::stod(tokens[3]);
        selectionIds[i] = tokens[tokens.size() - 1];

        boxMin = std::min({boxMin, positionsArray[i][0], positionsArray[i][1], positionsArray[i][2]});
        boxMax = std::max({boxMax, positionsArray[i][0], positionsArray[i][1], positionsArray[i][2]});
    }

    Positions positions(chemicalSymbols, positionsArray, selectionIds, boxMin, boxMax);
    auto uniqueElements = positions.getUniqueElements();
    if (uniqueElements.size() == 1)
    {
        positions.element = uniqueElements[0];
    }
    return positions;
}

Positions generateTestData(double lattice, int nRepeats, double noise, int seed, std::string element)
{
    std::mt19937 gen(seed);
    std::normal_distribution<double> dist(0.0, noise);
    int nAtoms = 4 * nRepeats * nRepeats * nRepeats;
    std::vector<std::array<double, 3>> positionsArray(nAtoms);
    std::vector<std::string> chemicalSymbols(nAtoms, element);
    std::vector<std::string> selectionIds(nAtoms, "0");
    double half_lattice = lattice / 2.0;
    size_t counter = 0;

    for (int x = 0; x < nRepeats; x++)
    {
        for (int y = 0; y < nRepeats; y++)
        {
            for (int z = 0; z < nRepeats; z++)
            {
                positionsArray[counter] = {x * lattice, y * lattice, z * lattice};
                counter++;
                positionsArray[counter] = {x * lattice + half_lattice, y * lattice + half_lattice, z * lattice};
                counter++;
                positionsArray[counter] = {x * lattice + half_lattice, y * lattice, z * lattice + half_lattice};
                counter++;
                positionsArray[counter] = {x * lattice, y * lattice + half_lattice, z * lattice + half_lattice};
                counter++;
            }
        }
    }

    double boxMin = 0;
    double boxMax = lattice * nRepeats;
    if (noise > 0)
    {
        for (auto & pos : positionsArray)
        {
            pos[0] += dist(gen);
            pos[1] += dist(gen);
            pos[2] += dist(gen);

            boxMin = std::min({boxMin, pos[0], pos[1], pos[2]});
            boxMax = std::max({boxMax, pos[0], pos[1], pos[2]});
        }
    }

    Positions positions(chemicalSymbols, positionsArray, selectionIds, boxMin, boxMax);
    positions.element = element;
    return positions;
}



Positions makePeriodic(const Positions &positions, int nx, int ny, int nz,
                      double offset_x, double offset_y, double offset_z)
{
    std::vector<std::array<double, 3>> newPositions;
    double boxSize = positions.getBoxSize();
    offset_x = offset_x + boxSize;
    offset_y = offset_y + boxSize;
    offset_z = offset_z + boxSize;
    newPositions.reserve(nx * ny * nz * positions.size());

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                for (size_t l = 0; l < positions.size(); l++)
                {
                    newPositions.push_back({positions.Xd[l] + i * offset_x,
                                            positions.Yd[l] + j * offset_y,
                                            positions.Zd[l] + k * offset_z});
                }
            }
        }
    }

    std::vector<std::string> chemicalSymbols(newPositions.size(), positions.chemicalSymbols[0]);
    return {chemicalSymbols, newPositions};
}