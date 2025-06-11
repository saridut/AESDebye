#pragma once

#include <DataTypes.hpp>
#include <Positions.hpp>
#include <iomanip>
#include <random>
#include <string>

/**
 * @brief Generates a random seed for the random number generator.
 *
 * This function generates a random seed for the random number generator.
 *
 * @return A random seed.
 */
unsigned generateRandomSeed();

/**
 * @brief Reads atomic positions from a file.
 *
 * This function reads atomic positions from an XYZ file.
 *
 * @param filename The name of the file to read.
 * @return The atomic positions.
 * @see Positions
 */
Positions readXYZ(const std::string &filename, std::string delimiter = " ",
                  uint skipLines = 2,
                  const std::string &typeMapping = "0:None");
/**
 * @brief Generates a random set of atomic positions.
 *
 * This function generates a random set of atomic positions.
 *
 * @param lattice The lattice constant.
 * @param nRepeats The number of repeats in each direction.
 * @param noise The std of the noise to add to the positions. This is in the
 * same units as the lattice constant.
 * @param seed The seed for the random number generator.
 * @param element The element to use.
 * @return The atomic positions.
 * @see Positions
 */
Positions generateTestData(double lattice, int nRepeats, double noise, int seed,
                           std::string element = "Pt");

Positions makePeriodic(const Positions &positions, int nx, int ny, int nz,
                       double offset_x, double offset_y, double offset_z);