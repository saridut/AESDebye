// Compile as:
// g++ -fPIC -shared -o debye_profile.so debye_profile.cpp -I/path/to/lammps/src
// -I/path/to/debye/include

#include "atom.h"
#include "comm.h"
#include "command.h"
#include "error.h"
#include "lammpsplugin.h"
#include "version.h"

#include <Calculator.hpp> 
#include <utils.hpp>      
#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace LAMMPS_NS {

class DebyeProfileCommand : public Command {
public:
  DebyeProfileCommand(class LAMMPS *lmp) : Command(lmp) {}
  void command(int, char **);
  Positions perparePositions(std::string typeMapping);
};

} // namespace LAMMPS_NS

using namespace LAMMPS_NS;

Positions DebyeProfileCommand::perparePositions(std::string typeMapping) {
  // --- Step 2: Build Positions using readXYZ-style logic ---
  int nlocal = atom->nlocal;
  double **x = atom->x;
  int *type = atom->type;
  int *tag = atom->tag;

  // Prepare arrays
  std::vector<std::array<double, 3>> positionsArray(nlocal);
  std::vector<std::string> selectionIds(nlocal);
  std::vector<std::string> chemicalSymbols(nlocal);

  // Fill arrays and track bounding box
  double boxMin = std::numeric_limits<double>::max();
  double boxMax = std::numeric_limits<double>::lowest();

  // Set up type mapping if provided
  // Set up type mapping if provided
  auto typeMap = setupTypeMap(typeMapping);
  bool hasTypeMap = !typeMapping.empty();

  for (int i = 0; i < nlocal; i++) {
    std::string element;
    if (hasTypeMap) {
      // If type map provided, use it
      std::string typeStr = std::to_string(type[i]);
      element = typeMap.count(typeStr) ? typeMap[typeStr] : "None";
    } else {
      // If no map, assign "None" explicitly
      element = "None";
    }

    chemicalSymbols[i] = element;
    positionsArray[i][0] = x[i][0];
    positionsArray[i][1] = x[i][1];
    positionsArray[i][2] = x[i][2];
    selectionIds[i] = "0";

    boxMin = std::min({boxMin, x[i][0], x[i][1], x[i][2]});
    boxMax = std::max({boxMax, x[i][0], x[i][1], x[i][2]});
  }

  Positions positions(chemicalSymbols, positionsArray, selectionIds, boxMin,
                      boxMax);
  auto uniqueElements = positions.getUniqueElements();
  if (uniqueElements.size() == 1) {
    positions.element = uniqueElements[0];
  }

  return positions;
}

void DebyeProfileCommand::command(int argc, char **argv) {
  // --- Step 1: Parameters with defaults ---
  int nThreads = -1;
  int nCells = 15;
  double start = 0.0, end = 20.0;
  int nSteps = 2000;
  double lambda = 0.4;
  bool twoTheta = false;
  bool verbose = false;
  std::string outputDir = "";
  std::string typeMapping = "0:None";

  // Total number of required arguments
  const int REQUIRED_ARGS = 10; // nThreads, nCells, start, end, nSteps, lambda,
                                // twoTheta, verbose, outputDir, typeMapping

  if (argc < REQUIRED_ARGS) {
    std::cerr << "Error: Missing arguments. "
              << "Expected " << REQUIRED_ARGS << " arguments but got " << argc
              << ".\n"
              << "Usage: debye_profile nThreads nCells start end nSteps lambda "
                 "twoTheta verbose outputDir typeMapping\n";
    return; // or error->all(FLERR, ...) in LAMMPS
  }

  // Parse arguments
  nThreads = std::atoi(argv[0]);
  std::cout << "nThreads = " << nThreads << std::endl;
  nCells = std::atoi(argv[1]);
  std::cout << "nCells = " << nCells << std::endl;
  start = std::atof(argv[2]);
  std::cout << "start = " << start << std::endl;
  end = std::atof(argv[3]);
  std::cout << "end = " << end << std::endl;
  nSteps = std::atoi(argv[4]);
  std::cout << "nSteps = " << nSteps << std::endl;
  lambda = std::atof(argv[5]);
  std::cout << "lambda = " << lambda << std::endl;
  twoTheta = (std::atoi(argv[6]) != 0);
  std::cout << "twoTheta = " << twoTheta << std::endl;
  verbose = (std::atoi(argv[7]) != 0);
  std::cout << "verbose = " << verbose << std::endl;
  outputDir = std::string(argv[8]);
  std::cout << "outputDir = " << outputDir << std::endl;
  typeMapping = std::string(argv[9]);
  std::cout << "typeMapping = " << typeMapping << std::endl;

  // --- Step 3: Initialize DebyeCalculator ---
  DebyeCalculator debye(nThreads, nCells, 1.0,
                        /*useMPI*/ false, /*useGPU*/ false,
                        /*useLocalHistogram*/ true,
                        /*smallBins*/ false, /*pseudoCoal*/ false,
                        /*fillGPU*/ false, verbose);
  debye.config.useGPUCellList = false;

  // --- Step 4: Prepare Positions and calculate profiles ---
  auto positions = perparePositions(typeMapping);
  auto results = debye.calculateProfile(positions, start, end, nSteps, twoTheta,
                                        lambda, "");

  // --- Step 5: Save results ---
  if (comm->me == 0) {
    for (auto &[name, pdfProfile] : results) {
      auto &pdf = pdfProfile.first;
      auto &profile = pdfProfile.second;
      std::string out_prefix = outputDir + name;
      pdf.toCSV(out_prefix + "_pdf.csv");
      profile.toCSV(out_prefix + "_profile.csv");
    }
    utils::logmesg(lmp, fmt::format("Debye profiles written with prefix {}\n", outputDir));
  }
}

// Creator function
void *myopcreator(LAMMPS *lmp) { return new DebyeProfileCommand(lmp); }

// Plugin initialization
extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc) {
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc)regfunc;

  plugin.version = ".1.0";
  plugin.style = "command";
  plugin.name = "debye_profile";
  plugin.info = "Compute Debye scattering profiles using AESDebye library";
  plugin.author = "Navid Panchi";
  plugin.creator.v1 = (lammpsplugin_factory1 *)&myopcreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
