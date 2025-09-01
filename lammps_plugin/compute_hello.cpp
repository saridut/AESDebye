#include "lammpsplugin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "command.h"
#include "version.h"

#include <cstring>
#include <fstream>
#include <iostream>

namespace LAMMPS_NS {
  
  // Hello class inherits from the Command base class
  class Hello : public Command {
   public:
    // Constructor
    Hello(class LAMMPS *lmp) : Command(lmp) {};
    // This is the function that gets called when the command is executed
    void command(int, char **);
  };
}

using namespace LAMMPS_NS;

// The main command function
void Hello::command(int argc, char **argv)
{
  // Check for the correct number of arguments.
  if (argc != 2) error->all(FLERR, "Illegal hello command");

  // The second argument is the output filename.
  char* filename = argv[1];

  // Open the output file.
  std::ofstream outfile;
  // This check ensures that only the main processor writes to the file.
  if (comm->me == 0) {
    outfile.open(filename);
    if (!outfile.is_open()) {
      error->all(FLERR, "Cannot open output file");
    }
  }

  // Get the total number of atoms from the 'atom' object.
  // This is a property of the Atom class and holds the number of atoms
  // local to each processor.
  int nlocal = atom->nlocal;
  
  // Access the atomic position and type data directly from the 'atom' object.
  // 'x' is a pointer to a 2D array of doubles for positions (x, y, z).
  double** x = atom->x;
  // 'type' is a pointer to an array of integers for atom types.
  int* type = atom->type;

  // Loop over all local atoms and perform the operation.
  // The loop iterates from 0 up to the number of local atoms.
  for (int i = 0; i < nlocal; ++i) {
    // Perform a simple operation, e.g., print to console and file.
    // For a real application, you would do a more complex calculation here.
    if (comm->me == 0) {
      outfile << atom->tag[i] << " " << type[i] << " "
              << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
      
      // Example of an operation: print atoms with a specific type to console
      if (type[i] == 1) {
          utils::logmesg(lmp, fmt::format("Found atom type 1 with ID {} at {}, {}, {}", atom->tag[i], x[i][0], x[i][1], x[i][2]));
      }
    }
  }

  // Close the output file.
  if (comm->me == 0) {
    outfile.close();
    utils::logmesg(lmp, fmt::format("Successfully wrote atomic data to {}", filename));
  }
}

// Creator function for the plugin
void *myopcreator(LAMMPS *lmp)
{
  return new Hello(lmp);
}

// Plugin initialization function
extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style   = "command";
  plugin.name    = "hello";
  plugin.info    = "Operations on atomic data";
  plugin.author  = "Your Name";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &myopcreator;
  plugin.handle  = handle;
  (*register_plugin)(&plugin,lmp);
}