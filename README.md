# AES-DEBYE: Accurate Efficient and Scalable Debye Scattering Calculation

AES-DEBYE is a high-performance software package designed for calculating Debye scattering patterns from atomic configurations. It is optimized for accuracy, speed and scalability, making it suitable for large-scale simulations in materials science, chemistry, and physics. If you end up using this code please cite the following paper:

```bibtex

```

### Installation

Current installation approach requires you to clone this repository and install it using `pip`. We recommend using a dedicated conda/mamba environment for this purpose.

```bash
git clone
cd aesdebye
pip install . -v
```

CMake will detect the availibility of MPI and CUDA automatically and enable those features. MPI is detected using `find_library` command from CMake, so please make sure your mpi library is accessible.

### Getting Started

#### Python Interface

Here is a simple example of how to use AES-DEBYE to calculate the Debye scattering pattern from a set of atomic positions.

```python
import aesdebye as debye

# Define atomic positions and types
elements = ['Pt', 'Pt', 'Pt', 'Pt']
coordinates = [
    [0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0],
    [2.0, 2.0, 2.0],
    [3.0, 3.0, 3.0]
]

# Create positions object
positions = debye.Positions(chemicalSymbols=elements,
                             coordinates=coordinates)

# Initialize Debye calculator - for more options see documentation
# print documentation using print(help(debye.DebyeCalculator))
calculator = debye.DebyeCalculator(nThreads=10,
                                   nCells = 15,
                                   useMPI=False,
                                   useGPU=False,
                                   verbose=True)

# compute and select the Pt-Pt profile
# This returns a dictionary with all the partials and a summed up one
# called "Full" - but only when we have more than one element
pdf, profile = calculator.calculateProfile(positions,
                                       start=0.0,
                                        stop=10.0,
                                        nsteps=1000)["Pt-Pt"]

# Save the results
profile.toCSV("profile_Pt.csv")

# plot
import matplotlib.pyplot as plt

plt.semilogy(profile.q, profile.intensity)
plt.xlabel("q (1/A)")
plt.ylabel("I(q)")
plt.show()
```

#### CLI interface

We provide a cli version of the calculator that can be used to calculate Debye scattering patterns directly from the command line.

Example usage:
```bash
aesdebye -f input_file -o output_dir -nc 15
```

For help use:

```bash
aesdebye -h

usage: aesdebye [-h] [-f INPUTFILENAME] [-tm TYPEMAPPING] [-nr NREPEATS] [-s STDDEV] [-sp] [-b BINSRESOLUTION] [-sb] [-nc NCELLS] [-nt NTHREADS] [-mpi] [-gpu] [-nlh] [-pc] [-fg] [-gcl] [-st START]
                [-e END] [-stps STEPS] [-wl WAVELENGTH] [-tt] [-o OUTPUTDIR] [-p] [-nv]

Debye implementation CLI

options:
  -h, --help            show this help message and exit

Input configuration:
  -f INPUTFILENAME, --inputFilename INPUTFILENAME
                        Input filename, any format supported by the ase.io.read function

Input configuration:
  -nr NREPEATS, --nRepeats NREPEATS
                        Number of repeats in the lattice
  -s STDDEV, --stdDev STDDEV
                        Std dev of noise in the lattice
  -sp, --shufflePositions
                        Shuffle the positions

Computation parameters:
  -b BINSRESOLUTION, --binsResolution BINSRESOLUTION
                        Resolution of the bins
  -sb, --smallBins      Use small bins for the PDF
  -nc NCELLS, --nCells NCELLS
                        Number of cells in the lattice
  -nt NTHREADS, --nThreads NTHREADS
                        Number of threads to use
  -mpi, --useMPI        Use MPI for parallelization
  -gpu, --useGPU        Use GPU for parallelization
  -nlh, --dontUseLocalHistogram
                        Dont use local histogram for noisy calcs
  -pc, --pseudoCoal     Use pseudo coal
  -fg, --fillGPU        Fill the GPU with threads
  -gcl, --useGPUCellList
                        Use GPU cell list

Range and physics settings:
  -st START, --start START
                        Start of the q/theta range
  -e END, --end END     End of the q/theta range
  -stps STEPS, --steps STEPS
                        Number of steps in the q/theta range
  -wl WAVELENGTH, --wavelength WAVELENGTH
                        Wavelength of the X-ray
  -tt, --twoThetaSpace  Use two theta instead of q

Output options:
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        Output directory to save the data
  -p, --plot            Plot the results
  -nv, --nonVerbose     Disable verbose output

```

### MPI usage

To run the CLI with MPI support, use the `mpirun` or `mpiexec` command. For example:
```bash
mpirun -n 4 aesdebye -f input_file -o output_dir -nc 15 -mpi
```
This will run the calculation using 4 MPI processes.

### GPU usage
To run the CLI with GPU support, simply add the `-gpu` flag to the command. For example:
```bash
aesdebye -f input_file -o output_dir -nc 15 -gpu
```

### Multi-GPU usage
You can use multiple GPUs by using one MPI process per GPU. For example, to use 2 GPUs:
```bash
mpirun -n 2 aesdebye -f input_file -o output_dir -nc 15 -mpi -gpu
```

### MPI and GPU with Python interface

The Python interface supports all of the options from the CLI. To use MPI and GPU support, simply set the `useMPI` and `useGPU` flags when initializing the `DebyeCalculator` class.
```python
calculator = debye.DebyeCalculator(nThreads=10,
                                   nCells = 15,
                                   useMPI=True,
                                   useGPU=True,
                                   verbose=True)
```

Let's say you save the above code in a file called `calculate.py`. You can then run it with MPI using the following command:
```bash
mpirun -n 2 python calculate.py
```
This will run the calculation using 2 MPI processes, each utilizing a GPU.

