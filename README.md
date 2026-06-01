# AES-DEBYE: Accurate Efficient and Scalable Debye Scattering Calculation

AES-DEBYE is a high-performance software package designed for calculating Debye scattering equation and Pair distrubtion function (PDF) from atomic configurations. It is optimized for accuracy, speed and scalability, making it suitable for large-scale simulations in materials science, chemistry, and physics. If you end up using this code please cite the following paper:

```bibtex
@article{Panchi2026AES-Debye,
    author    = {Panchi, Navid and Kuckuk, Sebastian and Wittmann, Markus and Engel, Michael and Leonardi, Alberto},
    title     = {AES-Debye: an Accurate, Efficient, and Scalable solver for the Debye scattering equation},
    journal   = {Submitted to Journal of Applied Crystallography (IUCrJ)}
    year      = {2026},
}
```

### Installation

Current installation approach requires you to clone this repository and install it using `pip`. We recommend using a dedicated conda/mamba environment for this purpose. For more information on creating environments using conda/mamba see [here](https://docs.conda.io/en/latest/miniconda.html).

```bash
# git clone --recursive git@gitlab.cs.fau.de:iq23adyz/debye.git
# cd debye
# pip install . -v
```

CMake will detect the availibility of MPI and CUDA automatically and enable those features. MPI is detected using `find_library` command from CMake, so please make sure your mpi library is accessible. 

For example, on a HPC system, you can load the cuda and mpi modules using the following command (or similar, based on your HPC provider):

```bash
module load mpi cuda
```

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

You can provide `None` for the chemical symbols to avoid multiplication with the atomic scattering factors.

#### CLI interface

We provide a cli version of the calculator that can be used to calculate Debye scattering patterns directly from the command line.

Example usage:
```bash
aesdebye -f input_file -o output_dir -nc 15
```

For help use:

```bash
aesdebye -h
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

### Building Documentation

To build the HTML documentation, install the documentation requirements and run Sphinx:

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

The generated HTML files will be available in `docs/_build/html/`.

### Building as a C++ Library (Without Pip)

If you wish to use AES-DEBYE as a C++ library in another project rather than installing the Python bindings, you can build and install it directly via CMake:

```bash
cmake -B build -S . -DBUILD_PYTHON_BINDINGS=OFF -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build -j
cmake --install build
```

When `BUILD_PYTHON_BINDINGS` is disabled, CMake skips compiling the Python bindings (`_core` target) and instead installs the shared C++ libraries (`libaesdebye.so` and `libaesdebyeGPU.so`), the headers, and the C++ executable (`aesdebye_cpp`) to standard folders under the specified installation path.

#### Example: LAMMPS Plugin Integration

The [lammps_plugin](file:./lammps_plugin) directory contains an example of how to integrate the C++ library with other projects (specifically as a LAMMPS command extension plugin).

To build a plugin or command linking to the library, you link your target against the compiled `aesdebye` library target:

```cmake
add_library(lammps_debye_plugin SHARED compute_debye.cpp)
target_link_libraries(lammps_debye_plugin PRIVATE LAMMPS::lammps)
target_link_libraries(lammps_debye_plugin PRIVATE aesdebye)
```


