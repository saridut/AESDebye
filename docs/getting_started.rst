Getting Started
===============

This guide provides examples on how to use AES-DEBYE via the Python interface and the Command Line Interface (CLI).

Python Interface
----------------

Here is a simple example of how to use AES-DEBYE to calculate the Debye scattering pattern from a set of atomic positions.

.. code-block:: python

   import aesdebye as debye
   import matplotlib.pyplot as plt

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

   # Initialize Debye calculator
   calculator = debye.DebyeCalculator(nThreads=10,
                                      nCells=15,
                                      useMPI=False,
                                      useGPU=False,
                                      verbose=True)

   # Compute and select the Pt-Pt profile
   # This returns a dictionary with all the partials and a summed up one
   # called "total" (or "Pt-Pt" in this multi-element example)
   results = calculator.calculateProfile(positions,
                                         start=0.0,
                                         stop=10.0,
                                         nsteps=1000)
   pdf, profile = results["Pt-Pt"]

   # Save the results
   profile.toCSV("profile_Pt.csv")

   # Plot the scattering pattern
   plt.semilogy(profile.q, profile.intensity)
   plt.xlabel("q (1/A)")
   plt.ylabel("I(q)")
   plt.show()

Command Line Interface (CLI)
----------------------------

AES-DEBYE provides a command-line tool `aesdebye` to perform calculations directly from coordinate files.

Basic Usage
~~~~~~~~~~~

.. code-block:: bash

   aesdebye -i input_file.xyz -o output_dir -nc 15

To see all available CLI options, run:

.. code-block:: bash

   aesdebye -h

MPI and GPU Parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run the CLI with MPI support, use `mpirun` or `mpiexec`:

.. code-block:: bash

   mpirun -n 4 aesdebye -i input_file.xyz -o output_dir -nc 15 -mpi

To run with GPU support:

.. code-block:: bash

   aesdebye -i input_file.xyz -o output_dir -nc 15 -gpu

To run with multi-GPU support using one MPI process per GPU (e.g., 2 GPUs):

.. code-block:: bash

   mpirun -n 2 aesdebye -i input_file.xyz -o output_dir -nc 15 -mpi -gpu
