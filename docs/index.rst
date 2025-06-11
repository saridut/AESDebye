.. PyAESDebye documentation master file, created by
   sphinx-quickstart on Mon Jun 10 21:58:57 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AESDebye's documentation!
======================================


This a high performance library to calculate accurate (full) pair distribution functions (PDFs) from atomic structues to be used for PDF analysis or to be used as a basis for the calculation of powder diffraction profiles using Debye scattering equation. This is the implementation of the algorithm described in the paper "AES-Debye: Accurate, Efficient and Scalable Implementation of Debye Scattering Equation". Please cite this paper if you use this library in your work. The algorithm was developed earlier by Leonardi et. al. (ref), please cite this paper as well.

It can be used in the following ways:

1. From command line (ref)
2. As a Python library (ref)
3. As a C++ library (ref)

The library has been parallelized using OpenMP for single node parallelization and MPI for multi-node parallelization. It is also GPU accelerated using CUDA.

You can follow the instructions here(link) to install the library. 

We support XYZ, and LAMMPSTRJ file formats natively.

Usage examples:
===============


.. currentmodule:: PyAESDebye

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   notebooks/SimpleExample

.. automodule:: PyAESDebye
   :members:
   :undoc-members:
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
