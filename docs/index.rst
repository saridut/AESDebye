.. AES-DEBYE documentation master file

Welcome to AES-DEBYE's documentation!
=====================================

AES-DEBYE is a high-performance software package designed for calculating the Debye scattering equation (DSE) and Pair Distribution Function (PDF) from atomic configurations. It is optimized for accuracy, speed, and scalability, making it suitable for large-scale simulations in materials science, chemistry, and physics.

Key Features & Methodology
--------------------------

AES-Debye presents an accuracy-preserving DSE framework with the following highlights:

* **Accuracy Preserving:** Aggregates pair distances into a pair distribution function (PDF) using corrected bin centers to suppress discretization artifacts.
* **Numerically Robust:** Uses robust accumulation algorithms to suppress floating-point summation errors.
* **Data Locality Aware:** Features a domain-decomposition-based design for highly local and predictable memory access, avoiding random-access cache degradation.
* **Hybrid Parallelization:** Harnesses OpenMP, MPI, and CUDA (CPU and GPU support) to scale to massive configurations (e.g., strong scalability demonstrated up to 90 million atoms).

If you use this code, please cite the following paper:

.. code-block:: bibtex

   @article{Panchi2026AES-Debye,
       author    = {Panchi, Navid and Kuckuk, Sebastian and Wittmann, Markus and Engel, Michael and Leonardi, Alberto},
       title     = {AES-Debye: an Accurate, Efficient, and Scalable solver for the Debye scattering equation},
       journal   = {Submitted to Journal of Applied Crystallography (IUCrJ), Under review}
       year      = {2026},
   }

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
