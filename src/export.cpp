#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/Calculator.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {

  py::module version = m.def_submodule("version", "Version information");
#ifdef USE_GPU
  version.attr("gpu_enabled") = true;
#else
  version.attr("gpu_enabled") = false;
#endif

#ifdef USE_MPI
  version.attr("mpi_enabled") = true;
#else
  version.attr("mpi_enabled") = false;
#endif
  version.attr("version") = VERSION_INFO;
  version.attr("compiler") = CXX_COMPILER;

  // options to make the docstring more informative
  m.doc() = R"pbdoc(
        aesdebye - High-Performance Debye Scattering Calculations

        aesdebye is a high-performance library for calculating accurate pair distribution
        functions (PDFs) from atomic structures using the Debye scattering equation. This
        implementation is based on the AES-Debye algorithm for efficient and scalable
        computation of scattering profiles.

        Key Features:
        - Fast PDF calculations using optimized algorithms
        - GPU acceleration support (CUDA)
        - MPI parallelization for multi-node computing
        - Support for XYZ and LAMMPS trajectory files
        - Efficient memory management and data structures
        - Python bindings for easy integration

        The library provides classes for atomic positions (Positions), pair distribution
        functions (PDF), scattering profiles (Profile), and the main calculator
        (DebyeCalculator) along with parallel computing utilities (ParallelHelper).

        Example:
            >>> import aesdebye
            >>> calc = aesdebye.DebyeCalculator()
            >>> positions = aesdebye.readXYZ("structure.xyz")
            >>> pdf = calc.calculatePDF(positions)
            >>> profile = calc.calculateIntensity(pdf, 0, 10, 1000)
        )pbdoc";

  py::options options;
  options.disable_function_signatures();

  py::class_<Profile>(m, "Profile", py::module_local(), R"pbdoc(
        Profile - Scattering intensity profile data structure

        The Profile class represents a calculated scattering intensity profile containing
        q-values, two-theta angles, and corresponding intensity values. This class is
        typically returned by intensity calculation methods and provides utilities for
        data export and visualization.

        Note:
            Profile objects cannot be directly instantiated in Python. They are created
            as results of calculations performed by DebyeCalculator methods.

        Attributes:
            q (List[float]): Momentum transfer values (Å⁻¹)
            twoTheta (List[float]): Scattering angles in degrees
            intensity (List[float]): Calculated intensity values
            calculationTime (float): Time taken for the calculation in seconds
            testPassed (bool): Flag indicating if internal validation tests passed
        )pbdoc")
      .def("toCSV", &Profile::toCSV, py::arg("filename"), R"pbdoc(
        toCSV(filename: str) -> None

        Export the scattering profile to a CSV file.

        The CSV file will contain columns for q-values, two-theta angles, and
        intensity values, making it suitable for further analysis or plotting
        with external tools.

        Parameters
        ----------
        filename : str
            Path to the output CSV file. The file will be created or overwritten
            if it already exists.

        Example
        -------
        >>> profile = calc.calculateIntensity(pdf, 0, 10, 1000)
        >>> profile.toCSV("scattering_profile.csv")
        )pbdoc")

      // operators
      .def(py::self + py::self, R"pbdoc(
        Add two Profile objects element-wise.

        Returns
        -------
        Profile
            New Profile with summed intensity values at corresponding q-points.
        )pbdoc")
      .def(py::self - py::self, R"pbdoc(
        Subtract two Profile objects element-wise.

        Returns
        -------
        Profile
            New Profile with subtracted intensity values at corresponding q-points.
        )pbdoc")

      // __str__ method
      .def("__repr__", &Profile::toString)

      // props
      .def_readonly("q", &Profile::q,
                    "Momentum transfer values in Å⁻¹. These are the q-values "
                    "at which the intensity was calculated.")
      .def_readonly("twoTheta", &Profile::twoTheta,
                    "Scattering angles in degrees (2θ). Converted from "
                    "q-values using the specified wavelength.")
      .def_readonly("intensity", &Profile::intensity,
                    "Calculated scattering intensity values corresponding to "
                    "each q-value.")
      .def_readonly(
          "calculationTime", &Profile::calculationTime,
          "Computation time in seconds required to calculate this profile.")
      .def_readonly("testPassed", &Profile::testPassed,
                    "Boolean flag indicating whether internal validation tests "
                    "passed during calculation.");

  // Classes
  py::class_<ParallelHelper>(m, "ParallelHelper", py::module_local())
      .def(py::init<int, bool, bool>(), py::arg("nThreads") = -1,
           py::arg("useMPI") = false, py::arg("verbose") = true)
      .def("setNumThreads", &ParallelHelper::setNumThreads, py::arg("nThreads"),
           R"pbdoc(
                setNumThreads(nThreads: int) -> None

                Set the number of OpenMP threads.

                Parameters
                ----------
                nThreads : int
                    The number of OpenMP threads to use.
                )pbdoc")
      .def("wait", &ParallelHelper::wait, R"pbdoc(
                wait() -> None

                Wait for all MPI ranks to synchronize.
                )pbdoc")
      .def_readonly("world_size", &ParallelHelper::worldSize,
                    "Total size of the MPI world")
      .def_readonly("world_rank", &ParallelHelper::worldRank,
                    "Rank of the current process in the MPI world")
      .def_readonly("ompThreads", &ParallelHelper::ompThreads,
                    "Number of OpenMP threads")
      .def_readonly("totalThreads", &ParallelHelper::totalThreads,
                    "Total number of threads")
      .doc() = R"pbdoc(
                ParallelHelper(nThreads: int = -1, useMPI: bool = False, verbose: bool = True) -> None

                Construct a new ParallelHelper object.

                Parameters
                ----------
                nThreads : int, optional
                        The number of threads to use. The default is -1.
                useMPI : bool, optional
                        Flag to indicate if MPI should be used. The default is False.
                verbose : bool, optional
                        Flag to indicate if verbose output should be enabled. The default is True.
                )pbdoc";

  py::class_<Positions>(m, "Positions", py::module_local(), R"pbdoc(
        Positions - Atomic structure data container

        The Positions class represents atomic positions and associated metadata for a
        molecular or crystalline system. It stores atomic coordinates, chemical symbols,
        selection IDs, and simulation box information. This class is the primary input
        for PDF and intensity calculations.

        The class provides methods for data manipulation including filtering by element
        type or selection ID, subsampling, and file I/O operations. It also handles
        periodic boundary conditions and spatial partitioning for efficient calculations.

        Attributes:
            element (str): Default element type for all atoms
            selectionIds (List[str]): Unique identifiers for each atom
            chemicalSymbols (List[str]): Chemical element symbols for each atom

        Example:
            >>> positions = aesdebye.readXYZ("structure.xyz")
            >>> print(f"Number of atoms: {positions.size()}")
            >>> pt_atoms = positions.filterByElement("Pt")
        )pbdoc")
      .def(
          py::init<std::vector<std::string>, std::vector<std::array<double, 3>>,
                   std::vector<std::string>, double, double>(),
          py::arg("chemicalSymbols"), py::arg("coordinates"),
          py::arg("selectionIds") = std::vector<std::string>(),
          py::arg("boxMin") = 1e100, py::arg("boxMax") = -1e100, R"pbdoc(
        Construct a new Positions object from atomic data.

        Parameters
        ----------
        chemicalSymbols : List[str]
            Chemical element symbols for each atom (e.g., ['C', 'O', 'H'])
        coordinates : List[Tuple[float, float, float]]
            Cartesian coordinates for each atom in Angstroms
        selectionIds : List[str], optional
            Unique identifiers for each atom. If empty, indices will be used.
        boxMin : float, optional
            Minimum coordinate value for automatic box size detection
        boxMax : float, optional
            Maximum coordinate value for automatic box size detection

        Example
        -------
        >>> symbols = ['C', 'C', 'O']
        >>> coords = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (0.75, 1.3, 0.0)]
        >>> positions = aesdebye.Positions(symbols, coords)
        )pbdoc")
      .def_readwrite(
          "element", &Positions::element,
          "Default element type for atoms when not specified individually")
      .def_readwrite(
          "selectionIds", &Positions::selectionIds,
          "List of unique atom identifiers for selection and filtering")
      .def_readwrite(
          "chemicalSymbols", &Positions::chemicalSymbols,
          "List of chemical element symbols corresponding to each atom")
      // static method for fromIndexList
      .def_static("fromIndexList", &Positions::fromIndexList,
                  py::arg("positions"), py::arg("indexList"))
      .def("filterBySelectionId", &Positions::filterBySelectionId,
           py::arg("atomId"), R"pbdoc(
                filterBySelectionId(atomId: int) -> Positions

                Filter the positions by atom id.

                Parameters
                ----------
                atomId : int
                    The atom id to filter by.

                Returns
                -------
                Positions
                    The filtered positions object.
                )pbdoc")
      .def("filterByElement", &Positions::filterByElement,
           py::arg("elementName"), R"pbdoc(
                filterByElement(atomType: int) -> Positions

                Filter the positions by atom type.

                Parameters
                ----------
                atomType : int
                    The atom type to filter by.

                Returns
                -------
                Positions
                    The filtered positions object.
                )pbdoc")
      .def("getUniqueSelectionIds", &Positions::getUniqueSelectionIds, R"pbdoc(
                getUniqueSelectionIds() -> List[int]

                Get the unique atom ids in the positions object.

                Returns
                -------
                List[int]
                    List of unique atom ids.
                )pbdoc")
      .def("getUniqueElements", &Positions::getUniqueElements, R"pbdoc(
                getUniqueElements() -> List[int]

                Get the unique atom types in the positions object.

                Returns
                -------
                List[int]
                    List of unique atom types.
                )pbdoc")
      .def("size", &Positions::size, R"pbdoc(
                size() -> int

                Get the number of atoms in the positions object.

                Returns
                -------
                int
                    The number of atoms in the positions object.
                )pbdoc")
      .def("getBoxSize", &Positions::getBoxSize, R"pbdoc(
                getBoxSize() -> float

                Get the box size of the positions object.

                Returns
                -------
                float
                    The box size of the positions object.
                )pbdoc")
      .def("setBoxSize", &Positions::setBoxSize, py::arg("boxSize"), R"pbdoc(
                setBoxSize(boxSize: float) -> None

                Set the box size of the positions object.

                Parameters
                ----------
                boxSize : float
                    The box size to set.
                )pbdoc")
      .def("subsample", &Positions::subsample, py::arg("subSampleRatio"),
           py::arg("seed") = 1, R"pbdoc(
                subsample(subSampleRatio: float, seed: int = 1) -> Positions

                Subsample the positions object to a given number of atoms.

                Parameters
                ----------
                subSampleRatio : float
                    The ratio of atoms to keep.
                seed : int, optional
                        The seed for the random number generator. The default is 1.

                Returns
                -------
                Positions
                    The subsampled positions object.
                )pbdoc")
      .def("sliceCellList", &Positions::sliceCellList, py::arg("start"),
           py::arg("end"), R"pbdoc(
                slice(start: int, end: int) -> Positions

                Slice the positions object.

                Parameters
                ----------
                start : int
                    The start index of the slice.
                end : int
                    The end index of the sliceCellList.

                Returns
                -------
                Positions
                    The sliced positions object.
                )pbdoc")
      .def("toList", &Positions::toStdVector, R"pbdoc(
                toList() -> List[Tuple[float, float, float]]

                Create a list of atomic positions.

                Returns
                -------
                List[Tuple[float, float, float]]
                    List of positions.
                )pbdoc")
      .def("toXYZ", &Positions::toXYZ, py::arg("filename"), R"pbdoc(
                toXYZ(filename: str) -> None

                Write the positions to an XYZ file.

                Parameters
                ----------
                filename : str
                    The path to the output file.
                )pbdoc")
      .doc() = R"pbdoc(
                Positions(positions: List[Tuple[float, float, float]], boxMin: Tuple[float, float, float], boxMax: Tuple[float, float, float]) -> None

                Construct a new Positions object.

                Parameters
                ----------
                positions : List[Tuple[float, float, float]]
                    List of positions
                boxMin : Tuple[float, float, float]
                    Minimum box size
                boxMax : Tuple[float, float, float]
                    Maximum box size
                )pbdoc";

  py::class_<CellList>(m, "CellList", py::module_local())
      .def(py::init<uint>(), py::arg("nCells"))
      .def("createCellList", &CellList::createCellList, py::arg("positions"),
           R"pbdoc(
                createCellList(positions: Positions) -> None

                Create the cell list for the given positions.

                Parameters
                ----------
                positions : Positions
                    The positions object.
                )pbdoc");

  py::class_<PDF>(m, "PDF", py::module_local(), R"pbdoc(
        PDF - Pair Distribution Function data structure

        The PDF class represents a calculated pair distribution function containing
        histogram data with bin centers, counts, and corrected centers. This class
        stores the results of PDF calculations and provides methods for data export,
        visualization, and further analysis.

        The PDF contains both corrected and uncorrected bin centers, where corrected
        centers account for the finite bin size effects in the histogram. The class
        supports arithmetic operations for combining multiple PDFs and provides
        efficient storage using internal histogram structures.

        Attributes:
            centers (List[float]): Corrected bin centers accounting for finite bin effects
            uncorrectedCenters (List[float]): Original bin centers without corrections
            counts (List[float]): Histogram counts for each bin

        Note:
            PDF objects are typically created as results of DebyeCalculator.calculatePDF()
            methods rather than being instantiated directly.

        Example:
            >>> calc = aesdebye.DebyeCalculator()
            >>> pdf = calc.calculatePDF(positions)
            >>> pdf.toCSV("output.csv")
            >>> print(f"PDF has {len(pdf.centers)} bins")
        )pbdoc")
      .def(py::init<double, bool>(), py::arg("boxSize") = 1,
           py::arg("samePositions") = true, R"pbdoc(
        Construct a new PDF object.

        Parameters
        ----------
        boxSize : float, optional
            Size of the simulation box in Angstroms. Used for proper normalization
            and bin center calculations. Default is 1.0.
        samePositions : bool, optional
            Flag indicating whether this PDF represents correlations within the same
            set of positions (self-correlation) or between different position sets
            (cross-correlation). Default is True.

        Note:
            Direct instantiation is rarely needed. PDFs are typically created by
            DebyeCalculator methods.
        )pbdoc")
      .def(py::self += py::self, R"pbdoc(
        Add another PDF to this PDF in-place.

        The two PDFs must have compatible box sizes and bin structures.
        This operation combines the histogram counts and updates the
        internal data structures accordingly.

        Parameters
        ----------
        other : PDF
            Another PDF object to add to this one

        Returns
        -------
        PDF
            Reference to this PDF object after addition
        )pbdoc")
      .def("__repr__", &PDF::toString)
      .def("save", &PDF::save, py::arg("filename"), R"pbdoc(
        save(filename: str) -> None

        Save the PDF to a binary file for efficient storage and loading.

        The binary format preserves all internal data structures and metadata,
        allowing for exact reconstruction of the PDF object. This is the
        recommended format for intermediate storage during calculations.

        Parameters
        ----------
        filename : str
            Path to the output binary file. The file will be created or
            overwritten if it already exists.

        See Also
        --------
        load : Load a PDF from a binary file
        toCSV : Export to human-readable CSV format

        Example
        -------
        >>> pdf.save("my_pdf.bin")
        )pbdoc")
      .def("load", &PDF::load, py::arg("filename"), R"pbdoc(
        load(filename: str) -> float

        Load the PDF from a binary file created with save().

        This method reconstructs the complete PDF object from the binary
        file, including all internal histogram data and metadata.

        Parameters
        ----------
        filename : str
            Path to the input binary file created with save()

        Returns
        -------
        float
            The box size of the loaded PDF

        Example
        -------
        >>> pdf = aesdebye.PDF()
        >>> box_size = pdf.load("my_pdf.bin")
        )pbdoc")
      .def("toCSV", &PDF::toCSV, py::arg("filename"),
           py::arg("complete") = false, R"pbdoc(
        toCSV(filename: str, complete: bool = False) -> None

        Export the PDF to a CSV file for analysis and visualization.

        The CSV file contains columns for bin centers, corrected centers,
        and counts, making it suitable for plotting and further analysis
        with external tools like Excel, MATLAB, or Python plotting libraries.

        Parameters
        ----------
        filename : str
            Path to the output CSV file. The file will be created or
            overwritten if it already exists.
        complete : bool, optional
            If False (default), only non-empty bins are exported to reduce
            file size. If True, all bins including empty ones are exported.
            Setting to True can result in very large files.

        Example
        -------
        >>> pdf.toCSV("pdf_data.csv")
        >>> pdf.toCSV("complete_pdf.csv", complete=True)
        )pbdoc")
      .def("readCSV", &PDF::readCSV, py::arg("filename"), R"pbdoc(
        readCSV(filename: str) -> None

        Read PDF data from a CSV file.

        This method loads PDF data from a CSV file with the expected format
        (bin centers, corrected centers, counts). The CSV should match the
        format produced by toCSV().

        Parameters
        ----------
        filename : str
            Path to the input CSV file containing PDF data

        Note:
            The CSV file must have the correct format with appropriate columns
            for bin centers, corrected centers, and counts.
        )pbdoc")
      // properties
      .def_property_readonly(
          "uncorrectedCenters", &PDF::getUncorrectedCenters,
          "List of uncorrected bin centers in Angstroms. These are the nominal "
          "bin positions without finite-size corrections.")
      .def_property_readonly(
          "centers", &PDF::getCenters,
          "List of corrected bin centers in Angstroms. These account for "
          "finite bin size effects and should be used for accurate analysis.")
      .def_property_readonly(
          "counts", &PDF::getCounts,
          "List of histogram counts for each bin. These represent the number "
          "of atom pairs found at each distance.")
      .def_readonly("calculationTime", &PDF::calculationTime,
                    "Time taken for the PDF calculation in seconds.")
      .def_readonly("testPassed", &PDF::testPassed,
                    "Boolean flag indicating whether internal "
                    "validation tests passed during PDF calculation.");

  py::class_<DebyeCalculator>(m, "DebyeCalculator", py::module_local(), R"pbdoc(
        DebyeCalculator - Main computational engine for Debye scattering calculations

        The DebyeCalculator class is the primary interface for performing pair distribution
        function (PDF) and scattering intensity calculations using the AES-Debye algorithm.
        It provides high-performance computation with support for parallel processing via
        OpenMP and MPI, as well as GPU acceleration through CUDA.

        The calculator handles the complete workflow from atomic positions to scattering
        profiles, including spatial partitioning via cell lists, histogram generation,
        and intensity calculations with atomic scattering factors.

        Key Features:
        - Efficient PDF calculations using optimized algorithms
        - Cell list spatial partitioning for O(N) scaling
        - OpenMP parallelization for multi-core systems
        - MPI support for distributed computing
        - GPU acceleration (CUDA) when available
        - Flexible bin resolution and histogram management

        Example:
            >>> calc = aesdebye.DebyeCalculator(nThreads=8, useGPU=True)
            >>> pdf = calc.calculatePDF(positions)
            >>> profile = calc.calculateIntensity(pdf, 0, 10, 1000)
        )pbdoc")
      .def(py::init<int, int, double, bool, bool, bool, bool, bool, bool,
                    bool>(),
           py::arg("nThreads") = -1, py::arg("nCells") = 15,
           py::arg("binsResolution") = 1.0, py::arg("useMPI") = false,
           py::arg("useGPU") = false, py::arg("useLocalHist") = true,
           py::arg("smallBin") = false, py::arg("pseudoCoal") = true,
           py::arg("fillGPU") = false, py::arg("verbose") = true, R"pbdoc(
        Construct a new DebyeCalculator with specified computational parameters.

        Parameters
        ----------
        nThreads : int, optional
            Number of OpenMP threads to use. Default is -1 (use all available cores).
            Set to 1 for single-threaded execution.
        nCells : int, optional
            Number of cells per dimension for spatial partitioning. Default is 15.
            Higher values improve performance for large systems but use more memory.
            Set to 1 to disable cell lists (only for very small systems).
        binsResolution : float, optional
            Histogram bin resolution factor in range (0, 1]. Default is 1.0.
            Lower values increase bin resolution but require more memory.
        useMPI : bool, optional
            Enable MPI parallelization. Default is False. When True, the program
            must be launched with an MPI runner (e.g., mpirun).
        useGPU : bool, optional
            Enable GPU acceleration via CUDA. Default is False. Requires GPU-enabled
            build and compatible hardware.
        useLocalHist : bool, optional
            Use local histograms for thread-level parallelization. Default is True.
            Set to False for crystalline structures with high symmetry.
        smallBin : bool, optional
            Use smaller bin sizes for higher resolution. Default is False.
        pseudoCoal : bool, optional
            Enable pseudo-coalescence optimization. Default is True.
        fillGPU : bool, optional
            Fill GPU memory completely for maximum performance. Default is False.
        verbose : bool, optional
            Enable verbose output for debugging and monitoring. Default is True.

        Raises
        ------
        RuntimeError
            If GPU support is requested but not available in the build.
        ValueError
            If binsResolution is not in the valid range (0, 1].

        Example
        -------
        >>> # Basic calculator for small systems
        >>> calc = aesdebye.DebyeCalculator()
        >>>
        >>> # High-performance setup for large systems
        >>> calc = aesdebye.DebyeCalculator(
        ...     nThreads=16,
        ...     nCells=20,
        ...     useGPU=True,
        ...     verbose=False
        ... )
        )pbdoc")
      .def("calculatePDF",
           py::overload_cast<Positions &, Positions &>(
               &DebyeCalculator::calculatePDF),
           py::arg("positionsI"), py::arg("positionsJ").none(false), R"pbdoc(
        calculatePDF(positionsI: Positions, positionsJ: Positions) -> PDF

        Calculate the pair distribution function between two sets of atomic positions.

        This method computes the cross-correlation PDF between two different position
        sets, which is useful for analyzing correlations between different atom types
        or different regions of a system.

        Parameters
        ----------
        positionsI : Positions
            First set of atomic positions. This defines the reference atoms for
            the pair correlation calculation.
        positionsJ : Positions
            Second set of atomic positions. This defines the target atoms for
            the pair correlation calculation.

        Returns
        -------
        PDF
            Calculated pair distribution function containing histogram data with
            bin centers, counts, and corrected centers accounting for finite
            bin size effects.

        Example
        -------
        >>> calc = aesdebye.DebyeCalculator()
        >>> pt_atoms = positions.filterByElement("Pt")
        >>> au_atoms = positions.filterByElement("Au")
        >>> pdf = calc.calculatePDF(pt_atoms, au_atoms)  # Pt-Au correlations
        )pbdoc")

      .def("calculatePDF",
           py::overload_cast<Positions &>(&DebyeCalculator::calculatePDF),
           py::arg("positionsI"), R"pbdoc(
        calculatePDF(positionsI: Positions) -> PDF

        Calculate the pair distribution function for a single set of atomic positions.

        This method computes the self-correlation PDF within a single position set,
        which represents the standard radial distribution function showing how
        atomic density varies as a function of distance from any given atom.

        Parameters
        ----------
        positionsI : Positions
            Atomic positions for which to calculate the PDF. The method will
            compute all pairwise distances within this set.

        Returns
        -------
        PDF
            Calculated pair distribution function containing histogram data with
            bin centers, counts, and corrected centers. The samePositions flag
            will be set to True for proper normalization.

        Example
        -------
        >>> calc = aesdebye.DebyeCalculator()
        >>> positions = aesdebye.readXYZ("structure.xyz")
        >>> pdf = calc.calculatePDF(positions)
        >>> pdf.toCSV("output_pdf.csv")
        )pbdoc")

      .def("calculateIntensity",
           py::overload_cast<
               std::vector<double> const &, std::vector<double> const &,
               std::vector<double> const &, std::string, std::string>(
               &DebyeCalculator::calculateIntensity),
           py::arg("centers"), py::arg("counts"), py::arg("qVector"),
           py::arg("elementI") = "", py::arg("elementJ") = "", R"pbdoc(
        calculateIntensity(centers: List[float], counts: List[float], qVector: List[float], elementI: str = "", elementJ: str = "") -> List[float]

        Calculate scattering intensity from PDF data at specified q-values.

        This method computes the scattering intensity using the Debye scattering equation
        from pre-calculated PDF histogram data. It applies atomic scattering factors
        and performs the Fourier transform to convert from real space (PDF) to
        reciprocal space (intensity).

        Parameters
        ----------
        centers : List[float]
            Corrected bin centers from PDF calculation in Angstroms. Use PDF.centers
            for best accuracy as these account for finite bin size effects.
        counts : List[float]
            Histogram counts from PDF calculation. Use PDF.counts to get the
            pair correlation data.
        qVector : List[float]
            Momentum transfer values (q) in Å⁻¹ at which to calculate intensity.
            Typical range is 0.5 to 25 Å⁻¹ for X-ray scattering.
        elementI : str, optional
            Chemical symbol for first element type (e.g., "Pt", "Au"). If empty,
            uses generic scattering factors. Default is "".
        elementJ : str, optional
            Chemical symbol for second element type. For self-correlation PDFs,
            this should match elementI. Default is "".

        Returns
        -------
        List[float]
            Calculated scattering intensity values corresponding to each q-value.
            Units are arbitrary but proportional to scattered intensity.

        Example
        -------
        >>> calc = aesdebye.DebyeCalculator()
        >>> pdf = calc.calculatePDF(positions)
        >>> q_values = [i * 0.1 for i in range(1, 101)]  # 0.1 to 10.0 Å⁻¹
        >>> intensity = calc.calculateIntensity(pdf.centers, pdf.counts, q_values, "Pt")
        )pbdoc")

      .def("calculateIntensity",
           py::overload_cast<PDF &, double, double, int, bool, double>(
               &DebyeCalculator::calculateIntensity),
           py::arg("pdf"), py::arg("start"), py::arg("end"), py::arg("nSteps"),
           py::arg("twoThetaSpace") = false, py::arg("wavelength") = 0.4)

      .def("calculateProfile",
           py::overload_cast<Positions &, double, double, int, bool, double,
                             std::string>(&DebyeCalculator::calculateProfile),
           py::arg("positions"), py::arg("start") = 0, py::arg("end") = 10,
           py::arg("steps") = 1000, py::arg("twoTheta") = false,
           py::arg("wavelength") = 0.4, py::arg("filter") = "", R"pbdoc(
        calculateProfile(
            positions: Positions,
            start: float = 0,
            end: float = 10,
            steps: int = 1000,
            twoTheta: bool = False,
            wavelength: float = 0.4,
            filter: str = ""
        ) -> Tuple[PDF, Profile]

        Calculate the scattering intensity profile for a given set of positions.

        Parameters
        ----------
        positions : Positions
            The positions object containing atomic coordinates.
        start : float, optional
            The starting q (or angle) value. Default is 0.
        end : float, optional
            The ending q (or angle) value. Default is 10.
        steps : int, optional
            Number of q steps. Default is 1000.
        twoTheta : bool, optional
            If True, interpret the range as 2θ values instead of q. Default is False.
        wavelength : float, optional
            Wavelength of the incident beam in Å. Default is 0.4.
        filter : str, optional
            Optional filter for scattering contributions (e.g., element selection). Default is "".

        Returns
        -------
        Tuple[PDF, Profile]
            Tuple containing:
            - PDF object with pair distribution data,
            - Profile object with the scattering intensity profile.

        Example
        -------
        >>> calc = aesdebye.DebyeCalculator()
        >>> positions = aesdebye.readXYZ("structure.xyz")
        >>> results = calc.calculateProfile(positions, start=0, end=10, steps=1000, wavelength=0.4)
        >>> pdf, profile = results["Pt-Pt"]  # Accessing the PDF and Profile for Pt-Pt correlations
        >>> pdf_total, profile_total = results["total"]  # Accessing the total PDF and Profile
     )pbdoc")
      .def(
          "calculateProfile",
          [](DebyeCalculator &self,
             const std::vector<std::array<double, 3>> &coords,
             std::optional<std::vector<std::string>> elements, double start,
             double end, int steps, bool twoTheta, double wavelength,
             std::string filter) {
            std::vector<std::string> el_list;
            if (!elements || elements->empty()) {
              el_list.assign(coords.size(), "None");
            } else {
              el_list = *elements;
            }

            // Initialize Positions and call existing method
            Positions pos(el_list, coords);
            return self.calculateProfile(pos, start, end, steps, twoTheta,
                                         wavelength, filter);
          },
          py::arg("positions"), py::arg("elements") = py::none(),
          py::arg("start") = 0, py::arg("end") = 10, py::arg("steps") = 1000,
          py::arg("twoTheta") = false, py::arg("wavelength") = 0.4,
          py::arg("filter") = "")

      .def_static("calculateASFProfile", &DebyeCalculator::calculateASFProfile,
                  py::arg("qVector"), py::arg("element"), R"pbdoc(
        calculateASFProfile(qVector: List[float], element: str) -> List[float]

        Calculate the atomic scattering factor (ASF) profile for a given element
        at specified q-values [values have to be Å⁻¹].

        This method evaluates the atomic scattering factor as a function of
        momentum transfer q using tabulated scattering factor coefficients.

        Parameters
        ----------
        qVector : List[float]
            Momentum transfer values (q) in Å⁻¹ at which to evaluate the ASF.
        element : str
            Chemical symbol of the element (e.g., "Pt", "Au").

        Returns
        -------
        List[float]
            Atomic scattering factor values corresponding to each q-value.

        Example
        -------
        >>> q_values = [i * 0.1 for i in range(1, 101)]  # 0.1 to 10.0 Å⁻¹
        >>> asf_profile_values = DebyeCalculator.calculateASFProfile(q_values, "Pt")
    )pbdoc")
      .def_readonly("parallelHelper", &DebyeCalculator::parallelHelper,
                    "ParallelHelper object created by the DebyeCalculator",
                    R"pbdoc(
        parallelHelper: ParallelHelper
        ParallelHelper object created by the DebyeCalculator for managing
        parallel execution settings and MPI communication.
        )pbdoc")
      .def_readonly("cellList", &DebyeCalculator::cellList, R"pbdoc(
        cellList: CellList
        CellList object used for spatial partitioning of atomic positions.
        )pbdoc");

  m.def("generateData", &generateTestData, py::arg("lattice"),
        py::arg("nRepeats"), py::arg("noise") = 0.0, py::arg("seed") = 1,
        py::arg("element") = "Pt",
        R"pbdoc(
        generateData(lattice: float, nRepeats: int, noise: float = 0.0, seed: int = 1, element: str = "Pt") -> Positions

        Generate synthetic FCC (face-centered cubic) crystal structures for testing and benchmarking.

        This utility function creates perfect FCC lattice structures with optional Gaussian
        noise for realistic atomic positions. It's useful for testing algorithms, benchmarking
        performance, and creating reference structures for validation.

        Parameters
        ----------
        lattice : float
            Lattice constant in Angstroms. Typical values are 3.5-4.0 Å for metals like Pt, Au.
        nRepeats : int
            Number of unit cell repetitions in each direction (x, y, z). Total atoms will be
            approximately 4 * nRepeats³ for FCC structure.
        noise : float, optional
            Standard deviation of Gaussian noise added to atomic positions in Angstroms.
            Default is 0.0 (perfect crystal). Typical values: 0.01-0.1 Å for thermal motion.
        seed : int, optional
            Random seed for reproducible noise generation. Default is 1.
        element : str, optional
            Chemical element symbol for all atoms. Default is "Pt" (platinum).

        Returns
        -------
        Positions
            Generated atomic positions with FCC structure and specified parameters.

        Example
        -------
        >>> # Perfect platinum crystal
        >>> positions = aesdebye.generateData(lattice=3.92, nRepeats=10)
        >>>
        >>> # Realistic structure with thermal motion
        >>> positions = aesdebye.generateData(lattice=3.92, nRepeats=5, noise=0.05, element="Au")
        >>> print(f"Generated {positions.size()} atoms")
        )pbdoc");

  m.def("readXYZ", &readXYZ, py::arg("filename"), py::arg("delimiter") = " ",
        py::arg("skipLines") = 2, py::arg("typeMapping") = "0:None", R"pbdoc(
        readXYZ(filename: str, delimiter: str = " ", skipLines: int = 2, typeMapping: str = "0:None") -> Positions

        Read atomic positions from an XYZ format file.

        This function parses standard XYZ files containing atomic coordinates and element
        information. It supports various XYZ formats including those with additional
        columns for atom types or selection IDs.

        Parameters
        ----------
        filename : str
            Path to the XYZ file to read. The file should contain atomic coordinates
            in standard XYZ format.
        delimiter : str, optional
            Column delimiter used in the file. Default is " " (space). Use "," for
            comma-separated files or "\t" for tab-separated files.
        skipLines : int, optional
            Number of header lines to skip at the beginning of the file. Default is 2
            (standard XYZ format with atom count and comment lines).
        typeMapping : str, optional
            Mapping string for atom type columns. Default is "0:None". Format is
            "column_index:element_symbol" for custom element assignment.

        Returns
        -------
        Positions
            Positions object containing the atomic coordinates, element symbols,
            and metadata from the XYZ file.

        Raises
        ------
        FileNotFoundError
            If the specified file does not exist.
        ValueError
            If the file format is invalid or cannot be parsed.

        Example
        -------
        >>> # Read standard XYZ file
        >>> positions = aesdebye.readXYZ("structure.xyz")
        >>> print(f"Loaded {positions.size()} atoms")
        >>>
        >>> # Read CSV-formatted coordinate file
        >>> positions = aesdebye.readXYZ("coords.csv", delimiter=",", skipLines=1)
        )pbdoc");

  m.def("makePeriodic", &makePeriodic);
}
