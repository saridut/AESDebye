#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/Calculator.hpp"

namespace py = pybind11;

#define STRINGIFY(x) #x

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
  version.attr("version") = STRINGIFY(VERSION_INFO);
  version.attr("compiler") = STRINGIFY(CXX_COMPILER);

  // options to make the docstring more informative
  m.doc() = R"pbdoc(
        PyAESDebye
        )pbdoc";

  py::options options;
  options.disable_function_signatures();

  py::class_<Profile>(
      m, "Profile", py::module_local()) // cannot initialize in python directly
      .def("toCSV", &Profile::toCSV, py::arg("filename"), R"pbdoc(
        toCSV(filename: str) -> None

        Save the profile to a CSV file.

        Parameters
        ----------
        filename : str
            The filename to save the profile to.
        )pbdoc")

      // operators
      .def(py::self + py::self)
      .def(py::self - py::self)

      // __str__ method
      .def("__repr__", &Profile::toString)

      // props
      .def_readonly("q", &Profile::q, "List of q values")
      .def_readonly("twoTheta", &Profile::twoTheta, "List of two theta values")
      .def_readonly("intensity", &Profile::intensity,
                    "List of intensity values")
      .def_readonly("calculationTime", &Profile::calculationTime,
                    "Time taken for the calculation")
      .def_readonly("testPassed", &Profile::testPassed,
                    "Flag to indicate if the test passed");

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

  py::class_<Positions>(m, "Positions", py::module_local())
      .def(
          py::init<std::vector<std::string>, std::vector<std::array<double, 3>>,
                   std::vector<std::string>, double, double>(),
          py::arg("chemicalSymbols"), py::arg("positions"),
          py::arg("selectionIds") = std::vector<std::string>(),
          py::arg("boxMin") = 1e100, py::arg("boxMax") = -1e100)
      .def_readwrite("element", &Positions::element)
      .def_readwrite("selectionIds", &Positions::selectionIds,
                     "List of atom ids")
      .def_readwrite("chemicalSymbols", &Positions::chemicalSymbols,
                     "List of atom elementss")
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

  py::class_<PDF>(m, "PDF", py::module_local())
      .def(py::init<double, bool>(), py::arg("boxSize") = 1,
           py::arg("samePositions") = true)
      .def(py::self += py::self)
      .def("__repr__", &PDF::toString)
      .def("save", &PDF::save, py::arg("filename"), R"pbdoc(
                save(filename: str) -> None

                Save the PDF to a binary file.

                Parameters
                ----------
                filename : str
                    The path to the output file.
                )pbdoc")
      .def("load", &PDF::load, py::arg("filename"), R"pbdoc(
                load(filename: str) -> float

                Load the PDF from a binary file.

                Parameters
                ----------
                filename : str
                    The path to the input file.

                Returns
                -------
                float
                    The box size.
                )pbdoc")
      .def("toCSV", &PDF::toCSV, py::arg("filename"),
           py::arg("complete") = false, R"pbdoc(
                toCSV(filename: str, complete: bool = False) -> None

                Save the PDF to a CSV file.

                Parameters
                ----------
                filename : str
                    The path to the output file.
                complete : bool, optional
                    Flag to indicate if the complete PDF should be saved. Saving complete PDF can take up more memory. The default is False.
                )pbdoc")
      .def("readCSV", &PDF::readCSV, R"pbdoc(
                readCSV(filename: str) -> Tuple[List[float], List[float], List[int], List[float]]

                Read a PDF from a CSV file.

                Parameters
                ----------
                filename : str
                    The path to the input file.

                Returns
                -------
                Tuple[List[float], List[float], List[int], List[float]]
                    Tuple containing the bin centers, corrected bin centers, bin counts, and delta.
                )pbdoc")
      // properties
      .def_property_readonly("uncorrectedCenters", &PDF::getUncorrectedCenters,
                             "List of uncorrected bin centers")
      .def_property_readonly("centers", &PDF::getCenters, "List of bin centers")
      .def_property_readonly("counts", &PDF::getCounts, "List of bin counts")

      .doc() = R"pbdoc(
                PDF(boxSize: float, samePositions: bool = True) -> None

                Construct a new PDF object.

                Parameters
                ----------
                boxSize : float
                        The size of the box.
                samePositions : bool, optional
                        Flag to indicate if this is a self-PDF are the same. The default is True.
                )pbdoc";

  py::class_<DebyeCalculator>(m, "DebyeCalculator", py::module_local())
      .def(py::init<int, int, double, bool, bool, bool, bool, bool, bool,
                    bool>(),
           py::arg("nThreads") = -1, py::arg("nCells") = 15,
           py::arg("binsResolution") = 1.0, py::arg("useMPI") = false,
           py::arg("useGPU") = false, py::arg("useLocalHist") = true,
           py::arg("smallBin") = false, py::arg("pseudoCoal") = true,
           py::arg("fillGPU") = false, py::arg("verbose") = true)
      .def("calculatePDF",
           py::overload_cast<Positions &, Positions &>(
               &DebyeCalculator::calculatePDF),
           py::arg("positionsI"), py::arg("positionsJ").none(false), R"doc(
                calculatePDF(positionsI: Positions, positionsJ: Positions = None) -> PDF

                Calculate the PDF for the given positions. If two positions are given, the PDF between the two positions is calculated.
                Else, the self-PDF is calculated.

                Warning
                -------
                Do not provide None for the second argument, this will result in an error.

                Note
                -------
                For self-PDF, you can simply provide a single positions object.

                Parameters
                ----------
                positionsI : Positions
                    The first positions object.

                positionsJ : Positions, optional
                    The second positions object.

                Returns
                -------
                PDF
                    The PDF calculated between the two positions objects.
                )doc")

      .def("calculatePDF",
           py::overload_cast<Positions &>(&DebyeCalculator::calculatePDF),
           py::arg("positionsI"))

      .def("calculateIntensity",
           py::overload_cast<
               std::vector<double> const &, std::vector<double> const &,
               std::vector<double> const &, std::string, std::string>(
               &DebyeCalculator::calculateIntensity),
           py::arg("centers"), py::arg("counts"), py::arg("qVector"),
           py::arg("elementI") = "", py::arg("elementJ") = "", R"pbdoc(
                calculateIntensity(start: float, end: float, nSteps: int, pdf:PDF=None, centers: List[float]=None, counts: List[int]=None) -> Tuple[List[float], List[float], List[int], List[float]

                Calculate the intensity for a given set of centers and counts.

                Parameters
                ----------
                start : float
                    The starting q/theta value.
                end : float
                        The ending q/theta value.
                nSteps : int
                        The number of q/theta steps.
                pdf : PDF, optional
                        The PDF object. Either provide this, or provide the centers and counts. The default is None.
                centers : List[float], optional
                        List of bin centers. These bin centers are generally the corrected bin centers. Only provide this if pdf is not provided.
                counts : List[int], optional
                        List of bin counts. Only provide this if pdf is not provided.

                Warning
                -------
                Either provide the PDF object or the centers and counts. Providing both will result in an error. Recommended way is to use to PDF.
                This calculateIntensity function has two overloads, one for each case.


                Returns
                -------
                Tuple[List[float], List[float], bool]
                    Tuple containing a list of q values, intensity values, and a flag indicating if the intensity is > 0.

                )pbdoc")

      .def("calculateIntensity",
           py::overload_cast<PDF &, double, double, int, bool, double>(
               &DebyeCalculator::calculateIntensity),
           py::arg("pdf"), py::arg("start"), py::arg("end"), py::arg("nSteps"),
           py::arg("twoThetaSpace") = false, py::arg("wavelength") = 0.4)

      .def("calculateProfile",
           py::overload_cast<Positions &, double, double, int, bool, double,
                             std::string>(&DebyeCalculator::calculateProfile),
           py::arg("Positions"), py::arg("start") = 0, py::arg("end") = 10,
           py::arg("steps") = 1000, py::arg("twoTheta") = false,
           py::arg("wavelength") = .4, py::arg("filter") = "", R"pbdoc(
                calculateProfile(positions: Positions, start: float = 0, end: float = 10, nSteps: int = 1000) -> Tuple[List[float], List[float], List[int], List[float]

                Calculate the intensity profile for a given set of positions.

                Parameters
                ----------
                positions : Positions
                    The positions object.
                start : float, optional
                    The starting q value. The default is 0.
                end : float, optional
                    The ending q value. The default is 10.
                nSteps : int, optional
                    The number of q steps. The default is 1000.

                Returns
                -------
                Tuple[List[float], List[float], PDF]
                    Tuple containing a list of q values, intensity values, and the PDF object.
                )pbdoc")
      .def_static("calculateASFProfile", &DebyeCalculator::calculateASFProfile,
                  py::arg("qVector"), py::arg("element"))
      .def_readonly("parallelHelper", &DebyeCalculator::parallelHelper,
                    "ParallelHelper object created by the DebyeCalculator")
      .def_readonly("cellList", &DebyeCalculator::cellList,
                    "CellList object created by the DebyeCalculator")
      .doc() = R"pbdoc(
                DebyeCalculator(nThreads: int = -1, nCells: int = 15, useMPI: bool = False, useGPU: bool = False, useLocalHist: bool = True, verbose: bool = True) -> None

                The main calculator object. All the computations with the library are done using this object.

                Parameters
                ----------
                nThreads : int, optional
                        The number of threads to use. The default is -1, which means all available threads.
                nCells : int, optional
                        The number of cells (in each dim) to use to partition the space. The default is 15.
                binsResolution : float, optional
                        The resolution of the bins. The default is 1.0. Has to be in range (0, 1.0].
                useMPI : bool, optional
                        Flag to indicate if MPI should be used, When this is True, the Python script must be called with an MPI runner. The default is False.
                useGPU : bool, optional
                        Flag to indicate if GPU should be used. This will result in a runtime error if the GPU version is not compiled. The default is False.
                useLocalHist : bool, optional
                        Flag to indicate if local histograms should be used. Keep this True except in cases where the computations need to be done on a crystalline sample. The default is True.
                verbose : bool, optional
                        Flag to indicate if verbose output should be enabled. The default is True.
                )pbdoc";

  m.def("generateData", &generateTestData, py::arg("lattice"),
        py::arg("nRepeats"), py::arg("noise") = 0.0, py::arg("seed") = 1,
        py::arg("element") = "Pt",
        R"pbdoc(
        generateData(lattice: float, nRepeats: int, noise: float, seed: int) -> Positions
        Generate FCC lattices with a given number of repeats in in each direction.

        Args:
            lattice (float): Lattice constant
            nRepeats (int): Number of repeats in each direction
            noise (float): Noise level
            seed (int): Seed for random number generator
            element (str): Element type to use for the positions, default is "Pt"

        Returns:
                Positions: Positions object containing the generated data
        )pbdoc");

  m.def("readXYZ", &readXYZ, py::arg("filename"), py::arg("delimiter") = " ",
        py::arg("skipLines") = 2, py::arg("typeMapping") = "0:None", R"pbdoc(
        readXYZ(filename: str) -> Positions
        Read a single frame XYZ file. Expects to have 4 columns with atom id, x, y, z and atom type.

        Args:
            filename (str): Path to the XYZ file

        Returns:
                Positions: Positions object containing the read data
        )pbdoc");

  m.def("makePeriodic", &makePeriodic);
}
