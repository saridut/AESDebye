/**
 * @file Calculator.hpp
 * @brief Contains the declaration of the DebyeCalculator class.
 *
 * This file contains the declaration of the DebyeCalculator class, which is the main calculator class for DSE calculations.
 * It provides methods for calculating the pair distribution function (PDF) and intensity profile using DSE.
 */

#pragma once

#include <optional>

#include <DataTypes.hpp>
#include <IntensityKernels.hpp>
#include <IntensityKernelsGPU.hpp>
#include <ParallelHelper.hpp>
#include <PDF.hpp>
#include <PDFKernels.hpp>
#include <PDFKernelsGPU.hpp>
#include <Positions.hpp>
#include <Profile.hpp>
#include <utils.hpp>
#include <ASF.hpp>

/**
 * @class DebyeCalculator
 * @brief Main calculator class for DSE calculations
 *
 * The DebyeCalculator class provides methods for calculating the pair distribution function (PDF),
 * intensity profile using DSE. This class is calls the appropriate computation kernels.
 */
class DebyeCalculator
{
private:
    bool verbose = false;        //< Flag indicating whether to enable verbose output.
    double binsResolution = 1.0; //< The resolution of the bins in the histogram. goes from (0, 1.0]

    /**
     * @brief Updates the box size of the positions based on binsResolution.
     *
     * Scales the box size of positionsI (and positionsJ if different) to adjust the binning resolution.
     *
     * @param positionsI The first positions of the particles.
     * @param positionsJ The second positions of the particles.
     * @param samePositions Flag indicating if positionsI and positionsJ refer to the same object.
     * @param scaleUp Flag indicating whether to scale the box size (true) or restore it to its original size (false).
     */
    void updateBoxSize(Positions &positionsI, Positions &positionsJ, bool samePositions, bool scaleUp);


public:
    ParallelHelper parallelHelper; //< needs to be public for export to python @see ParallelHelper
    CellList cellList; //< CellList object for cell list computation. @see CellList
    CalculationConfig config; //< Configuration for the histogram computation. @see CalculationConfig, just a boolean container
    ASFTable asfTable; //< Atomic scattering factor table for element form factors

    /**
     * @brief Constructs a DebyeCalculator object.
     *
     * This constructor initializes a DebyeCalculator object with the specified parameters.
     *
     * @param nThreads The number of threads to use for parallel computation. Default is -1, which means to use the maximum available threads.
     * @param nCells The number of cells to use for cell list computation. Default is 15. If you don't want to use cell list (only recommended for very small MD systems), set this to 1.
     * @param binsResolution The resolution of the bins in the histogram.
     * @param useMPI Flag indicating whether to use MPI for parallel computation. Default is false. When set to true, the runtime will need to use an mpi-runner to run the code. MPI library that was used for compilation should be loaded in the runtime environment.
     * @param useGPU Flag indicating whether to use GPU for parallel computation. Default is false.
     * @param useLocalHistogram Flag indicating whether to use local histogram for parallel computation. Default is true. Set to false for crystalline structures.
     * @param smallBins Use smaller bins for PDF. Default is false.
     * @param pseudoCoal Enable pseudo-coalescence optimization. Default is true.
     * @param fillGPU Fill GPU memory completely. Default is false.
     * @param verbose Flag indicating whether to enable verbose output. Default is true.
     * @param asfFilePath Path to custom ASF table file (TSV/CSV/text). If empty, standard data/asf.tsv is loaded.
     * @param asfFormulation Formulation standard to use: "WaasmaierKirfel5" (default, 5 Gaussians) or "CromerMann4" (4 Gaussians).
     */
    explicit DebyeCalculator(int nThreads = -1, int nCells = 15,
                             double binsResolution = 1.0,
                             bool useMPI = false,
                             bool useGPU = false,
                             bool useLocalHistogram = true,
                             bool smallBins = false,
                             bool pseudoCoal = true,
                             bool fillGPU = false,
                             bool verbose = true,
                             const std::string& asfFilePath = "",
                             const std::string& asfFormulation = "WaasmaierKirfel5") : binsResolution(binsResolution), verbose(verbose),
                                                    parallelHelper(nThreads, useMPI, verbose)
    {

#ifndef USE_GPU
        if (useGPU)
        {
            throw std::runtime_error("GPU support is not enabled in the current build. Please recompile with GPU support.");
        }
#endif

        if (binsResolution > 1.0 || binsResolution <= 0)
        {
            throw std::invalid_argument("Bins resolution should be between 0 and 1.0");
        }

        config.useCellList = nCells > 1;
        config.useMPI = useMPI;
        config.useGPU = useGPU;
        config.useLocalHistogram = useLocalHistogram;
        config.smallBins = smallBins;
        config.pseudoCoal = pseudoCoal;
        config.fillGPU = fillGPU;
        config.useGPUCellList = false;

        // Initialize ASF Table with the explicit formulation
        ASFFormulation form = stringToASFFormulation(asfFormulation);
        std::string targetASFFile = asfFilePath.empty() ? ASFTable::DEFAULT_ASF_FILE : asfFilePath;
        try
        {
            asfTable.loadFromFile(targetASFFile, form);
        }
        catch (const std::exception& e)
        {
            parallelHelper << ">> [DebyeCalculator] Error loading ASF file (" << targetASFFile 
                           << "): " << e.what() << "\n";
            throw;
        }

        parallelHelper.printSection("DebyeCalculator initialized");

        parallelHelper << ">> Please cite:\n\n";
        parallelHelper << ">> Panchi, N., Kuckuk, S., Wittmann, M., Engel, M., & Leonardi, A. (2026).\n"
                       << ">> AES-Debye: An Accurate, Efficient and Scalable Engine for Debye Scattering Calculations.\n"
                       << ">> Journal of Applied Crystallography, 59(5), 1478–1490. https://doi.org/10.1107/S1600576726007429\n\n";
        parallelHelper << ">> BibTeX:\n\n";
        parallelHelper << "@article{panchi_aes-debye_2026,\n"
                       << "  title   = {{AES-Debye}: An Accurate, Efficient and Scalable Engine for {Debye} Scattering Calculations},\n"
                       << "  author  = {Panchi, Navid and Kuckuk, Sebastian and Wittmann, Markus and Engel, Michael and Leonardi, Alberto},\n"
                       << "  journal = {Journal of Applied Crystallography},\n"
                       << "  volume  = {59},\n"
                       << "  number  = {5},\n"
                       << "  pages   = {1478--1490},\n"
                       << "  year    = {2026},\n"
                       << "  month   = {oct},\n"
                       << "  issn    = {1600-5767},\n"
                       << "  doi     = {10.1107/S1600576726007429}\n"
                       << "}\n";
        parallelHelper << " ------------------------------------------------\n";
        parallelHelper << config;
        parallelHelper << ">> ASF Formulation: " << asfFormulationToString(asfTable.getFormulation())
                       << " (" << asfTable.size() << " elements loaded)\n";
        if (!asfTable.getLoadedFilePath().empty())
        {
            parallelHelper << ">> ASF File: " << asfTable.getLoadedFilePath() << "\n";
        }
        parallelHelper << " ------------------------------------------------\n";

        // if cell list is enabled, initialize it
        if (config.useCellList)
        {
            cellList = CellList(nCells);

            parallelHelper << ">> CellList initialized with " << nCells << " cells\n";
            parallelHelper << ">> Setting up cell pairs ...";
            double start_time = helpers::get_wall_time();

            // create the cell pairs - this is a one-time operation, 
            // will be used for all PDF calculations
            cellList.createCellPairs(true);
            parallelHelper << " complete! Took " << helpers::get_wall_time() - start_time << " s\n";
        }
    };

    /**
     * @brief Calculate the pair distribution function (PDF) for the given positions.
     *
     * Creates two identical references to the `Positions` object and forwards it to the `calculatePDF` method.
     * Acts as a helper in CPP, but very useful in Python bindings since when calling from Python, `Positions` object is passed by value - so a copy is created.
     *
     * @param positions The positions of the particles.
     * @return A tuple containing the PDF, the time taken to calculate the PDF, and a boolean indicating whether the test passed.
     */
    PDF calculatePDF(Positions &positions);

    /**
     * @brief Calculate the pair distribution function (PDF) for the given positions.
     *
     * This method calculates the pair distribution function (PDF) for the given positions.
     * The core setup logic is done here, the actual computations are done by the kernels.
     * Look at the `PDFCPUKernels` class for the actual computation logic.
     *
     * @param positionsI The first positions of the particles.
     * @param positionsJ The second positions of the particles.
     * @return A tuple containing the PDF, the time taken to calculate the PDF, and a boolean indicating whether the test passed.
     */
    PDF calculatePDF(Positions &positionsI, Positions &positionsJ);



    /**
     * @brief Calculates the intensity for a given range of q values.
     *
     * @param qVector A vector of q values.
     * @param centers A vector of centers. Use the corrected centers from the PDF by default to get the corrections.
     * @param counts A vector of counts.
     * @return A tuple containing the calculated intensities, the corresponding q values, and a boolean indicating success.
     */
    std::vector<double> calculateIntensity(std::vector<double> const &centers, const std::vector<double> &counts,
                                           std::vector<double> const &qVector,
                                           std::string elementI= "", std::string elementJ = "");

    /**
     * @brief Calculates the intensity for a given range of q values.
     *
     * @param qStart The starting value of q.
     * @param qEnd The ending value of q.
     * @param nqSteps The number of nSteps between qStart and qEnd.
     * @param histogram The PDF histogram used for intensity calculation.
     * @return A tuple containing two vectors of doubles representing the calculated intensity values and a boolean indicating the success of the calculation.
     */
    Profile
    calculateIntensity(PDF &pdf, double start, double end, int nSteps, 
    bool twoThetaSpace = false, double lambda = 0.4);

    /**
     * @brief Calculates both the PDF and intensity profile for different atom-pair combinations.
     *
     * If filter is empty, it automatically computes profiles for all unique element pairs.
     * Otherwise, it processes only the pairs specified in the filter string (comma-separated, e.g., "Fe-Fe,Fe-O").
     * Returns a map containing the calculated PDF and Profile pairs, including a "total" combined profile.
     *
     * @param positions The positions of the particles.
     * @param start The starting value of the range (q or 2theta).
     * @param end The ending value of the range (q or 2theta).
     * @param nSteps The number of steps between start and end.
     * @param twoThetaSpace Flag indicating if the range is in 2theta space (true) or q space (false). Default is false.
     * @param lambda The X-ray wavelength to use for 2theta conversions. Default is 0.4.
     * @param filter Filter string to specify element pairs to calculate. Default is "".
     * @return A map with pair names (e.g., "Fe-Fe", "total") as keys and their corresponding PDF and Profile as values.
     */
    std::map<std::string, std::pair<PDF, Profile>>
    calculateProfile(Positions &positions, double start, double end, int nSteps,
                     bool twoThetaSpace = false, double lambda = 0.4, std::string filter = "");

    /**
     * @brief Calculates the Atomic Scattering Factor (ASF) profile for a given element.
     *
     * Computes the ASF for each q value in the provided qVector using the calculator's
     * configured ASF table.
     *
     * @param qVector A vector of q values.
     * @param elementI The chemical symbol of the element.
     * @return A vector containing the calculated ASF values.
     */
    std::vector<double>
    calculateASFProfile(const std::vector<double> &qVector,
                        std::string elementI) const;

    /**
     * @brief Validates that a given element exists in the ASF table, or is 'None' / empty.
     * Throws std::invalid_argument if the element is not found.
     *
     * @param element Chemical symbol of the element.
     */
    void validateElement(const std::string& element) const
    {
        if (element.empty() || element == "None" || element == "none" || element == "NONE")
        {
            return;
        }
        if (!asfTable.contains(element))
        {
            throw std::invalid_argument("Element '" + element + "' not found in ASF table (file: " +
                                        asfTable.getLoadedFilePath() + ")");
        }
    }

    /**
     * @brief Validates that all elements present in a Positions object exist in the ASF table.
     *
     * @param positions Positions object to check.
     */
    void validateElements(const Positions& positions) const
    {
        for (const auto& elem : positions.getUniqueElements())
        {
            validateElement(elem);
        }
    }

    /**
     * @brief Sets the verbosity flag for the calculator and parallelHelper.
     *
     * @param verbose Flag indicating whether to enable verbose output.
     */
    void setVerbosity(bool verbose)
    {
        this->verbose = verbose;
        parallelHelper.verbose = verbose;
    }

};
