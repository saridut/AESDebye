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

    void updateBoxSize(Positions &positionsI, Positions &positionsJ, bool samePositions, bool scaleUp);

public:
    ParallelHelper parallelHelper; //< needs to be public for export to python @see ParallelHelper
    CellList cellList; //< CellList object for cell list computation. @see CellList
    CalculationConfig config; //< Configuration for the histogram computation. @see CalculationConfig, just a boolean container

    /**
     * @brief Constructs a DebyeCalculator object.
     *
     * This constructor initializes a DebyeCalculator object with the specified parameters.
     *
     * @param nThreads The number of threads to use for parallel computation. Default is -1, which means to use the maximum available threads.
     * @param nCells The number of cells to use for cell list computation. Default is 15. If you don't want to use cell list (only recommended for very small MD systems), set this to 1.
     * @param useMPI Flag indicating whether to use MPI for parallel computation. Default is false. When set to true, the runtime will need to use an mpi-runner to run the code. MPI library that was used for compilation should be loaded in the runtime environment.
     * @param useGPU Flag indicating whether to use GPU for parallel computation. Default is false.
     * @param useLocalHistogram Flag indicating whether to use local histogram for parallel computation. Default is true. Set to false for crystalline structures.
     * @param verbose Flag indicating whether to enable verbose output. Default is true.
     */
    explicit DebyeCalculator(int nThreads = -1, int nCells = 15,
                             double binsResolution = 1.0,
                             bool useMPI = false,
                             bool useGPU = false,
                             bool useLocalHistogram = true,
                             bool smallBins = false,
                             bool pseudoCoal = true,
                             bool fillGPU = false,
                             bool verbose = true) : binsResolution(binsResolution), verbose(verbose),
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
            std::cout << "bins resolution: "<< binsResolution << std::endl;
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

        parallelHelper.printSection("DebyeCalculator initialized");
        parallelHelper << ">> bins resolution: " << binsResolution << "\n";
        parallelHelper << config;

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

    std::map<std::string, std::pair<PDF, Profile>>
    calculateProfile(Positions &positions, double start, double end, int nSteps,
                     bool twoThetaSpace = false, double lambda = 0.4, std::string filter = "");

    std::map<std::string, std::pair<PDF, Profile>>
    calculateProfileBruteForce(Positions &positions, double start, double end, int nSteps,
                               bool twoThetaSpace = false, double lambda = 0.4, std::string filter = "");

//    std::vector<double>
//    addASFContribution(const std::vector<double> &qVector,
//        const std::vector<double> &intensity,
//        std::string elementI="", std::string elementJ="");

    static std::vector<double>
    calculateASFProfile(const std::vector<double> &qVector,
                        std::string elementI);

    void setVerbosity(bool verbose)
    {
        this->verbose = verbose;
        parallelHelper.verbose = verbose;
    }

    // PDF reducePDF(const PDF &pdf, Profile &profile);
    // std::vector<double> reducePDF(const std::vector<double> &centers, const std::vector<double> &counts,
    //                               const std::vector<double> &qVector, const std::vector<double> &intensity,
    //                               uint64 size);



};