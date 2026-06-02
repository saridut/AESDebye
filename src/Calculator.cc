#include <utility>
#include <string>
#include <iostream>

#include <Calculator.hpp>
#include <ASF.hpp>

#pragma region mpi_reduction

/**
 * Adds the elements of two PDF vectors.
 *
 * This function takes PDF vectors, `in` and `inout`, and adds the corresponding elements
 * of `in` to `inout`. The result is stored in `inout`.
 *
 * @tparam T The type of the elements in the PDF vectors.
 * @param in Pointer to the input vectors.
 * @param inout Pointer to the vector where the result will be stored.
 * @param len Pointer to the length of the arrays.
 */
void addPDFVectorsMPI(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
    HistBin *inPtr = (HistBin *)in;
    HistBin *inoutPtr = (HistBin *)inout;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < *len; i++)
    {
        inoutPtr[i] = inPtr[i] + inoutPtr[i];
    }
}

/**
 * Reduces the PDF across all the ranks.
 *
 * This function reduces the local histogram across all the ranks using MPI_Allreduce.
 *
 * This is kept separate from the ParallelHelper so that we dont have ParallelHelper being
 * dependent on the PDF bin class.
 *
 * @param localHistogram The local histogram to be reduced.
 * @param parallelHelper The parallelHelper object.
 * @return The global histogram after reduction.
 */
std::vector<HistBin> reducePDF(std::vector<HistBin> &localHistogram, ParallelHelper &parallelHelper)
{
    // synchronize before starting the reduction
    parallelHelper.wait();

    // if not using MPI or there is only one rank: no need for reduction
    if (!(parallelHelper.worldSize > 1 && parallelHelper.useMPI))
        return localHistogram;

    double start = helpers::get_wall_time();

    // set up the MPI data types
    int numBlocks = sizeof(HistBin) / sizeof(int64);              // Number of blocks in the datatype
    MPI_Datatype MPI_Bin;
    MPI_Op MPI_binAdd;
    MPI_Datatype mpiType = MPI_LONG_LONG_INT; // MPI datatype

    MPI_Type_contiguous(numBlocks, mpiType, &MPI_Bin);        // create the contiguous datatype
    MPI_Type_commit(&MPI_Bin);                                // commit the datatype
    MPI_Op_create(addPDFVectorsMPI, 1, &MPI_binAdd); // create the custom MPI operation

    // create an empty global histogram
    std::vector<HistBin> globalHistogram(localHistogram.size());

    // reduce the local histogram across all the ranks
    MPI_Allreduce(localHistogram.data(), globalHistogram.data(), localHistogram.size(), MPI_Bin, MPI_binAdd, MPI_COMM_WORLD);

    // free the MPI data types
    MPI_Type_free(&MPI_Bin);
    MPI_Op_free(&MPI_binAdd);

    parallelHelper << "MPI reduction time: " << helpers::get_wall_time() - start << " s\n";

    // synchronize after the reduction - should be done implicitly but for good measure :)
    parallelHelper.wait();

    return globalHistogram;
}

#pragma endregion mpi_reduction

PDF DebyeCalculator::calculatePDF(Positions &positions)
{
    // Just pass on the same reference - self PDF calculation
    // this is important if the calculatePDF is being called from Python
    // as the positions will passed by value and copied into two different memory locations
    return calculatePDF(positions, positions);
}



PDF DebyeCalculator::calculatePDF(Positions &positionsI, Positions &positionsJ)
{
    bool samePositions = &positionsI == &positionsJ;
    parallelHelper.printSection("PDF calculation started");
    parallelHelper << ">> I(element:" << positionsI.element << ", Atoms:" << positionsI.size() <<
                        ") x " << (samePositions ?"I":"J") <<"(element:" << positionsJ.element
                        << ", Atoms:" << positionsJ.size() << ")\n";
    parallelHelper.wait();
    // make sure the box size is the same for the positions
    if (positionsI.getBoxSize() != positionsJ.getBoxSize())
    {
        parallelHelper << "Unit length mismatch! " << positionsI.getBoxSize() << " != "
                       << positionsJ.getBoxSize() << "\n";
        parallelHelper << "Setting the unit length to the larger of the two\n";

        // use the bigger box size
        double boxSize = std::max(positionsI.getBoxSize(), positionsJ.getBoxSize());

        // update the box size for both positions
        positionsI.setBoxSize(boxSize);
        positionsJ.setBoxSize(boxSize);
    }

    if (config.useCellList)
    {
        // setup cellList for the positions - this is required to be done **before** the binsResolution is set
        // otherwise the number of cells won't be appropriate!
        cellList.createCellList(positionsI);

        if (!samePositions)
            cellList.createCellList(positionsJ);
    }

    double originalBoxSize = positionsI.getBoxSize();

    // update the box size of both the positions if binsResolution is set
    if (binsResolution != 1.0)
    {
        parallelHelper << "Updating the box size for the bins resolution\n";
        parallelHelper << "Old box size: " << originalBoxSize << " , ";
        double boxSize = originalBoxSize * (1 / binsResolution);
        parallelHelper << "New box size: " << boxSize << "\n";
        positionsI.setBoxSize(boxSize);
        if (!samePositions)
            positionsJ.setBoxSize(boxSize);
    }

    // create the PDF object. Needs to know the boxSize and if both positions will be the same
    uint nBins = !config.smallBins ? N_BINS : N_SMALL_BINS;
    PDF pdf(positionsI.getBoxSize(), samePositions, nBins,
            positionsI.element, positionsJ.element);
    parallelHelper << "PDF has " << pdf.pdf_vector.size() << " bins\n";

    // Setup for dividing the work among the ranks
    // we use the triangle area to divide the work for same positions
    // for different positions, we simply divide the work equally
    int64 start, stop;
    if (samePositions)
    {
        stop = (int64)std::floor(std::sqrt((double)(parallelHelper.worldRank + 1) /
                                           (parallelHelper.worldSize)) *
                                 (double) positionsI.size());
        start = (int64)std::floor(std::sqrt((double)(parallelHelper.worldRank) /
                                            (parallelHelper.worldSize)) *
                                  (double) positionsI.size());
        parallelHelper << "thread: " << parallelHelper.worldRank << " PositionsI start: " << start << " , stop: " << stop << "\n";
    }
    else
    {
        int64 chunkSize = positionsJ.size() / parallelHelper.worldSize;
        start = parallelHelper.worldRank * chunkSize;
        stop = (parallelHelper.worldRank + 1) * chunkSize;
    }
    stop = parallelHelper.worldRank == parallelHelper.worldSize - 1 ? positionsI.size() : stop;

    double calculation_start = helpers::get_wall_time();
    if (config.useGPU)
    {
        calculatePDFGPU(positionsI, positionsJ, pdf.pdf_vector, cellList.getCellPairsList(positionsI, positionsJ, samePositions),
                        start, stop, samePositions, parallelHelper,
                        config);
    }
    else
    {
        if (config.useCellList)
        {
            // the pairs of cells that need to be calculated
            auto const cellPairsList = cellList.getCellPairsList(positionsI, positionsJ, samePositions);

            // setup the rank based divisions with efficient load balancing
            auto [rankStarts, rankEnds] = cellList.setupRankBasedDivisions(positionsI, positionsJ,
                                                                           cellPairsList, parallelHelper.worldSize,
                                                                           true);

            // get the start and stop for the current rank
            start = rankStarts[parallelHelper.worldRank];
            stop = rankEnds[parallelHelper.worldRank];
            if (verbose && parallelHelper.worldSize > 1)
            {
                std::cout << ">> world rank: " << parallelHelper.worldRank << " , ";
                std::cout << "start: " << start << " , stop: " << stop << " , ";
                std::cout << "cellPairs: " << cellPairsList.size() << " , ";
                std::cout << "cells: " << 100.0 * double(stop - start) / (double)cellPairsList.size() << "% \n";
            }
            calculation_start = helpers::get_wall_time();

            // calculate the PDF - notice how we only pass the vector datastructure
            if (config.useLocalHistogram)
            {
                calculatePDFKernel<true>(positionsI, positionsJ,
                                         cellPairsList, pdf.pdf_vector, start, stop, samePositions);
            }
            else
            {
                calculatePDFKernel<false>(positionsI, positionsJ,
                                          cellPairsList, pdf.pdf_vector, start, stop, samePositions);
            }
        }
        else
        {
            parallelHelper << "Non-cell list calculation\n";
            calculation_start = helpers::get_wall_time();
            // non cell list based PDF calculation
            if (config.useLocalHistogram)
            {
                calculatePDFKernel<true>(positionsI, positionsJ,
                                         pdf.pdf_vector, start, stop, samePositions);
            }
            else
            {
                calculatePDFKernel<false>(positionsI, positionsJ,
                                          pdf.pdf_vector, start, stop, samePositions);
            }
        }
    }
    double calculation_end = helpers::get_wall_time();
    pdf.calculationTime = calculation_end - calculation_start;

    // reduce the pdfs accross all the ranks
    pdf.pdf_vector = reducePDF(pdf.pdf_vector, parallelHelper);

    // self pairs - It will check if same positions and add self pairs
    pdf.addSelfPairs(positionsI.size()); // self pairs are equal to the number of atoms

    // create pdf vectors for access
    pdf.createPDFVectors();

    // restore box size if binsResolution is set
    updateBoxSize(positionsI, positionsJ, samePositions, true);

    // bool test = testHistogram(pdf.pdf_vector, positionsI, positionsJ);
    parallelHelper.wait();
    pdf.test(positionsI.size(), positionsJ.size());

    // output
    uint64 nPd = samePositions ? positionsI.size() * (positionsI.size() - 1) / 2 : positionsI.size() * positionsJ.size();
    auto pdPerSecond = (uint64) ((double) nPd / pdf.calculationTime / 1e6);
    uint64 pdPerSecondPerThread = pdPerSecond / parallelHelper.totalThreads;
    std::string output = "Complete! time: " + std::to_string(pdf.calculationTime) + " s" +
                         " ," + std::to_string(pdPerSecondPerThread) +"/" + std::to_string(pdPerSecond) + " MPd/s";
    parallelHelper.printSection(output);
    return pdf;
}

void
DebyeCalculator::updateBoxSize(Positions &positionsI, Positions &positionsJ, bool samePositions,
                                      bool scaleUp) {
    double originalBoxSize = positionsI.getBoxSize();
    if (binsResolution != 1.0)
    {
        double boxSize = originalBoxSize * (scaleUp ? binsResolution : (1 / binsResolution));
        positionsI.setBoxSize(boxSize);
        if (!samePositions)
            positionsJ.setBoxSize(boxSize);
        parallelHelper << ">> Changed box size from: " << originalBoxSize << " to " << boxSize << "\n";
    }
}

std::vector<double>
DebyeCalculator::calculateIntensity(std::vector<double> const &centers, const std::vector<double> &counts,
                                    std::vector<double> const &qVector, std::string elementI, std::string elementJ)
{
    // call the kernel to calculate the intensity
    std::vector<double> intensity;
    if (config.useGPU && (qVector.size() > 5000 || parallelHelper.ompThreads <= 4)) // either too large intensity, or too less cpus
    {
        parallelHelper << ">> Intensity calculation on GPU\n";
        intensity = calculateIntensityGPU(qVector, centers, counts);
    }
    else
    {
        parallelHelper << ">> Intensity calculation on CPU\n";
        intensity = calculateIntensityCPU(qVector, centers, counts);
    }

    // Add the ASF contribution if the elements are provided
    if (elementI.empty() && elementJ.empty())
    {
        return intensity;
    }
    else if (elementJ.empty())
    {
        elementJ = elementI;
    }

    auto asfProfileI = calculateASFProfile(qVector, elementI);
    auto asfProfileJ = elementI == elementJ? asfProfileI : calculateASFProfile(qVector, elementJ);

#pragma omp parallel for schedule(static) shared(qVector, intensity, asfProfileI, asfProfileJ) default(none)
    for (size_t i = 0; i < qVector.size(); i++)
    {
        intensity[i] *= asfProfileI[i] * asfProfileJ[i];
    }

    return intensity;
}

Profile
DebyeCalculator::calculateIntensity(PDF &pdf, double start, double end, int nSteps, bool twoThetaSpace, double lambda)
{
    parallelHelper.wait();
    Profile profile(start, end, nSteps, twoThetaSpace, lambda);
    double calculation_start = helpers::get_wall_time();
    profile.intensity = calculateIntensity(pdf.centers, pdf.counts, profile.q, pdf.elementI, pdf.elementJ);
    double calculation_end = helpers::get_wall_time();
    profile.calculationTime = calculation_end - calculation_start;
    profile.test();
    return profile;
}

std::map<std::string, std::pair<PDF, Profile>>
DebyeCalculator::calculateProfile(Positions &positions, double start, double end, int nSteps,
                                  bool twoThetaSpace, double lambda, std::string filter)
{
    std::vector<std::string> positionPairs;
    std::map<std::string, std::pair<PDF, Profile>> results;

    if (filter.empty())
    {
        auto uniqueElements = positions.getUniqueElements();
        for (size_t i = 0; i < uniqueElements.size(); i++)
        {
            for (size_t j = i; j < uniqueElements.size(); j++)
            {
                positionPairs.push_back(uniqueElements[i] + "-" + uniqueElements[j]);
            }
        }
    }
    else
    {
        parallelHelper << "Filtering with: " << filter << "\n";
        positionPairs = helpers::stringSplit(filter, ",");
    }

    for (const auto& pairName : positionPairs)
    {
        std::string iType = helpers::stringSplit(pairName, "-")[0];
        std::string jType = helpers::stringSplit(pairName, "-")[1];
        Positions positionsI = positions.filterByElement(iType);
        PDF pdf;
        if (iType == jType) // faster calculation for same positions due to symmetry
        {
            pdf = calculatePDF(positionsI);
        }
        else
        {
            Positions positionsJ = positions.filterByElement(jType);
            pdf = calculatePDF(positionsI, positionsJ);
        }
        Profile profile = calculateIntensity(pdf, start, end, nSteps, twoThetaSpace, lambda);
        results.emplace(pairName, std::make_pair(std::move(pdf), std::move(profile)));
    }

    parallelHelper.wait();

    auto fullPDF = results[positionPairs[0]].first;
    auto fullProfile = results[positionPairs[0]].second;
    for (size_t i = 1; i < positionPairs.size(); i++)
    {
        fullProfile = fullProfile + results[positionPairs[i]].second;
        fullPDF = fullPDF + results[positionPairs[i]].first;
    }
    results.emplace("total", std::make_pair(std::move(fullPDF), std::move(fullProfile)));

    return results;
}



std::vector<double>
DebyeCalculator::calculateASFProfile(const std::vector<double> &qVector,
                                     std::string elementI)
{
    ASFCoeffs asf = ASFTable[elementI];
    std::vector<double> intensity(qVector.size(), 1.0);
#pragma omp parallel for schedule(static) shared(qVector, intensity, asf) default(none)
    for (size_t qIdx = 0; qIdx < qVector.size(); qIdx++)
    {
        double q = qVector[qIdx];
        intensity[qIdx] = asf.ASFq(q);
    }
    return intensity;
}


