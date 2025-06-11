#include <Calculator.hpp>
#include <argparse.hpp>

int main(int argc, char *argv[])
{
    int nRepeats, nCells, nThreads, nSteps;
    double stdDev, binsResolution, start=0, end=20, lambda=0.4;
    bool useMPI, useGPU, dontUseLocalHistogram, nonVerbose, smallBins, pseudoCoal, 
    fillGPU, useGPUCellList, twoTheta;
    bool shufflePositions;
    std::string inputFilename = "", outputDir = "", typeMapping = "0:None";
    argparse::ArgumentParser program("AESDebye");
    program.add_argument("-nr", "--nRepeats")
        .scan<'d', int>() // 'd' is the type of the argument, int in this case
        .default_value(-1)
        .help("Number of repeats in the lattice").store_into(nRepeats);

    program.add_argument("-s", "--stdDev")
        .scan<'g', double>()
        .help("Std dev of noise in the lattice")
        .default_value(0.0).store_into(stdDev);

    program.add_argument("-b", "--binsResolution")
        .scan<'g', double>()
        .help("Resolution of the bins")
        .default_value(1.0).store_into(binsResolution);

    program.add_argument("-sb", "--smallBins")
        .help("Use small bins for the PDF")
        .default_value(false)
        .implicit_value(true).store_into(smallBins);

    program.add_argument("-nc", "--nCells")
        .scan<'d', int>()
        .help("Number of cells in the lattice")
        .default_value(15).store_into(nCells);

    program.add_argument("-nt", "--nThreads")
        .scan<'d', int>() // 'd' is the type of the argument, int in this case
        .default_value(-1)
        .help("Number of threads to use").store_into(nThreads);


    program.add_argument("-mpi", "--useMPI")
        .help("Use MPI for parallelization")
        .default_value(false)
        .implicit_value(true).store_into(useMPI);

    program.add_argument("-gpu", "--useGPU")
        .help("Use GPU for parallelization")
        .default_value(false)
        .implicit_value(true).store_into(useGPU);

    program.add_argument("-nlh", "--dontUseLocalHistogram")
        .help("dont use local histogram for noisy calcs")
        .default_value(false)
        .implicit_value(true).store_into(dontUseLocalHistogram);

    program.add_argument("-v", "--nonVerbose")
        .help("Dont enable verbose output")
        .default_value(true)
        .implicit_value(true).store_into(nonVerbose);

    program.add_argument("-f", "--inputFilename")
        .help("Filename to save the data")
        .default_value("").store_into(inputFilename);

    program.add_argument("-o", "--outputDir")
        .help("Output directory to save the data")
        .default_value("").store_into(outputDir);

    program.add_argument("-sp", "--shufflePositions")
        .help("Shuffle the positions")
        .default_value(false)
        .implicit_value(true).store_into(shufflePositions);

    program.add_argument("-pc", "--pseudoCoal")
        .help("Use pseudo coal")
        .default_value(false)
        .implicit_value(true).store_into(pseudoCoal);

    program.add_argument("-fg", "--fillGPU")
        .help("Fill the GPU with threads")
        .default_value(false)
        .implicit_value(true).store_into(fillGPU);

    program.add_argument("-gcl", "--useGPUCellList")
        .help("Use GPU cell list")
        .default_value(false)
        .implicit_value(true).store_into(useGPUCellList);

    program.add_argument("-tm", "--typeMapping")
        .help("Type mapping for the elements")
        .default_value("0:None").store_into(typeMapping);

    program.add_argument("-st", "--start")
        .scan<'g', double>()
        .help("Start of the q/theta range")
        .default_value(0.0).store_into(start);

    program.add_argument("-e", "--end")
        .scan<'g', double>()
        .help("End of the q/theta range")
        .default_value(20.0).store_into(end);

    program.add_argument("-stps", "--steps")
        .scan<'d', int>()
        .help("Number of steps in the q/theta range")
        .default_value(2000).store_into(nSteps);

    program.add_argument("-l", "--lambda")
        .scan<'g', double>()
        .help("Wavelength of the X-ray")
        .default_value(0.4).store_into(lambda);

    program.add_argument("-tt", "--twoThetaSpace")
        .help("Use two theta instead of q")
        .default_value(false)
        .implicit_value(true).store_into(twoTheta);

    bool benchmark;
    program.add_argument("--benchmark")
        .help("Run the benchmark")
        .default_value(false)
        .implicit_value(true).store_into(benchmark);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error &err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        exit(-1);
    }

    // check if either nRepeats > 0 or inputFilename is not empty
    if (nRepeats < 0 && inputFilename.empty())
    {
        std::cout << "Either nRepeats or inputFilename must be provided" << std::endl;
        std::cout << program;
        exit(-1);
    }

//    std::cout << "Shuffling positions: " << shufflePositions << std::endl;

    DebyeCalculator debye(nThreads, nCells, binsResolution, useMPI, useGPU, !dontUseLocalHistogram, 
    smallBins, pseudoCoal, fillGPU, !nonVerbose);
    debye.parallelHelper << "Shuffling positions: " << shufflePositions << "\n";
    debye.config.useGPUCellList = useGPUCellList;

    Positions positions;
    if (nRepeats > 0){
        positions = generateTestData(3.89070, nRepeats, stdDev, 10);
    } else {
        positions = readXYZ(inputFilename, " ", 2, typeMapping);
    }

    if (shufflePositions){
        std::cout << "Size: " << positions.size() << std::endl;
        positions = positions.subsample(1.0, 0);
        std::cout << "Size: " << positions.size() << std::endl;
    }


    auto results = debye.calculateProfile(positions, start, end, nSteps, twoTheta, lambda, "");

    if (!nonVerbose)
    {
        auto [pdf, profile] = std::prev(results.end())->second; // get last element
        debye.parallelHelper << profile.toString();
    }

    if (!outputDir.empty() && debye.parallelHelper.worldRank == 0){

        std::string prefix;
        if (!inputFilename.empty()){
            prefix = inputFilename.substr(0, inputFilename.find_last_of("."));
        }
        else
        {
            prefix = std::to_string(nRepeats) + "_repeats_" + std::to_string(stdDev) + "_stdDev_";
        }

        for (auto &[name, pdfProfile] : results)
        {   
            auto &pdf = pdfProfile.first;
            auto &profile = pdfProfile.second;
            std::string output_prefix = outputDir + "/" + prefix + "_" + name + "_";
            pdf.toCSV(output_prefix + "_pdf.csv");
            profile.toCSV(output_prefix + "_profile.csv");
        }
    }


    auto [pdf, profile] = std::prev(results.end())->second;
    debye.parallelHelper << "PDF calculation took: " << pdf.calculationTime << " s\n";
    debye.parallelHelper << "Intensity calculation took: " << profile.calculationTime << " s\n";

    return 0;
}