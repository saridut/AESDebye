#include <PDF.hpp>


/**
 * @brief Operator overloading for adding two HistBin
 * 
 * @return HistBin - result of the addition
 */
HistBin HistBin::operator+(const HistBin &rhsBin) const {
    // create a new bin
    HistBin result = *this;

    // add the counts
    result.count += rhsBin.count;

    // check for overflow/underflow
    if ((result.delta > 0) && (rhsBin.delta > 0)) {
        if (LLONG_MAX - result.delta < rhsBin.delta) {
            result.delta -= LLONG_MAX;
            result.deltaOverflowCount++;
        }
    } else if ((result.delta < 0) && (rhsBin.delta < 0)) {
        if (-LLONG_MAX - result.delta > rhsBin.delta) {
            result.delta += LLONG_MAX;
            result.deltaOverflowCount--;
        }
    }

    // update the delta and overflow count
    result.delta += rhsBin.delta;
    result.deltaOverflowCount += rhsBin.deltaOverflowCount;
    return result;
}


void PDF::addSelfPairs(uint64 nPairs) {
    // return if the positions are not the same
    if (!samePositions)
        return;

    // double the count of each bin
#pragma omp parallel for schedule(static) default(none) shared(pdf_vector)
    for (uint binId = 0; binId < nBins; binId++) {
        pdf_vector[binId].count *= 2;
    }

    // add the self pairs to the zeroth bin
    pdf_vector[0].count += nPairs;
}

void
PDF::createPDFVectors() {
    // create vectors to store the bin centers, corrected centers, counts and deltas

    uncorrectedCenters = std::vector<double>(nBins);
    centers = std::vector<double>(nBins);
    counts = std::vector<double>(nBins);

    // calculate the bin centers, corrected centers, counts and deltas
#pragma omp parallel for schedule(static) default(none) shared(uncorrectedCenters, centers, counts)
    for (uint binId = 0; binId < nBins; binId++) {
        counts[binId] = (double) pdf_vector[binId].count;
        uncorrectedCenters[binId] = boxSize * (double) (binId) * centerMultiplicationFactor;
        centers[binId] = calculateCorrectedCenter(uncorrectedCenters[binId], (double) counts[binId],
                                                  calculateDelta(pdf_vector[binId]));
    }
}

bool PDF::test(uint64 positionsISize, uint64 positionsJSize) {
    // function to be called after the PDF has been finalized
    // self pairs are also added and double counted

    uint64 nPd = positionsISize * positionsJSize; // expected
    uint64 sum = 0; // total
    for (size_t binId = 0; binId < nBins; binId++) {
        uint64 count = pdf_vector[binId].count;
        sum += count;
    }
    testPassed = sum == nPd;
    return testPassed;
}


double PDF::calculateDelta(HistBin const &bin) const {
    double delta = (boxSize * boxSize * (0.0000000000000000010 *
                                         ((double) bin.deltaOverflowCount * (double) LLONG_MAX + (double) bin.delta)));
    if (!samePositions) delta = delta / 2.0;
    return delta;
}

double PDF::calculateCorrectedCenter(double center, double count, double delta) {
    double correctedCenter = center;
    if (fabs(delta) > FLT_EPSILON) {
        if (center > FLT_EPSILON) {
            double centerInv = 1.0 / center;
            double Delta = delta * centerInv /
                           count; // @@@ dCount has been already counted 2 times by the binning algorithm, whereas dDelta had been caunted only one time
            correctedCenter += Delta * (1.0 - 0.50 * Delta * centerInv * (1.0 - Delta * centerInv));
        } else {
            correctedCenter += sqrt(delta /
                                    count); // @@@ dCount has been already counted 2 times by the binning algorithm, whereas dDelta had been caunted only one time
        }
    }
    return correctedCenter;
}


PDF &PDF::operator+=(const PDF &other) {
    if (boxSize != other.boxSize) {
        throw std::runtime_error("Unit length mismatch");
    }
#pragma omp parallel for schedule(static) default(none) shared(pdf_vector, other)
    for (uint binId = 0; binId < nBins; binId++) {
        pdf_vector[binId] = pdf_vector[binId] + other.pdf_vector[binId];
    }
    return *this;
}


PDF PDF::operator+(const PDF &other) const {
    PDF result = *this;
    result += other;
    return result;
}

// Read write and other helpers

[[maybe_unused]] std::string PDF::toString() {
    uint64 sum = 0;
    uint64 filledBins = 0;
    int maxFilledBin = 0;
    for (uint64 binId = 0; binId < nBins; binId++) {
        uint64 count = pdf_vector[binId].count;
        sum += count;
        if (count > 0) {
            filledBins++;
            maxFilledBin = binId;
        }
    }
    std::string output = "Hist(filled=" + std::to_string(filledBins) + ", total counts="
                         + std::to_string(sum) + ", max filled=" + std::to_string(maxFilledBin) + " at " +
                         std::to_string(maxFilledBin * boxSize * centerMultiplicationFactor) + ")";
    return output;
}

// binary file dump for pdf_vector and histogram_gpu
void PDF::save(std::string const &filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    file.write(reinterpret_cast<char *>(&boxSize), sizeof(double));

    for (int i = 0; i < nBins; ++i) {
        file.write(reinterpret_cast<char *>(&pdf_vector[i]), sizeof(HistBin));
    }

    file.close();
}

double PDF::load(std::string const &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }
    file.read(reinterpret_cast<char *>(&boxSize), sizeof(double));

    for (int i = 0; i < nBins; ++i) {
        file.read(reinterpret_cast<char *>(&pdf_vector[i]), sizeof(HistBin));
    }

    file.close();
}

void PDF::toCSV(std::string const &filename, bool complete) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }
    file << "boxSize," << boxSize << ",nBins," << nBins << "\n";
    file << "binIdx,uncorrectedCenter,center,count,delta,deltaOverflowCount\n";
    for (int i = 0; i < nBins; ++i) {
        if ((!complete) && (counts[i] == 0))
            continue;

        file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10)
             << i << "," << uncorrectedCenters[i] << "," << centers[i] << "," << counts[i] << "," << pdf_vector[i].delta << "," << pdf_vector[i].deltaOverflowCount << "\n";
    }
}

void
PDF::readCSV(std::string const &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    std::string line;
    std::getline(file, line); // read the box size and nBins
    auto values = helpers::stringSplit(line, ",");
    if (values.size() != 4) {
        throw std::runtime_error("Error reading CSV file");
    }
    boxSize = std::stod(values[1]);
    nBins = std::stoi(values[3]);
    pdf_vector = std::vector<HistBin>(nBins);

    std::getline(file, line); // read the header and ignore

    while (std::getline(file, line)) {
        auto values = helpers::stringSplit(line, ",");
        if (values.size() != 6) {
            throw std::runtime_error("Error reading CSV file");
        }
        int binIdx = std::stoi(values[0]);
        uncorrectedCenters[binIdx] = std::stod(values[1]);
        centers[binIdx] = std::stod(values[2]);
        counts[binIdx] = std::stoull(values[3]);
        pdf_vector[binIdx].count = counts[binIdx];
        pdf_vector[binIdx].delta = std::stoll(values[4]);
        pdf_vector[binIdx].deltaOverflowCount = std::stoi(values[5]);
    }
}
