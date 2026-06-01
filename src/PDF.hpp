/**
 * This file contains the PDF class. This is the data structure that holds the PDF.
 */

#pragma once
#include <DataTypes.hpp>


/**
 * @brief Histogram Bin struct - Mainly used for CPU
 */
struct HistBin // 64 * 3
{
    uint64 count{0}; //< Counts
    int64 delta{0};    //< The delta value.
    int deltaOverflowCount{0}; //< The delta overflow count value.
    int tillOverFlowReset{MAX_DELTA_UPDATE_COUNT}; // keep it int and not uint, the check is for <= 0 //< when to reset the overflow count
    HistBin operator+(const HistBin &rhsBin) const;
};
#define HIST_BIN_SIZE 3

/**
 * @brief Pair distribution function data structure. This holds the PDF data.
 */
class PDF
{
private:
   // these vector data structures are required for efficiency reasons!
   // All the actual calculations are done using these data structures
   // Finally when the data is to be used, it can be converted into vectors using createPDFVectors() function
   std::vector<HistBin> pdf_vector; // CPU pdf vector
   uint nBins{N_BINS};              // Number of bins in the PDF
   double centerMultiplicationFactor{0.0000010};

   // Box size, this is used to get the correct centers of the bin
   // two PDFs can only if added if they have the same box size
   double boxSize{1.0};


   /**
    * @brief Calculate a single delta value for a given bin
    *
    * @param bin bin for which the delta value is to be calculated
    * @return double delta value
    */
   [[nodiscard]] double calculateDelta(HistBin const &bin) const;

   /**
    * @brief Calculate the corrected center of a bin
    *
    * This applies the formula given in the paper to calculate the corrected center
    *
    * @param center center of the bin
    * @param count count of the bin
    * @param delta delta of the bin
    * @return double corrected center
    */
   static double calculateCorrectedCenter(double center, double count, double delta);

   // empty constructor only available to Calculator class

    /**
  * @brief Get the actual PDF data
  *
  * This will create the uncorrectedCenters, centers, counts
  */
    void
    createPDFVectors();

   std::vector<double> uncorrectedCenters, centers;
   std::vector<double> counts;

    /**
  * @brief Double the count of the PDF and adds self pairs
  *
  * This function checks if the positions are the same and then doubles the count of the PDF
  */
    void addSelfPairs(uint64 nPairs);

public:
   // samePositions flag, this is used to divide to update the delta values
   bool samePositions{true};
   std::string elementI = "";
   std::string elementJ = "";
   double calculationTime = 0.0;
   bool testPassed = false;

   // default constructor, required for creating maps.
   PDF() {};

   /**
    * @brief Construct a new PDF object
    *
    * @param boxSize box size of the PDF
    * @param samePositions flag to indicate if the positions are the same, default is true
    */
   explicit PDF(double boxSize, bool samePositions = true, uint nBins = N_BINS,
               std::string elementI="", std::string elementJ="") : boxSize(boxSize),
                                                            samePositions(samePositions),
                                                            nBins(nBins),
                                                            pdf_vector(nBins),
                                                            elementI(elementI),
                                                            elementJ(elementJ)
   {
      if (nBins == 173210)
      {
         centerMultiplicationFactor = 0.000010;
      }
   }

   /**
    * @brief Helper to print the PDF object
    *
    * This binds to the __repr__ function in Python - if the python binding is used
    *
    * @return std::string
    */
   std::string toString();

   uint getNBins() const { return nBins; }
   [[nodiscard]] const std::vector<double>& getUncorrectedCenters() const { return uncorrectedCenters; }
   [[nodiscard]] const std::vector<double>& getCenters() const { return centers; }
    [[nodiscard]] const std::vector<double>& getCounts() const { return counts; }

    /**
     * @brief Get the box size of the PDF
     *
     * @return double box size
     */
    [[nodiscard]] double getBoxSize() const { return boxSize; }

   bool test(uint64 positionsISize, uint64 positionsJSize);

   std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
   getNormalizedPDF();

   std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
   getRDF();
   /**
    * @brief Add another PDF to this PDF
    *
    * This function adds another PDF to this PDF. The two PDFs must have the same box size.
    *
    * @param other PDF to be added
    * @return PDF& reference to this PDF
    */
   PDF &operator+=(const PDF &other);
   PDF operator+(const PDF &other) const;

   /**
    * @brief save the PDF to a binary file
    *
    * @param filename filename to save the PDF to
    */
   void save(std::string const &filename);

   /**
    * @brief load the PDF from a binary file
    *
    * @param filename filename to load the PDF from
    * @return double box size of the PDF
    */
   double load(std::string const &filename);

   /**
    * @brief save the PDF to a CSV file
    *
    * @param filename filename to save the PDF to
    * @param complete flag to indicate if the complete PDF should be saved, otherwise skip empty bins
    */
   void toCSV(std::string const &filename, bool complete = false);

   /**
    * @brief Read a CSV file and return the data
    *
    * @param filename filename to read the CSV file from
    */
   void readCSV(std::string const &filename);




   // Debugging functions

   // friend declarations
   friend class DebyeCalculator;
};

///**
// * @brief Calculate a single delta value for a given bin
// *
// * @param bin bin for which the delta value is to be calculated
// * @return double delta value
// */
//[[nodiscard]] double calculateDelta(UBin const &bin) const;