/**
 * @file Positions.hpp
 * @brief Contains the declaration of the Positions class.
 */
#pragma once
#include <DataTypes.hpp>

/**
 * @brief The Positions class represents the positions of atoms in a system.
 * 
 * This class stores the atom IDs, atom types, and coordinates (X, Y, Z) of atoms in a system.
 * 
 */
class Positions
{
private:
 
   std::vector<int64> X; //< X coordinates (int64)
   std::vector<int64> Y; //< Y coordinates (int64)
   std::vector<int64> Z; //< Z coordinates (int64)

   std::vector<double> Xd; //< X coordinates (double)
   std::vector<double> Yd; //< Y coordinates (double)
   std::vector<double> Zd; //< Z coordinates (double)
   double boxSize{1.0}; //< Size of the box


public:
   std::string element = "None"; //< Element of the atoms
   std::vector<int> cellHeads; //< Cell heads. Index of the first atom in each cell
   std::vector<int> atomCounts; //< Atom counts. Number of atoms in each cell
   std::vector<std::string> chemicalSymbols;//< Atom elements
   std::vector<std::string> selectionIds; //*< Atom IDs */
   /**
    * @brief Construct a new Positions object
    *
    * @param boxSize The unit length of the box. Default is 1.0
    */
   explicit Positions(double boxSize = 1.0) : boxSize(boxSize) {}

   /**
    * @brief Construct a new Positions object
    * 
    * @param positions (N, 3) A vector of  array of positions of atoms
    * @param boxMin
    * @param boxMax (3, ) array of maximum coordinates of the box (X, Y, Z)
    */
   Positions(std::vector<std::string> chemicalSymbols, std::vector<std::array<double, 3>> positions,
                std::vector<std::string> selectionIds={}, double boxMin=1e100, double boxMax=-1e100);
   /**
    * @brief Scale up the positions to int64
    */
   void scaleUp();
   
   /**
    * @brief Resize the positions containers
    * 
    * @param nAtoms new size of the positions
    */
   void resize(int nAtoms);

   /**
    * @brief Sample a subset of the positions
    *
    * @param subsampleRatio ratio of particles to choose from the total (0, 1]
    * @param seed seed for the random number generator
    * @return Sub-sampled positions
    */
   Positions subsample(double subsampleRatio, int seed=-1);

   /**
    * @brief Get the positions as a std::vector of std::array<double, 3>
    *
    * @return std::vector<std::array<double, 3>>
    */
   std::vector<std::array<double, 3>> toStdVector();

   /**
    * @brief Write the positions to a file in XYZ format
    * 
    * @param filename Name of the file
    */
   void toXYZ(std::string filename);


   [[nodiscard]] std::vector<std::string> getUniqueElements() const;
   [[nodiscard]] std::vector<std::string> getUniqueSelectionIds() const;

   /**
    * filter the positions based on the atom selectionIds, chemicalSymbols
    * 
    */
   Positions filterBySelectionId(std::string atomId);

   /**
    * filter the positions based on the atom selectionIds, chemicalSymbols
    * 
    */
   Positions filterByElement(std::string atomType);

   /**
    * @brief Get the size of the positions
    * 
    * @return int64 
    */
   [[nodiscard]] size_t size() const
   {
      return X.size();
   }

   /**
    * @brief Push back from another Positions object
    * 
    */
   void push_back(const Positions &positions, uint index);

   /**
    * @brief new positions from given positions and index list
    */
   static
   Positions fromIndexList(const Positions &positions, const std::vector<int> &indexList);

   /**
    * sliceCellList method
    */
   void sliceCellList(int start, int stop);


   size_t getCellIndex(size_t particleIdx) const;
   /**
    * @brief Get box size
    * 
    * @return double box size
    */
   double getBoxSize() const { return boxSize; }

   /**
    * @brief Set box size
    * 
    * @param size box size
    */
   void setBoxSize(double size)
   {
      // set the box size
      boxSize = size; 

      // scaling is bad now: scale up again!
      scaleUp();
   }


   // Friend declarations
   // for updating the cell heads and atom counts
   friend class CellList; 
   friend class PositionsGPU;
   friend Positions makePeriodic(const Positions &positions, int nx, int ny, int nz,
                      double offset_x, double offset_y, double offset_z);
   friend Positions readXYZ(const std::string& filename, std::string delimiter,
                  uint skipLines, const std::string& columns);

   // for calculating the actual bin index - a kernel function
   friend inline void findBinIndex(Positions const &positions1, 
                                   Positions const &positions2, 
                                   uint ii, uint jj, 
                                   int64 &binId, int64 &delta, float scaleDown,
                                   int64 scaleUp);
};
