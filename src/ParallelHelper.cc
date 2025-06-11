#include <ParallelHelper.hpp>

void ParallelHelper::setNumThreads(int numThreads)
{   
    // set number of OpenMP threads
    omp_set_num_threads(numThreads);
    
    // make sure the number of threads is set correctly
    assert(omp_get_max_threads() == numThreads);

    // update variables
    ompThreads = numThreads;
    totalThreads = worldSize * ompThreads;
}

int ParallelHelper::currentThread() const
{
    // get global thread number
    return omp_get_thread_num() + ompThreads * worldRank;
}

void ParallelHelper::wait(){
    // wait if MPI is being used and it is initialized with this object
    if (useMPI && mpiInitialized) {
        double start = helpers::get_wall_time();
        MPI_Barrier(MPI_COMM_WORLD);
        waitTime += helpers::get_wall_time() - start;
    }
}

void ParallelHelper::printSection(std::string output){
    const int totalLength = 80; // Total length of the output line
    const int outputLength = output.length(); // Length of the output string
    const int sideLength = (totalLength - outputLength - 4) / 2; // Calculate the length of each side of dashes
    std::string sideDashes(sideLength, '-'); // Create a string with 'sideLength' dashes

    *this << "\033[1;34m"; // Set the color to blue
    *this << sideDashes << " " << output << " " << sideDashes;
    if ((outputLength % 2) != (totalLength % 2)) { // Check if there's a mismatch in the total length
        *this << "-"; // Add an extra dash if needed
    }
    *this << "\033[0m\n"; // Reset the color
}


