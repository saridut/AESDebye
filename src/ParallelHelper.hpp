#pragma once
#include "DataTypes.hpp"
#include <streambuf>
#include <ostream>

#ifdef USE_MPI

#include <mpi.h>

#else
#define MPI_COMM_WORLD 0
#define MPI_Init(argc, argv)
#define MPI_Finalize()
#define MPI_Comm_size(comm, size)
#define MPI_Comm_rank(comm, rank)
#define MPI_Barrier(comm)
#define MPI_Wtime() 0.0
#define MPI_Datatype int
#define MPI_Op int
#define MPI_Type_contiguous(count, oldtype, newtype)
#define MPI_Type_commit(newtype)
#define MPI_Op_create(func, commute, op)
#define MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm)
#define MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
#define MPI_LONG_LONG_INT 1
#define MPI_UNSIGNED_LONG_LONG 1
#define MPI_Type_free(newtype)
#define MPI_Op_free(op)
#endif

#ifdef _OPENMP

#include <omp.h>

#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#define omp_set_num_threads(numThreads)
#endif

/**
 * @class ParallelHelper
 * @brief Helper class to deal with parallelism.
 * 
 * Provides an interface to deal with parallelism in the code.
 * It initializes MPI if required, sets up appropriate flags, and provides a consistent interface for printing to stdout.
 */
class ParallelHelper
{
public:
    int worldSize{1}; /**< Total size of the MPI world */
    int worldRank{0}; /**< Rank of the current process in the MPI world */
    int ompThreads{1}; /**< Number of OpenMP threads */
    int totalThreads{1}; /**< Total number of threads */
    bool useMPI{false}; /**< Flag to indicate if MPI is being used */
    bool verbose{true}; /**< Flag to indicate if verbose output is enabled */
    bool mpiInitialized{false}; /**< Flag to indicate if MPI has been initialized */
    double waitTime{0.0}; /**< Time spent waiting in MPI_Barrier */

    /**
     * @brief Construct a new ParallelHelper object
     * 
     * Default constructor.
     */
    ParallelHelper() = default;

    /**
     * @brief Construct a new ParallelHelper object
     * 
     * This constructor initializes MPI if required, sets up appropriate flags.
     * 
     * @param nThreads The number of threads to use.
     * @param useMPI Flag to indicate if MPI should be used.
     * @param verbose Flag to indicate if verbose output should be enabled.
     */
    ParallelHelper(int nThreads, bool useMPI, bool verbose) : useMPI(useMPI), verbose(verbose)
    {

#ifndef USE_MPI
            if (useMPI)
            {
                throw std::runtime_error("MPI support is not enabled in the current build. Please recompile with MPI support.");
            }
#endif

        // if MPI is being used, initialize it
        if (useMPI)
        {
            MPI_Init(NULL, NULL);
            MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
            MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
            mpiInitialized = true; // this is to ensure that MPI_Finalize is called
            *this << "MPI initialized: DO NOT initialize another calculator with MPI\n";
        }

        // set the number of OpenMP threads
        setNumThreads(nThreads > 0 ? nThreads : omp_get_max_threads());

        *this << "ParallelHelper: worldSize: " << worldSize << " worldRank: " << worldRank << " ompThreads: " << ompThreads << " total_threads: " << totalThreads << "\n";
    }

    /**
     * @brief Set the number of OpenMP threads.
     * 
     * @param numThreads The number of OpenMP threads to use.
     */
    void setNumThreads(int numThreads);

    /**
     * @brief Get the current thread number.
     * 
     * @return The current thread number in the global context.
     */
    int currentThread() const;

    /**
     * @brief Destroy the Parallel Helper object
     * 
     * It will finalize MPI if it was initialized.
     * 
     */
    ~ParallelHelper(){
        wait();
        if (useMPI && mpiInitialized) {
            *this << "Finalizing MPI\n";
            MPI_Finalize();
        }
    }


    /**
     * @brief Wait for all processes to reach this point.
     * 
     * This function will wait for all processes to reach this point only if MPI is being used.
     * 
     */
    void wait();

    /**
     * @brief Overloaded operator to print to stdout.
     * 
     * @tparam T Type of the object to print.
     * @param t Message to print.
     * @return ParallelHelper& 
     */
    template <typename T>
    ParallelHelper &operator<<(const T &t)
    {
        if ((worldRank == 0) && verbose)
        {
            std::cout << t;
        }
        return *this;
    }

    /**
     * @brief Print a section header to stdout.
     * 
     * Creates a nice section with blue color with text in the middle.
     * 
     * @param output The message to print.
     */
    void printSection(std::string output);
};


