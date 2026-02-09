cmake_minimum_required(VERSION 3.5)

# compilation tested with gcc12, intelmpi, cuda12.6 on fritz and alex

# if name and version are not set, set them, this happens when using cmake directly, without pip
if(NOT DEFINED PROJECT_NAME)
    set(PROJECT_NAME "aesdebye")
endif()

if(NOT DEFINED PROJECT_VERSION)
    set(PROJECT_VERSION "0.0.1")
endif()

project(
    ${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# automatically detect the compiler available, prefference is given to intel, clang and then gcc
if(NOT CMAKE_CXX_COMPILER)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_CXX_COMPILER "icpc")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_COMPILER "clang++")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_COMPILER "g++")
    else()
        message(FATAL_ERROR "No C++ compiler found")
    endif()
endif()

message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")

# Options
option(AUTOCONFIG "Enable automatic configuration" ON)
option(BUILD_LAMMPS_PLUGIN "Build LAMMPS plugin" OFF)

# OpenMP - always required for basic single node parallelism
find_package(OpenMP REQUIRED)

# Automatic Configuration
if(AUTOCONFIG)
    message(STATUS "Autoconfig is ON - Detecting MPI, GPU and Python support")

    # Detect MPI
    find_package(MPI QUIET)

    if(MPI_FOUND)
        set(ENABLE_MPI ON)
        message(STATUS "MPI detected and enabled")
    else()
        message(STATUS "MPI not detected")
    endif()

    # Detect CUDA/GPU
    find_package(CUDA QUIET)

    if(CUDA_FOUND)
        set(ENABLE_GPU ON)
        message(STATUS "CUDA/GPU detected and enabled")

        # Detect CUDA architecture
        include(FindCUDA/select_compute_arch)
        CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS)
        string(REPLACE " " ";" INSTALLED_GPU_CCS "${INSTALLED_GPU_CCS}")
        string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS}")
        # set(CMAKE_CUDA_ARCHITECTURES "native")
        message(STATUS "CUDA architecture: ${CUDA_ARCH_LIST}")
    else()
        message(STATUS "CUDA/GPU not detected")
    endif()
endif()

include_directories(include extern)

# Library setup
add_library(debyeLib SHARED src/Calculator.cc src/CellList.cc
    src/IntensityKernels.cc src/ParallelHelper.cc
    src/PDF.cc src/PDFKernels.cc
    src/Positions.cc src/Profile.cc
    src/utils.cc)

target_compile_options(debyeLib PUBLIC -O3 -march=native -w -g -ffast-math)
target_link_libraries(debyeLib PUBLIC OpenMP::OpenMP_CXX)

if(NOT APPLE)
    target_link_libraries(debyeLib PUBLIC stdc++fs)
endif()

# MPI Setup
if(ENABLE_MPI)
    find_package(MPI)

    if(MPI_FOUND)
        target_include_directories(debyeLib PUBLIC ${MPI_INCLUDE_PATH})
        target_link_libraries(debyeLib PUBLIC MPI::MPI_CXX)
        target_compile_definitions(debyeLib PUBLIC -DUSE_MPI)
        message(STATUS "MPI support enabled, Inc dir: ${MPI_INCLUDE_PATH}")
    else()
        message(STATUS "MPI not found, skipping MPI support")
        set(ENABLE_MPI OFF)
    endif()
endif()

# GPU Setup
if(ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    add_library(debyeLibGPU SHARED src/IntensityKernelsGPU.cu src/PDFKernelsGPU.cu)
    target_compile_options(debyeLibGPU PRIVATE -O3 -march=native)
    target_link_libraries(debyeLibGPU PRIVATE ${CUDA_LIBRARIES})
    target_compile_definitions(debyeLib PUBLIC -DUSE_GPU)
    target_compile_definitions(debyeLibGPU PUBLIC -DUSE_GPU)
    target_link_libraries(debyeLib PUBLIC debyeLibGPU)
    message(STATUS "GPU support enabled")
endif()

# Python Library Setup
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sys; print(sys.prefix)"
    OUTPUT_VARIABLE Python_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(BINARY_DIR ${Python_PREFIX}/bin)

# get site-packages directory
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE Python_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(FINAL_INSTALL_DIR ${Python_SITE_PACKAGES}/aesdebye)

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_core src/export.cpp)
target_compile_options(_core PRIVATE -O3 -march=native -w)
target_link_libraries(_core PUBLIC debyeLib)

# version and compiler info for the python module
target_compile_definitions(_core PUBLIC -DVERSION_INFO="${PROJECT_VERSION}")
target_compile_definitions(_core PUBLIC -DCXX_COMPILER="${CMAKE_CXX_COMPILER}")

if(ENABLE_GPU)
    target_compile_definitions(_core PUBLIC -DUSE_GPU)
endif()

if(ENABLE_MPI)
    target_compile_definitions(_core PUBLIC -DUSE_MPI)
endif()

message(STATUS "Building Python library")

# # CLI Setup
add_executable(debye_cpp src/cli.cpp)
target_link_libraries(debye_cpp PRIVATE debyeLib)

# Install
macro(install_symlink filepath sympath)
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})")
    install(CODE "message(\"-- Created symlink: ${sympath} -> ${filepath}\")")
endmacro(install_symlink)

install(TARGETS _core DESTINATION aesdebye)
install(TARGETS debyeLib DESTINATION aesdebye)
install(TARGETS debye_cpp DESTINATION aesdebye)
install(DIRECTORY include/ DESTINATION aesdebye/include)

# symlink binary into the bin directory, so it can be called from the command line
install_symlink(${FINAL_INSTALL_DIR}/debye_cpp ${BINARY_DIR}/debye_cpp)

if(ENABLE_GPU)
    install(TARGETS debyeLibGPU DESTINATION aesdebye)
endif()

# based on linux or mac set the rpath to $ORIGIN or @loader_path
if(APPLE)
    # ("Setting RPATH to @loader_path")
    message("Setting RPATH to @loader_path")
    set(RPATH "@loader_path")
else()
    # ("Setting RPATH to $ORIGIN")
    message("Setting RPATH to $ORIGIN")
    set(RPATH "$ORIGIN")
endif()


# rpath using $ORIGIN
set_target_properties(_core PROPERTIES INSTALL_RPATH ${RPATH})
set_target_properties(debyeLib PROPERTIES INSTALL_RPATH ${RPATH})
set_target_properties(debye_cpp PROPERTIES INSTALL_RPATH ${RPATH})

if(ENABLE_GPU)
    set_target_properties(debyeLibGPU PROPERTIES INSTALL_RPATH ${RPATH})
endif()

if (BUILD_LAMMPS_PLUGIN)
    add_subdirectory(lammps_plugin)
endif()
