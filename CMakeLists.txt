cmake_minimum_required(VERSION 3.5)

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

option(AUTOCONFIG "Enable automatic configuration" OFF)
option(BUILD_LAMMPS_PLUGIN "Build LAMMPS plugin" OFF)

find_package(OpenMP REQUIRED)

if(AUTOCONFIG)
    message(STATUS "Autoconfig is ON - Detecting MPI, GPU and Python support")

    find_package(MPI QUIET)
    if(MPI_FOUND)
        set(ENABLE_MPI ON)
        message(STATUS "MPI detected and enabled")
    endif()

    find_package(CUDA QUIET)
    if(CUDA_FOUND)
        set(ENABLE_GPU ON)
        message(STATUS "CUDA/GPU detected and enabled")

        include(FindCUDA/select_compute_arch)
        CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS)
        string(REPLACE " " ";" INSTALLED_GPU_CCS "${INSTALLED_GPU_CCS}")
        string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS}")
        message(STATUS "CUDA architecture: ${CUDA_ARCH_LIST}")
    endif()
endif()

include_directories(include extern)

# --- Centralized Settings (Flags, Defs, Libs) ---
add_library(common_settings INTERFACE)

target_compile_options(common_settings INTERFACE 
    -O3 -march=native -w -g -ffast-math
)
target_compile_definitions(common_settings INTERFACE 
    _GLIBCXX_USE_CXX11_ABI=0
)
target_link_libraries(common_settings INTERFACE 
    stdc++
)
if(NOT APPLE)
    target_link_libraries(common_settings INTERFACE stdc++fs)
endif()

# --- MPI Setup ---
if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        target_include_directories(common_settings INTERFACE ${MPI_INCLUDE_PATH})
        target_link_libraries(common_settings INTERFACE MPI::MPI_CXX)
        target_compile_definitions(common_settings INTERFACE USE_MPI)
        message(STATUS "MPI support enabled")
    else()
        set(ENABLE_MPI OFF)
    endif()
endif()

# --- Core Library ---
add_library(debyeLib SHARED src/Calculator.cc src/CellList.cc
    src/IntensityKernels.cc src/ParallelHelper.cc
    src/PDF.cc src/PDFKernels.cc
    src/Positions.cc src/Profile.cc
    src/utils.cc)

target_link_libraries(debyeLib PUBLIC common_settings OpenMP::OpenMP_CXX)

# --- GPU Setup ---
if(ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    
    add_library(debyeLibGPU SHARED src/IntensityKernelsGPU.cu src/PDFKernelsGPU.cu)
    target_link_libraries(debyeLibGPU PUBLIC common_settings ${CUDA_LIBRARIES})
    
    target_compile_definitions(common_settings INTERFACE USE_GPU)
    target_link_libraries(debyeLib PUBLIC debyeLibGPU)
    message(STATUS "GPU support enabled")
endif()

# --- Python Setup ---
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys; print(sys.prefix)" OUTPUT_VARIABLE Python_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
set(BINARY_DIR ${Python_PREFIX}/bin)

execute_process(COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE Python_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
set(FINAL_INSTALL_DIR ${Python_SITE_PACKAGES}/aesdebye)

find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(_core src/export.cpp)
target_link_libraries(_core PUBLIC debyeLib)

# Python-specific definitions
target_compile_definitions(_core PRIVATE 
    VERSION_INFO="${PROJECT_VERSION}"
    CXX_COMPILER="${CMAKE_CXX_COMPILER}"
)

# --- CLI Setup ---
add_executable(debye_cpp src/cli.cpp)
target_link_libraries(debye_cpp PRIVATE debyeLib)

# --- Installation ---
macro(install_symlink filepath sympath)
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})")
    install(CODE "message(\"-- Created symlink: ${sympath} -> ${filepath}\")")
endmacro()

install(TARGETS _core DESTINATION aesdebye)
install(TARGETS debyeLib DESTINATION aesdebye)
install(TARGETS debye_cpp DESTINATION aesdebye)
install(DIRECTORY include/ DESTINATION aesdebye/include)

install_symlink(${FINAL_INSTALL_DIR}/debye_cpp ${BINARY_DIR}/debye_cpp)

if(ENABLE_GPU)
    install(TARGETS debyeLibGPU DESTINATION aesdebye)
endif()

if(APPLE)
    set(RPATH "@loader_path")
else()
    set(RPATH "$ORIGIN")
endif()

set_target_properties(_core PROPERTIES INSTALL_RPATH ${RPATH})
set_target_properties(debyeLib PROPERTIES INSTALL_RPATH ${RPATH})
set_target_properties(debye_cpp PROPERTIES INSTALL_RPATH ${RPATH})

if(ENABLE_GPU)
    set_target_properties(debyeLibGPU PROPERTIES INSTALL_RPATH ${RPATH})
endif()

if (BUILD_LAMMPS_PLUGIN)
    add_subdirectory(lammps_plugin)
endif()
