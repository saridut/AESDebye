cmake_minimum_required(VERSION 3.15)
#cmake_policy(SET CMP0068 NEW)

if(SKBUILD)
    set(PROJECT_NAME ${SKBUILD_PROJECT_NAME})
    set(PROJECT_VERSION ${SKBUILD_PROJECT_VERSION})
else()
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" PYPROJECT_CONTENT)
        if(NOT DEFINED PROJECT_NAME)
            string(REGEX MATCH "name[ \t]*=[ \t]*\"([^\"]+)\"" _ "${PYPROJECT_CONTENT}")
            set(PROJECT_NAME ${CMAKE_MATCH_1})
        endif()
        if(NOT DEFINED PROJECT_VERSION)
            string(REGEX MATCH "version[ \t]*=[ \t]*\"([^\"]+)\"" _ "${PYPROJECT_CONTENT}")
            set(PROJECT_VERSION ${CMAKE_MATCH_1})
        endif()
    endif()

    if(NOT DEFINED PROJECT_NAME)
        set(PROJECT_NAME "aesdebye")
    endif()
    if(NOT DEFINED PROJECT_VERSION)
        set(PROJECT_VERSION "1.0.0")
    endif()
endif()

project(
    ${PROJECT_NAME}
    VERSION ${PROJECT_VERSION}
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
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

option(AUTOCONFIG "Enable automatic configuration" ON)
option(BUILD_LAMMPS_PLUGIN "Build LAMMPS plugin" OFF)

find_package(OpenMP REQUIRED)

if(AUTOCONFIG)
    message(STATUS "Autoconfig is ON - Detecting MPI, GPU and Python support")

    # Set MPI_HOME hint from Intel MPI or general MPI environment variables if available
    if(NOT DEFINED MPI_HOME)
        if(DEFINED ENV{I_MPI_ROOT})
            set(MPI_HOME "$ENV{I_MPI_ROOT}")
        elseif(DEFINED ENV{MPI_HOME})
            set(MPI_HOME "$ENV{MPI_HOME}")
        endif()
    endif()

    find_package(MPI QUIET)
    if(MPI_FOUND)
        set(ENABLE_MPI ON)
        message(STATUS "MPI detected and enabled")
    endif()

    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        set(ENABLE_GPU ON)
        message(STATUS "CUDA/GPU detected and enabled")
    endif()
endif()

include_directories(src extern)

# --- Centralized Settings (Flags, Defs, Libs) ---
add_library(common_settings INTERFACE)

target_compile_options(common_settings INTERFACE 
    -O3 -march=native -w -g -ffast-math
)

target_link_libraries(common_settings INTERFACE 
    stdc++
)
if(NOT APPLE)
    target_link_libraries(common_settings INTERFACE stdc++fs)
endif()

# --- MPI Setup ---
if(ENABLE_MPI)
    if(NOT DEFINED MPI_HOME)
        if(DEFINED ENV{I_MPI_ROOT})
            set(MPI_HOME "$ENV{I_MPI_ROOT}")
        elseif(DEFINED ENV{MPI_HOME})
            set(MPI_HOME "$ENV{MPI_HOME}")
        endif()
    endif()
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
add_library(aesdebye SHARED src/Calculator.cc src/CellList.cc
    src/IntensityKernels.cc src/ParallelHelper.cc
    src/PDF.cc src/PDFKernels.cc
    src/Positions.cc src/Profile.cc
    src/utils.cc)

target_link_libraries(aesdebye PUBLIC common_settings OpenMP::OpenMP_CXX)

# --- GPU Setup ---
if(ENABLE_GPU)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES) # Can be provided from outside, e.g. via -DCMAKE_CUDA_ARCHITECTURES=native or 80 90 whatever is supported
        if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0")
            set(CMAKE_CUDA_ARCHITECTURES 80) 
        elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "11.0")
            set(CMAKE_CUDA_ARCHITECTURES 80)
        elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
            set(CMAKE_CUDA_ARCHITECTURES 75)
        else()
            set(CMAKE_CUDA_ARCHITECTURES 60)
        endif()
        message(STATUS "CMAKE_CUDA_ARCHITECTURES set to latest supported: ${CMAKE_CUDA_ARCHITECTURES}")
    endif()
    
    add_library(aesdebyeGPU SHARED src/IntensityKernelsGPU.cu src/PDFKernelsGPU.cu)
    target_link_libraries(aesdebyeGPU PUBLIC common_settings)
    
    target_compile_definitions(common_settings INTERFACE USE_GPU)
    target_link_libraries(aesdebye PUBLIC aesdebyeGPU)
    message(STATUS "GPU support enabled")
endif()

# --- Python Setup ---
if(SKBUILD)
    option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
else()
    option(BUILD_PYTHON_BINDINGS "Build Python bindings" OFF)
endif()

if(BUILD_PYTHON_BINDINGS)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

    execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys; print(sys.prefix)" OUTPUT_VARIABLE Python_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(BINARY_DIR ${Python_PREFIX}/bin)

    execute_process(COMMAND ${Python_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])" OUTPUT_VARIABLE Python_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(FINAL_INSTALL_DIR ${Python_SITE_PACKAGES}/aesdebye)

    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_core src/export.cpp)
    target_link_libraries(_core PUBLIC aesdebye)

    # Python-specific definitions
    target_compile_definitions(_core PRIVATE 
        VERSION_INFO="${PROJECT_VERSION}"
        CXX_COMPILER="${CMAKE_CXX_COMPILER}"
    )
endif()

# --- CLI Setup ---
add_executable(aesdebye_cpp src/cli.cpp)
target_link_libraries(aesdebye_cpp PRIVATE aesdebye)

# add flags for CLI
target_compile_options(aesdebye_cpp PRIVATE -O3 -march=native -w -g -ffast-math)
# --- Installation ---
macro(install_symlink filepath sympath)
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})")
    install(CODE "message(\"-- Created symlink: ${sympath} -> ${filepath}\")")
endmacro()

if(BUILD_PYTHON_BINDINGS)
    install(TARGETS _core DESTINATION aesdebye)
    install(TARGETS aesdebye DESTINATION aesdebye)
    install(TARGETS aesdebye_cpp DESTINATION aesdebye)
    install(DIRECTORY src/ DESTINATION aesdebye/include FILES_MATCHING PATTERN "*.hpp")

    install_symlink(${FINAL_INSTALL_DIR}/aesdebye_cpp ${BINARY_DIR}/aesdebye_cpp)

    if(ENABLE_GPU)
        install(TARGETS aesdebyeGPU DESTINATION aesdebye)
    endif()
else()
    include(GNUInstallDirs)
    install(TARGETS aesdebye DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(TARGETS aesdebye_cpp DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(DIRECTORY src/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/aesdebye FILES_MATCHING PATTERN "*.hpp")

    if(ENABLE_GPU)
        install(TARGETS aesdebyeGPU DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
endif()

if(APPLE)
    set(RPATH "@loader_path")
else()
    set(RPATH "$ORIGIN")
endif()

#if(BUILD_PYTHON_BINDINGS)
#    set_target_properties(_core PROPERTIES INSTALL_RPATH ${RPATH})
#endif()
#set_target_properties(aesdebye PROPERTIES INSTALL_RPATH ${RPATH})
#set_target_properties(aesdebye_cpp PROPERTIES INSTALL_RPATH ${RPATH})
#
#if(ENABLE_GPU)
#    set_target_properties(aesdebyeGPU PROPERTIES INSTALL_RPATH ${RPATH})
#endif()

if (BUILD_LAMMPS_PLUGIN)
    add_subdirectory(lammps_plugin)
endif()
