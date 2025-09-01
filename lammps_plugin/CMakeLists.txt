cmake_minimum_required(VERSION 3.20)
project(ComputeHelloPlugin LANGUAGES CXX)

# cmake standard boilerplate
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(LAMMPS REQUIRED CONFIG)

add_library(compute_hello_plugin SHARED compute_hello.cpp)
target_link_libraries(compute_hello_plugin PRIVATE LAMMPS::lammps)

target_include_directories(compute_hello_plugin PRIVATE
    ~/temp/lammps/src/
    ~/.local/include/lammps/external/json  # <--- add this
)

# target_include_directories(compute_hello_plugin PRIVATE ${LAMMPS_INCLUDE_DIRS})

set_target_properties(compute_hello_plugin PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

