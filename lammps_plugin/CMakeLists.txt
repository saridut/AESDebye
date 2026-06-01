find_package(LAMMPS REQUIRED CONFIG)
message("lammps include dir is: ${LAMMPS_INCLUDE_DIRS}")

# if not given as an argument, then set it empty
if (NOT DEFINED LAMMPS_SRC_DIR)
    set(LAMMPS_SRC_DIR "")
    message(FATAL_ERROR "LAMMPS_SRC_DIR is not defined. Please set it to the LAMMPS source directory. This is required to get the headers!")
endif()

add_library(lammps_debye_plugin SHARED compute_debye.cpp)
target_link_libraries(lammps_debye_plugin PRIVATE LAMMPS::lammps)
target_link_libraries(lammps_debye_plugin PRIVATE aesdebye)
target_include_directories(lammps_debye_plugin PRIVATE ${LAMMPS_SRC_DIR})

# installation dir, need to use full path to load. or just use the build dir to load!
install(TARGETS lammps_debye_plugin DESTINATION aesdebye)
set_target_properties(lammps_debye_plugin PROPERTIES INSTALL_RPATH ${RPATH})

