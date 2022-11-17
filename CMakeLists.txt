#------------------------------------------------------------------------------#
# © 2021. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.24)

set(SPINER_VERSION 1.4.0)
project(spiner VERSION ${SPINER_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# bring in some helpful CMake scripts 
# make cache variables for install destinations
include(GNUInstallDirs)
# package config file
include(CMakePackageConfigHelpers)

# Don't allow in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# If the user doesn't specify a build type, prefer RelWithDebInfo
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(CMakeDependentOption)

# use HDF5
option (SPINER_USE_HDF "Pull in hdf5" OFF)
# use Kokkos offloading in tests
option (SPINER_TEST_USE_KOKKOS "Use kokkos offloading for tests" OFF)
option (SPINER_TEST_USE_KOKKOS_CUDA "Use kokkos cuda offloading for tests" OFF)

# These didn't behave like I expected.
# TODO: fix
# cmake_dependent_option (SPINER_TEST_USE_KOKKOS "Pull in Kokkos" OFF "BUILD_TESTING" OFF)
# cmake_dependent_option (SPINER_TEST_USE_KOKKOS_CUDA "Use the Kokkos cuda backend" OFF "SPINER_TEST_USE_KOKKOS" OFF)

# CTEST
include(CTest)

# clang format
include(cmake/Format.cmake)

# Add a library
add_library(spiner INTERFACE)
add_library(spiner::spiner ALIAS spiner)

include(FetchContent)
# i like chatty cmakes
set(FETCHCONTENT_QUIET FALSE)
# poc
FetchContent_Declare(
  ports-of-call
  GIT_REPOSITORY    https://github.com/lanl/ports-of-call
  # NOTE: This commit is still pending merge (PR16) 
  GIT_TAG           50bf032aa1522aef6da317cd07a7b5dd8c26aa82
  #  GIT_TAG           20de024a3a389301a03d0e694a5cabfc6481f623
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(ports-of-call)
target_link_libraries(spiner INTERFACE ports-of-call::ports-of-call)

# TODO: does FetchContent work here? not sure
# HDF5
if(SPINER_USE_HDF)
  find_package(HDF5 COMPONENTS C HL)
  if (HDF5_FOUND)
    target_link_libraries(spiner INTERFACE HDF5::HDF5)
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -SPINER_USE_HDF=OFF")
  endif()
endif()

# Enables
# #include <spiner>
target_include_directories(spiner
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
# Option imported from `CTest`
if (BUILD_TESTING)
  message("\nConfiguring tests")
  add_subdirectory(test)
endif()

install(
  TARGETS spiner
  EXPORT spinerTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/spinerConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/spinerConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spiner
)
 
# ...and the version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/spinerConfigVersion.cmake
  VERSION ${SPINER_VERSION}
  COMPATIBILITY SameMajorVersion
)
  
# Install the cmake configuration files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/spinerConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/spinerConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spiner )
  
# Install header files
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spiner"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/spiner
  FILES_MATCHING PATTERN "*.hpp"
)
  
# Install the export target. This will define the CMake target
# for external projects when used with `find_package`
install(EXPORT spinerTargets
  NAMESPACE spiner::
  FILE "spinerTargets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/spiner"
) 

# Export configuration for external projects that reference
# just our build-tree; e.g. for submodules. To use, ensure 
# `CMAKE_PREFIX_PATH` points to this source directory.
# NOTE: This config will not be relocatable!
export(TARGETS spiner 
  NAMESPACE spiner::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/spinerTargets.cmake"
)
