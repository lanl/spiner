# ------------------------------------------------------------------------------#
# © 2021. Triad National Security, LLC. All rights reserved.  This program was
# produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S.  Department of Energy/National Nuclear Security Administration.
# All rights in the program are reserved by Triad National Security, LLC, and
# the U.S. Department of Energy/National Nuclear Security Administration. The
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to
# reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
# ------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.19)

set(SPINER_VERSION 1.6.1)
project(spiner VERSION ${SPINER_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# bring in some helpful CMake scripts make cache variables for install
# destinations
include(GNUInstallDirs)
# package config file
include(CMakePackageConfigHelpers)
# dependent options
include(CMakeDependentOption)

# Don't allow in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
      "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles."
  )
endif()

# If the user doesn't specify a build type, prefer RelWithDebInfo
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# use HDF5
option(SPINER_USE_HDF "Use HDF5 for I/O" OFF)
option(SPINER_BUILD_TESTS "Compile Tests" OFF)
# use Kokkos offloading in tests
cmake_dependent_option(SPINER_TEST_USE_KOKKOS "Use kokkos offloading for tests"
                       ON "SPINER_BUILD_TESTS" OFF)
cmake_dependent_option(
  SPINER_TEST_USE_KOKKOS_CUDA "Use kokkos cuda offloading for tests" ON
  "SPINER_TEST_USE_KOKKOS" ON)

# CTest
include(CTest)

# clang format
include(cmake/Format.cmake)

# Add a library
add_library(spiner INTERFACE)
add_library(spiner::spiner ALIAS spiner)

# ##############################################################################
# Dependencies           #
# ##############################################################################

include(FetchContent)

# If we are on version 3.24+, then set FetchContent to always try to
# `find_package` before trying a download method
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)
  message(STATUS "FetchContent routines will try `find_package` first")
else()
  message(
    DEPRECATION
      "Detected cmake version (${CMAKE_VERSION}) older then 3.24. `spiner`"
      "will begin requiring this version soon.  The current depedency "
      "resolution is a soft-copy of the pattern of FetchContent introduced"
      "in cmake 3.24.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# content is a wrapper to `FetchContent` calls
include(content)
spiner_content_declare(
  ports-of-call
  NAMESPACE
  spinerDeps
  GIT_REPO
  https://github.com/lanl/ports-of-call
  # most recent relase as of April 05, 2023
  GIT_TAG
  v1.5.1)

if(SPINER_USE_HDF)
  spiner_content_declare(
    HDF5
    NO_FETCH
    NAMESPACE
    spinerDeps
    COMPONENTS
    C
    HL
    EXPECTED_TARGETS
    hdf5::hdf5
    hdf5::hdf5_hl)
  target_compile_definitions(spiner INTERFACE SPINER_USE_HDF)
endif()

spiner_content_populate(NAMESPACE spinerDeps)

# We don't know about this until HDF5 is populated, so we need to delay until
# it's been (NB: newer versions of HDF5/CMake appear to fix this)
if(HDF5_IS_PARALLEL)
  spiner_content_declare(
    MPI
    NO_FETCH
    COMPONENTS
    CXX
    EXPECTED_TARGETS
    MPI::MPI_CXX
    NAMESPACE
    spinerMPI)
  spiner_content_populate(NAMESPACE spinerMPI)
endif()

target_link_libraries(spiner INTERFACE ${spinerDeps_POPULATED_TARGETS}
                                       ${spinerMPI_POPULATED_TARGETS})

# Enables #include <spiner>
target_include_directories(
  spiner INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# Option imported from `CTest`
if(SPINER_BUILD_TESTS)
  message(STATUS "\nConfiguring tests")
  add_subdirectory(test)
endif()

# ----------------------------------------------------------------------------#
# Generate config settings
# ----------------------------------------------------------------------------#

# get all available `SPINER_` cmake cache variables set during configuration
get_cmake_property(_variableNames VARIABLES)
string(REGEX MATCHALL "(^|;)SPINER_[A-Za-z0-9_]*" _matchedVars
             "${_variableNames}")

# use config template to generate the configuration of the build not sure why
# CMake doesn't do this automatically, but w/e
foreach(_variableName ${_matchedVars})
  set(SPINER_CONFIG_CODE
      "${SPINER_CONFIG_CODE}\nset(${_variableName} \"${${_variableName}}\")")
endforeach()

install(
  TARGETS spiner
  EXPORT spinerTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/spinerConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/spinerConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spiner)

# ...and the version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/spinerConfigVersion.cmake
  VERSION ${SPINER_VERSION}
  COMPATIBILITY SameMajorVersion)

# Install the cmake configuration files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/spinerConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/spinerConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spiner)

# Install header files
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spiner"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.hpp")

# Install the export target. This will define the CMake target for external
# projects when used with `find_package`
install(
  EXPORT spinerTargets
  NAMESPACE spiner::
  FILE "spinerTargets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/spiner")

# Export configuration for external projects that reference just our build-tree;
# e.g. for submodules. To use, ensure `CMAKE_PREFIX_PATH` points to this source
# directory. NOTE: This config will not be relocatable!
export(
  TARGETS spiner
  NAMESPACE spiner::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/spinerTargets.cmake")
