# ------------------------------------------------------------------------------#
# © 2025. Triad National Security, LLC. All rights reserved.  This program was
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

cmake_minimum_required(VERSION 3.24)

set(SPINER_VERSION 1.7.0)
project(
  spiner
  VERSION ${SPINER_VERSION}
  LANGUAGES CXX)

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
  SPINER_TEST_USE_KOKKOS_CUDA
  "Use kokkos cuda offloading for tests (affects submodule-build only)" ON
  "SPINER_TEST_USE_KOKKOS" ON)

# CTest
include(CTest)

# clang format
include(cmake/Format.cmake)

# Add a library
add_library(spiner INTERFACE)
add_library(spiner::spiner ALIAS spiner)
target_compile_features(spiner INTERFACE cxx_std_17)

# ##############################################################################
# Dependencies           #
# ##############################################################################

include(FetchContent)

# always use find_package() for FetchContent calls first NB: not an error if
# find_package() fails (will try to download if it can)
set(FETCHCONTENT_TRY_FIND_PACKAGE_MODE ALWAYS)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

if(NOT TARGET ports-of-call::ports-of-call)
  FetchContent_Declare(
    ports-of-call
    GIT_REPOSITORY https://github.com/lanl/ports-of-call
    GIT_TAG v2.0.0 # TODO: update if needed
  )
set(SPINER_DEP_LIST "ports-of-call")
endif()

# list of dependencies
set(SPINER_DEP_TARGETS "ports-of-call::ports-of-call")


if(SPINER_USE_HDF)
  if(NOT TARGET hdf5::hdf5)
    FetchContent_Declare(HDF5 FIND_PACKAGE_ARGS COMPONENTS C HL)
    list(APPEND SPINER_DEP_LIST "HDF5")
  endif()
  target_compile_definitions(spiner INTERFACE SPINER_USE_HDF)
  list(APPEND SPINER_DEP_TARGETS "hdf5::hdf5;hdf5::hdf5_hl")
endif()

# spiner_content_populate(NAMESPACE spinerDeps)
FetchContent_MakeAvailable(${SPINER_DEP_LIST})

# MPI::MPI_C dependency isn't propagated through hdf5::hdf5
# link to MPI::MPI_CXX explicitly
if(HDF5_IS_PARALLEL)
  find_package(MPI COMPONENTS CXX REQUIRED)
  target_link_libraries(spiner INTERFACE MPI::MPI_CXX)
endif()

target_link_libraries(spiner INTERFACE ${SPINER_DEP_TARGETS})

# Enables #include <spiner>
target_include_directories(
  spiner INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

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
