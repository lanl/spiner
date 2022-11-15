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
cmake_dependent_option (SPINER_TEST_USE_KOKKOS "Pull in Kokkos" OFF "BUILD_TESTING" OFF)
# only provide this option if `SPINER_USE_KOKKOS` is enabled
cmake_dependent_option (SPINER_TEST_USE_KOKKOS_CUDA "Use the Kokkos cuda backend" OFF "SPINER_TEST_USE_KOKKOS" OFF)

option (SPINER_FORCE_INTERNAL_PORTS "Force use of internal ports of call" OFF)

# MAUNEYC
# With dependent option change above, this is not
# possible.
#if (SPINER_USE_CUDA AND NOT SPINER_USE_KOKKOS)
# message(FATAL_ERROR "Currently cuda requires Kokkos")
#endif()

# CTEST
include(CTest)

# clang format
include(cmake/Format.cmake)

# Add a library
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

include(FetchContent)
# poc
FetchContent_Declare(
  ports-of-call
  GIT_REPOSITORY    https://github.com/lanl/ports-of-call
  GIT_TAG           20de024a3a389301a03d0e694a5cabfc6481f623
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(ports-of-call)
target_link_libraries(${PROJECT_NAME} INTERFACE ports-of-call::ports-of-call)

# HDF5
if(SPINER_USE_HDF)
  find_package(HDF5 COMPONENTS C HL)
  if (HDF5_FOUND)
    target_link_libraries(${PROJECT_NAME} INTERFACE HDF5::HDF5)
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -SPINER_USE_HDF=OFF")
  endif()
endif()

# Enables
# #include <spiner>
target_include_directories(${PROJECT_NAME}
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
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)
 
# ...and the version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${SPINER_VERSION}
  COMPATIBILITY SameMajorVersion
)
  
# Install the cmake configuration files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME} )
  
# Install header files
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spiner"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp"
)
  
# Install the export target. This will define the CMake target
# for external projects when used with `find_package`
install(EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${PROJECT_NAME}Targets.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
) 

# Export configuration for external projects that reference
# just our build-tree; e.g. for submodules. To use, ensure 
# `CMAKE_PREFIX_PATH` points to this source directory.
# NOTE: This config will not be relocatable!
export(TARGETS ${PROJECT_NAME} 
  NAMESPACE ${PROJECT_NAME}::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
)
