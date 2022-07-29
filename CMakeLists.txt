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

cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option (SPINER_USE_HDF
  "Pull in hdf5" OFF)
option (SPINER_USE_KOKKOS "Pull in Kokkos" OFF)
option (SPINER_USE_CUDA "Use the Kokkos cuda backend" OFF)
option (SPINER_FORCE_INTERNAL_PORTS "Force use of internal ports of call" OFF)

if (SPINER_USE_CUDA AND NOT SPINER_USE_KOKKOS)
  message(FATAL_ERROR "Currently cuda requires Kokkos")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

set(SPINER_VERSION 1.4.0)
project(spiner VERSION ${SPINER_VERSION})

# bring in some helpful CMake scripts 
# make cache variables for install destinations
include(GNUInstallDirs)
# package config file
include(CMakePackageConfigHelpers)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
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

# Disable cmake's package registry.
message(STATUS "CMake's package registry can cause problems with ports-of-call. "
  "It is therefore disabled by default. If you would like to re-enable it, "
  "you can use "
  "-DCMAKE_FIND_USE_PACKAGE_REGISTRY=ON and -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=ON")
set(CMAKE_FIND_USE_PACKAGE_REGISTRY OFF CACHE BOOL "" FORCE)
set(CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY OFF CACHE BOOL "" FORCE)

# CTEST
include(CTest)

# clang format
include(cmake/Format.cmake)

# Add a library
set(SPLIB "spiner")
add_library(${SPLIB} INTERFACE)
add_library(${SPLIB}::${SPLIB} ALIAS ${SPLIB})

# HDF5
if(SPINER_USE_HDF)
  if(SPINER_HDF5_INSTALL_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SPINER_HDF5_INSTALL_DIR}")
  endif()
  find_package(HDF5 COMPONENTS C HL)
  if (HDF5_FOUND)
    target_compile_definitions(${SPLIB} INTERFACE SPINER_USE_HDF)
    target_link_libraries(${SPLIB} INTERFACE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    target_include_directories(${SPLIB} INTERFACE ${HDF5_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "HDF5 was requested but not found. Can be disabled with -SPINER_USE_HDF=OFF")
  endif()
endif()

if (SPINER_USE_KOKKOS)
  set(PORTABILITY_STRATEGY "Kokkos" CACHE STRING "Portability strategy" FORCE)
  # Import Kokkos if not already available as a build target
  if (NOT TARGET Kokkos::kokkos)
    if (SPINER_USE_KOKKOS_SRC)
      set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON CACHE BOOL
        "Kokkos aggressive vectorization" FORCE)
      if (SPINER_USE_CUDA)
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL
          "Cuda backend for Kokkos" FORCE)
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL
          "Enable cuda lambda" FORCE)
        target_compile_options(
          ${SPLIB}
          INTERFACE # Generator expression shamelessly copied from EAP
          "$<$<COMPILE_LANGUAGE:CXX>:--expt-relaxed-constexpr;>"
          )
      endif()
      if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        set(Kokkos_ENABLE_DEBUG ON CACHE BOOL
          "Most general debug settings" FORCE)
        set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK ON CACHE BOOL
          "Bounds checking on Kokkos views" FORCE)
        set(Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK ON CACHE BOOL
          "Sanity checks on Kokkos DualView" FORCE)
      endif()
      if (EXISTS ${SPINER_USE_KOKKOS_SRC})
        message(STATUS "Kokkos source directory found. Adding source.")
        add_subdirectory(${SPINER_USE_KOKKOS_SRC} ${PROJECT_BINARY_DIR}/kokkos)
      else()
        message(STATUS "No Kokkos found. Getting it from web.")
        Include(FetchContent)
        FetchContent_Declare(
          Kokkos
          GIT_REPOSITORY https://github.com/kokkos/kokkos.git
          GIT_TAG        3.4.00)
        FetchContent_MakeAvailable(Kokkos)
      endif()
    else()
      message(STATUS "Looking for Kokkos via find package")
      if (SPINER_KOKKOS_INSTALL_DIR)
        list(APPEND CMAKE_PREFIX_PATH "${SPINER_KOKKOS_INSTALL_DIR}")
        set(Kokkos_ROOT "${SPINER_KOKKOS_INSTALL_DIR}/lib64/cmake/Kokkos")
      endif()
      find_package(Kokkos REQUIRED)
      if(NOT Kokkos_FOUND)
        message(FATAL_ERROR
          "Kokkos was requested but not found. You may have luck specifying the Kokkos install directory via -DSPINER_KOKKOS_INSTALL_DIR=/path/to/dir. Kokkos can be disabled with -DSPINER_USE_KOKKOS=OFF. You can also specify a source directory via -DSPINER_USE_KOKKOS_SRC=/path/to/kokkos or tell the build system to download Kokkos from the web via -DSPINER_USE_KOKKOS_SRC=ON.")
      endif()
    endif()
  endif()
  target_link_libraries(${SPLIB} INTERFACE Kokkos::kokkos)
endif ()

# ports of call
# TODO: Integrate with find-package
if (NOT TARGET ports-of-call::ports-of-call)
  message(STATUS "Looking for ports-of-call")
  if (SPINER_FORCE_INTERNAL_PORTS)
    message(STATUS "Using in-system ports-of-call")
    add_subdirectory(ports-of-call)
  else()
    find_package(ports-of-call CONFIG)
    if (NOT ports-of-call_FOUND)
      message(STATUS "Ports of call not available in-system. Using submodule.")
      add_subdirectory(ports-of-call)
    else()
      message(STATUS "Ports of call is available. Using in-system installation.")
    endif()
  endif()
else()
  message(STATUS "Ports of call already in environment")
endif()
target_link_libraries(${SPLIB} INTERFACE ports-of-call::ports-of-call)

# Enables
# #include <spiner>
target_include_directories(${SPLIB}
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

# Option imported from `CTest`
if (BUILD_TESTING)
  message("\nConfiguring tests")
  add_subdirectory(test)
endif()

# Coordinate external CMAKE export with targets
if (SPINER_USE_KOKKOS_SRC)
  message(WARNING
    "Installation not supported when Kokkos built from source. Disabling.")
else()
  install(TARGETS ${SPLIB}
    EXPORT ${SPLIB}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
  
  configure_package_config_file(${SPLIB}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${SPLIB}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${SPLIB}
    )
  
  # ...and the version file
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${SPLIB}ConfigVersion.cmake
    VERSION ${SPINER_VERSION}
    COMPATIBILITY SameMajorVersion )
  
  # Install the cmake configuration files
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SPLIB}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${SPLIB}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${SPLIB} )
  
  # Install header files
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spiner"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
    )
  
  # Install the export target. This will define the CMake target
  # for external projects when used with `find_package`
  install(EXPORT ${SPLIB}Targets
    NAMESPACE ${SPLIB}::
    FILE "${SPLIB}Targets.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${SPLIB}"
    COMPONENT dev)
  
  # Export configuration for external projects that reference
  # just our build-tree; e.g. for submodules. To use, ensure 
  # `CMAKE_PREFIX_PATH` points to this source directory.
  # NOTE: This config will not be relocatable!
  export(TARGETS ${SPLIB} NAMESPACE ${SPLIB}::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${SPLIB}Targets.cmake")
  # see https://cmake.org/cmake/help/latest/policy/CMP0090.html
  # since v3.15, this command is a no-op. prior, it populates
  # the user package directory.
  # `spiner` is commonly used "in-tree" downstream, so this 
  # *may* confuse downstream builds that still search the 
  # `~/.cmake` when doing `find_package()`
  #  export(PACKAGE ${SPLIB})
endif()
