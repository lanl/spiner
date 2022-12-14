#------------------------------------------------------------------------------#
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
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

cmake_minimum_required(VERSION 3.23)

project(spiner-cmake-test)

# Don't allow in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

find_package(spiner CONFIG REQUIRED)
if (NOT spiner_FOUND)
  message(FATAL_ERROR "Spiner not found.")
endif()

if (TARGET spiner::spiner)
  message(STATUS "The spiner library was succesfully provided")
  get_target_property(_sp_incdir spiner::spiner INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "ports of call include directory: ${_sp_incdir}")
else()
  message(FATAL_ERROR "The spiner library is missing!")
endif()

add_library(mytest libtest.cpp)
target_link_libraries(mytest
  PRIVATE
  spiner::spiner
  )
