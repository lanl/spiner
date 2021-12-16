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

# Unit tests
if (NOT TARGET Catch2::Catch2)
  find_package(Catch2 QUIET)
endif()
if(NOT Catch2_FOUND)
  # Get catch2
  message(STATUS "Fetching Catch2 as needed")
  Include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.1)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
endif()

add_executable(test.bin test.cpp)
target_link_libraries(test.bin
  PRIVATE
  spiner::spiner
  Catch2::Catch2)

include(Catch)
catch_discover_tests(test.bin)

# Convergence test
include(FindPythonInterp)
add_executable(convergence.bin convergence.cpp)
target_link_libraries(convergence.bin
  PRIVATE spiner::spiner)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_convergence.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT convergence.dat convergence.png
  DEPENDS convergence.bin plot_convergence.py
  COMMAND ${CMAKE_COMMAND} -E time ./convergence.bin
  COMMAND ${PYTHON_EXECUTABLE} plot_convergence.py convergence.dat convergence.png)

add_custom_target(convergence
  DEPENDS convergence.png)
