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

spiner_content_declare(Catch2
  GIT_REPO       https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1 # or any later release
  EXPECTED_TARGETS Catch2::Catch2
  NAMESPACE spinerTest
)

add_library(spiner_test INTERFACE)
# instructions for local content
# https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_SOURCE_DIR_%3CuppercaseName%3E
# TLDR: set FETCHCONTENT_SOURCE_DIR_KOKKOS=path/to/local/checkout
if (SPINER_TEST_USE_KOKKOS)
  # put together a list of kokkos configure options
  list(APPEND _spiner_content_opts "Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION")
  if (SPINER_TEST_USE_KOKKOS_CUDA)
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA")
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA_LAMBDA")
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA_CONSTEXPR")
  endif()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_DEBUG")
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_DEBUG_BOUNDS_CHECK")
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK")
  endif()
 
  spiner_content_declare(
    Kokkos
    GIT_REPO https://github.com/kokkos/kokkos.git
    GIT_TAG        61d7db55fceac3318c987a291f77b844fd94c165
    EXPECTED_TARGETS Kokkos::kokkos
    ENABLE_OPTS ${_spiner_content_opts}
    NAMESPACE spinerTest
  )
 
  target_compile_definitions(spiner_test INTERFACE PORTABILITY_STRATEGY_KOKKOS)
endif()

# populate test content
# places targets in `spinerTest_POPULATED_TARGETS`
spiner_content_populate(
  NAMESPACE spinerTest
)
target_link_libraries(spiner_test
  INTERFACE
    ${spinerTest_POPULATED_TARGETS}
)
# add unit tests
add_executable(test.bin test.cpp)
target_link_libraries(test.bin
  PRIVATE
    spiner::spiner
    spiner_test
)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
include(Catch)
catch_discover_tests(test.bin)

# Benchmarking
add_executable(benchmark.bin benchmark.cpp)
target_link_libraries(benchmark.bin
  PRIVATE 
    spiner::spiner
    spiner_test
)

# Convergence test
add_executable(convergence.bin convergence.cpp)
target_link_libraries(convergence.bin
  PRIVATE 
    spiner::spiner
    spiner_test
)

find_package(Python3 COMPONENTS Interpreter)
if (Python3_Interpreter_FOUND)
  message(STATUS "Checking for numpy, matplotlib.")
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, matplotlib"
    RESULT_VARIABLE MATPLOTLIB_FOUND)
endif()
if (Python3_Interpreter_FOUND AND MATPLOTLIB_FOUND EQUAL 0)
  message(STATUS "Convergence plotting available. Setting up convergence test.")
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_convergence.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})  
  add_custom_command(OUTPUT convergence.dat convergence.png
    DEPENDS convergence.bin plot_convergence.py
    COMMAND ${CMAKE_COMMAND} -E time ./convergence.bin
    COMMAND ${Python3_EXECUTABLE} plot_convergence.py convergence.dat convergence.png)
  add_custom_target(convergence
    DEPENDS convergence.png)
else()
  message(WARNING
    "Python, numpy, or matplotlib is unavailable. "
    "Convergence can be run, but no plots will be generated.")
  add_custom_command(OUTPUT convergence.dat
    DEPENDS convergence.bin
    COMMAND ${CMAKE_COMMAND} -E time ./convergence.bin)
  add_custom_target(convergence
    DEPENDS convergence.dat)
endif()
