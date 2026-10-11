# © 2021-2026. Triad National Security, LLC. All rights reserved.  This program was
# produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC
# for the U.S.  Department of Energy/National Nuclear Security Administration.
# All rights in the program are reserved by Triad National Security, LLC, and
# the U.S. Department of Energy/National Nuclear Security Administration. The
# Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to
# reproduce, prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.

# Generative AI was used to assist with modifications to this file.

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.7.1 # TODO: make sure this is correct
)

list(APPEND SPINER_TEST_DEPS Catch2)
list(APPEND SPINER_TEST_DEP_TARGETS Catch2::Catch2)

add_library(spiner_test INTERFACE)
# instructions for local content
# https://cmake.org/cmake/help/latest/module/FetchContent.html#variable:FETCHCONTENT_SOURCE_DIR_%3CuppercaseName%3E
# TLDR: set FETCHCONTENT_SOURCE_DIR_KOKKOS=path/to/local/checkout
if(SPINER_TEST_USE_KOKKOS)
  # If Kokkos exists, ignore
  if(NOT TARGET Kokkos::kokkos)
    # put together a list of kokkos configure options TODO: put this into fetch
    list(APPEND _spiner_content_opts "Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION")
    if(SPINER_TEST_USE_KOKKOS_CUDA)
      list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA")
      list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA_LAMBDA")
      list(APPEND _spiner_content_opts "Kokkos_ENABLE_CUDA_CONSTEXPR")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      list(APPEND _spiner_content_opts "Kokkos_ENABLE_DEBUG")
      list(APPEND _spiner_content_opts "Kokkos_ENABLE_DEBUG_BOUNDS_CHECK")
      list(APPEND _spiner_content_opts
           "Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK")
    endif()

    FetchContent_Declare(
      Kokkos
      GIT_REPOSITORY https://github.com/kokkos/kokkos.git
      GIT_TAG 5.1.0 # TODO: make sure this is correct
    )

    foreach(ext_opt ${_spiner_content_opts})
      message(DEBUG "setting \"${ext_opt}\"=ON")
      set(${ext_opt}
          ON
          CACHE INTERNAL "")
    endforeach()

    list(APPEND SPINER_TEST_DEPS Kokkos)
  endif()
  
  list(APPEND SPINER_TEST_DEP_TARGETS Kokkos::kokkos)
  target_compile_definitions(spiner_test INTERFACE PORTABILITY_STRATEGY_KOKKOS)
endif()

# populate test content
FetchContent_MakeAvailable(${SPINER_TEST_DEPS})

target_link_libraries(spiner_test INTERFACE ${SPINER_TEST_DEP_TARGETS})
# add unit tests
add_executable(test.bin test.cpp)
target_link_libraries(test.bin PRIVATE spiner::spiner spiner_test)

list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
include(Catch)
catch_discover_tests(test.bin)

# Benchmarking
add_executable(benchmark.bin benchmark.cpp)
target_link_libraries(benchmark.bin PRIVATE spiner::spiner spiner_test)

add_executable(nonuniform_benchmark.bin nonuniform_benchmark.cpp)
target_link_libraries(nonuniform_benchmark.bin PRIVATE spiner::spiner spiner_test)

# Metadata is refreshed at build time, including source hashes for dirty builds.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_executable(rational_benchmark.bin rational_benchmark.cpp)
target_link_libraries(rational_benchmark.bin PRIVATE spiner::spiner spiner_test)
target_include_directories(rational_benchmark.bin PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" _rational_config)
set(_rational_flags "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_rational_config}}")
add_custom_target(rational_benchmark_metadata
  COMMAND ${CMAKE_COMMAND}
    "-DSOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/rational_benchmark_metadata.hpp"
    "-DBINARY_DIR=${PROJECT_BINARY_DIR}"
    "-DCOMPILER=${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})"
    "-DBUILD_TYPE=$<CONFIG>"
    "-DFLAGS=${_rational_flags}"
    "-DPOC_DIR=${ports-of-call_SOURCE_DIR}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/rational_benchmark_metadata.cmake"
  BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/rational_benchmark_metadata.hpp"
  VERBATIM)
add_dependencies(rational_benchmark.bin rational_benchmark_metadata)

# Small correctness/CLI checks only; no performance thresholds.
foreach(distribution interior-shuffled interior-ordered boundary)
  add_test(NAME rational_benchmark_${distribution}
    COMMAND rational_benchmark.bin --queries 256 --passes 1 --trials 1
      --warmup 0 --nx1d 7 --nx2d 7 --ny2d 4
      --distribution ${distribution}
      --output "${CMAKE_CURRENT_BINARY_DIR}/rational-smoke-${distribution}.csv")
  add_test(NAME rational_benchmark_fast_${distribution}
    COMMAND rational_benchmark.bin --grid fast-nonuniform
      --queries 256 --passes 1 --trials 1 --warmup 0
      --nx1d 7 --nx2d 7 --ny2d 4 --distribution ${distribution}
      --output "${CMAKE_CURRENT_BINARY_DIR}/rational-fast-smoke-${distribution}.csv")
endforeach()
add_test(NAME rational_benchmark_invalid_queries
  COMMAND rational_benchmark.bin --queries 0)
set_tests_properties(rational_benchmark_invalid_queries PROPERTIES WILL_FAIL TRUE)

# Convergence test
add_executable(convergence.bin convergence.cpp)
target_link_libraries(convergence.bin PRIVATE spiner::spiner spiner_test)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  message(STATUS "Checking for numpy, matplotlib.")
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, matplotlib"
                  RESULT_VARIABLE MATPLOTLIB_FOUND)
endif()
if(Python3_Interpreter_FOUND AND MATPLOTLIB_FOUND EQUAL 0)
  message(STATUS "Convergence plotting available. Setting up convergence test.")
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_convergence.py
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_command(
    OUTPUT convergence.dat convergence.png
    DEPENDS convergence.bin plot_convergence.py
    COMMAND ${CMAKE_COMMAND} -E time ./convergence.bin
    COMMAND ${Python3_EXECUTABLE} plot_convergence.py convergence.dat
            convergence.png)
  add_custom_target(convergence DEPENDS convergence.png)
else()
  message(WARNING "Python, numpy, or matplotlib is unavailable. "
                  "Convergence can be run, but no plots will be generated.")
  add_custom_command(
    OUTPUT convergence.dat
    DEPENDS convergence.bin
    COMMAND ${CMAKE_COMMAND} -E time ./convergence.bin)
  add_custom_target(convergence DEPENDS convergence.dat)
endif()
