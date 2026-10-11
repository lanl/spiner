# Copyright 2026 Triad National Security, LLC. All rights reserved.
# Generative AI was used to assist with this file.
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${SOURCE_DIR}"
  OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
execute_process(COMMAND git status --porcelain --untracked-files=no
  WORKING_DIRECTORY "${SOURCE_DIR}" OUTPUT_VARIABLE dirty ERROR_QUIET)
if(dirty)
  set(dirty "yes")
else()
  set(dirty "no (tracked files only; see source hashes)")
endif()
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${POC_DIR}"
  OUTPUT_VARIABLE poc_revision OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
set(info "compiler=${COMPILER}\nbuild_type=${BUILD_TYPE}\ncmake_cxx_flags=${FLAGS}\nsource_commit=${revision}\ndirty=${dirty}\nports_of_call_revision=${poc_revision}\n")
foreach(path test/rational_benchmark.cpp test/CMakeLists.txt test/rational_benchmark_metadata.cmake spiner/interpolation_stencils.hpp spiner/databox.hpp spiner/regular_grid_1d.hpp spiner/nonuniform_grid_1d.hpp spiner/fast_nonuniform_grid_1d.hpp spiner/piecewise_grid_1d.hpp)
  file(SHA256 "${SOURCE_DIR}/${path}" hash)
  string(APPEND info "sha256_${path}=${hash}\n")
endforeach()
if(EXISTS "${BINARY_DIR}/compile_commands.json")
  file(READ "${BINARY_DIR}/compile_commands.json" commands)
  string(JSON count LENGTH "${commands}")
  math(EXPR last "${count} - 1")
  foreach(i RANGE 0 ${last})
    string(JSON source GET "${commands}" ${i} file)
    if(source STREQUAL "${SOURCE_DIR}/test/rational_benchmark.cpp")
      string(JSON command GET "${commands}" ${i} command)
      string(APPEND info "compile_command=${command}\n")
    endif()
  endforeach()
endif()
# copy_if_different avoids needless recompilation when the metadata is unchanged.
file(WRITE "${OUTPUT}.tmp" "// Copyright 2026 Triad National Security, LLC. All rights reserved.\n// Generative AI assisted this generated file.\nconstexpr const char *benchmarkBuildMetadata = R\"metadata(${info})metadata\";\n")
execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${OUTPUT}.tmp" "${OUTPUT}")
