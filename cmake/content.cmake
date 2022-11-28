include(FetchContent)

###############################################################
# 
# For dependency management, use the `FetchContent` pattern
#
# NOTE: We seek to replicate the implementation of `FetchContent`
###############################################################

macro(spiner_content_declare pkg_name)
  set(options
    NO_FETCH
  )
  set(one_value_args
    GIT_REPO
    GIT_TAG
    NOTFOUND_MSG
  )
  set(multi_value_args
    COMPONENTS
    EXPECTED_TARGETS
  )

  cmake_parse_arguments(fp "${options}" "${one_value_args}" "${multi_value_args}" "${ARGN}")

  if (NOT CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
      find_package(${pkg_name} COMPONENTS ${fp_COMPONENTS} QUIET )
    if(${pkg_name}_FOUND)
      message(VERBOSE
            "${pkg_name} located with `find_package`"
            "${pkg_name}_DIR: ${${pkg_name}_DIR}"
        )
    else()
      message(VERBOSE
            "${pkg_name} NOT located with `find_package`"
        )
    endif()

    if(fp_NO_FETCH AND NOT ${pkg_name}_FOUND)
      string(JOIN "\n" _scd_error
        "${pkg_name} is requested, but it was not located and is not declared as fetchable."
        "You may want to try:"
        "\t - if ${pkg_name} is installed, set \"-D${pkg_name}_ROOT=<install-dir>\""
      )
      if(fp_NOTFOUND_MSG)
        string(JOIN "\n" _scd_error 
          "${_scd_error}"
          "\t - ${fp_NOTFOUND_MSG}"
        )
      endif()
      message(FATAL_ERROR "${_scd_error}")
      unset(_spc_error)
    endif()
  endif()

  if(fp_NO_FETCH)
    message(VERBOSE
      "\"${pkg_name}\" is specified not fetchable, will rely on `find_package` for population"
    )
    FetchContent_Declare(${pkg_name}
      DOWNLOAD_COMMAND ":"
    )
  else()
    message(VERBOSE 
      "\"${pkg_name}\" is fetchable, will fall-back to git clone [${fp_GIT_REPO}] if other population methods fail"
    )

    FetchContent_Declare(${pkg_name}
      GIT_REPOSITORY ${fp_GIT_REPO}
      GIT_TAG  ${fp_GIT_TAG}
    )
  endif()
  
  list(APPEND SPINER_DECLARED_EXTERNAL_CONTENT ${pkg_name})
  if(fp_EXPECTED_TARGETS)
      list(APPEND SPINER_DECLARED_EXTERNAL_TARGETS ${fp_EXPECTED_TARGETS})
  else()
      list(APPEND SPINER_DECLARED_EXTERNAL_TARGETS "${pkg_name}::${pkg_name}")
  endif()

endmacro()

macro(spiner_content_populate)
  set(options)
  set(one_value_args
    NAMESPACE
  )
  set(multi_value_args
    ENABLE_OPTS
  )

  cmake_parse_arguments(fp "${options}" "${one_value_args}" "${multi_value_args}" "${ARGN}")

  if(fp_ENABLE_OPTS)
    foreach(ext_opt IN LISTS fp_ENABLE_OPTS)
      message(DEBUG "setting \"${ext_opt}\"=ON")
      set(${ext_opt} ON CACHE INTERNAL "")
    endforeach()
  endif()

  message(STATUS
      "Populating dependency targets ${SPINER_DECLARED_EXTERNAL_TARGETS}\n"
      "Calling `FetchContent_MakeAvailable` with ${SPINER_DECLARED_EXTERNAL_CONTENT}\n"
      "This may take a few moments if a download is required...\n"
  )
  FetchContent_MakeAvailable(${SPINER_DECLARED_EXTERNAL_CONTENT})

  foreach(expected_target IN LISTS SPINER_DECLARED_EXTERNAL_TARGETS)
    if(NOT TARGET ${expected_target})
        message(FATAL_ERROR
            "target \"${expected_target}\" was expected, but does not exist after population!"
        )
    endif()
  endforeach()

  set(${fp_NAMESPACE}_POPULATED_TARGETS ${SPINER_DECLARED_EXTERNAL_TARGETS})

  unset(SPINER_DECLARED_EXTERNAL_CONTENT)
  unset(SPINER_DECLARED_EXTERNAL_TARGETS)
endmacro()


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
