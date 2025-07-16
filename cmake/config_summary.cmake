# -------------------------------------------------------------------------------------------------
# Set some colors

string(ASCII 27 Esc)
set(color_reset       "${Esc}[m")
set(color_boldblue    "${Esc}[1;34m")
set(color_boldcyan    "${Esc}[1;36m")
set(color_boldgreen   "${Esc}[1;32m")
set(color_boldgrey    "${Esc}[1;30m")
set(color_boldmagenta "${Esc}[1;35m")
set(color_boldplain   "${Esc}[1m")
set(color_boldred     "${Esc}[1;31m")
set(color_boldyellow  "${Esc}[1;33m")
set(color_cyan        "${Esc}[36m")
set(color_yellow      "${Esc}[33m")

# -------------------------------------------------------------------------------------------------
# A macro for internal use -- see below for the user macros

# Print out a key-value pair (colored differently than a user option)
macro(config_summary_kvc key value color)
    # Figure out spacing
    string(LENGTH ${key} _slen)
    math(EXPR _tlen "${_width} - ${_slen}")
    string(REPEAT " " ${_tlen} _tail)
    # Key
    string(APPEND _summary
        "${color_boldplain}"
        "   ${key}${_tail} : "
        "${color_reset}"
    )
    # Value
    string(APPEND _summary
        "${color}"
        "${value}"
        "${color_reset}"
    )
    string(APPEND _summary "\n")
endmacro()

# -------------------------------------------------------------------------------------------------
# Macros to write the configuration summary.
# -- config_summary_header should always be called before all other config_summary_*
# -- config_summary_print should always be called after all other config_summary_*

# Print a header and do some setup
# -- Arguments:
#    0) display_name: The name of your project for printouts
#    1) cmake_name: The name for your project as known to CMake (to find variables)
#    2) optional: number of characters for variable names (must be large enough for the longest
#       variable name, defaults to 48)
macro(config_summary_header display_name cmake_name)
    # save the terminal width (`tput cols` was not portable, so hard-code)
    set(_termwidth 80)
    # pretty-print a header bar
    math(EXPR _bwid "${_termwidth} - 1")
    string(REPEAT "-" ${_bwid} _bar)
    set(_bar "${color_boldcyan}${_bar}${color_reset}")
    string(APPEND _summary "${_bar}\n")
    # header message
    string(APPEND _summary
        "${color_boldcyan}"
        "${display_name} configuration summary (version ${${cmake_name}_VERSION})\n"
        "${color_reset}"
        )
    # pretty-print a header bar
    string(APPEND _summary "${_bar}\n")
    # Set the width for variable names
    if (${ARGC} GREATER 2)
        set(_width "${ARGV2}")
    else()
        set(_width 48) #default length if none provided
    endif()
endmacro()

# Start a block within the configuration summary
# -- Arguments:
#    0) The title of the block
macro(config_summary_block title)
    string(APPEND _summary
        "${color_boldcyan}"
        "${title}"
        "${color_reset}"
        "\n"
    )
endmacro()

# Print out a key-value pair for information
# -- Arguments:
#    0) the key of the key-value pair
#    1) the value of the key-value pair
macro(config_summary_keyval key value)
    config_summary_kvc("${key}" "${value}" "${color_yellow}")
endmacro()

# Shorthand to print out the name and value of an internal variable
# -- Arguments:
#    0) The name of the variable
macro(config_summary_variable varname)
    config_summary_keyval("${varname}" "${${varname}}")
endmacro()

# Print out a dependency: not found or version
# -- Arguments:
#    0) The display name of the dependency (for printing)
#    1) The CMake name of the dependency (for finding variables)
#    2) optional: condition for whether or not to print this dependency (e.g., if a dependency is
#       only conditionally included, you may only want to conditionally print out whether or not
#       the dependency was found)
macro(config_summary_dependency display_name cmake_name)
    if (${ARGC} GREATER 2)
        set(_condition "${ARGV2}")
    else()
        set(_condition ON) # always print if no condition provided
    endif()
    if (${_condition})
        if (NOT DEFINED ${cmake_name}_VERSION)
            set(_color "${color_boldred}")
            set(_version "not found")
        else()
            set(_color "${color_boldgreen}")
            set(_version "${${cmake_name}_VERSION}")
        endif()
        config_summary_kvc("${display_name}" "${_version}" "${_color}")
    endif()
endmacro()

# Print out the setting used for a user-settable variable.  This differs from
# config_summary_variable mostly in how the output is colored in order to highlight what the user
# turned on or off at a quick glance.
# -- Arguments:
#    0) The name of the variable
macro(config_summary_option varname)
    set(value "${${varname}}")
    if (${value})
        set(_color "${color_boldgreen}")
    else()
        set(_color "${color_boldgrey}")
    endif()
    config_summary_kvc("${varname}" "${value}" "${_color}")
endmacro()

# Finalize the configuration summary and do the actual printing to the terminal
macro(config_summary_print)
    # pretty-print a footer bar to match the header
    string(APPEND _summary "${_bar}")
    # actually print the summary, along with the footer bar
    message(NOTICE
        "${_summary}"
        "${_footerbar}"
    )
endmacro()
