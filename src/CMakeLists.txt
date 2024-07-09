#
# This script is based on the example used to create static executables
# linking to the static OpenBabel2 library.
#

set(CMAKE_LEGACY_CYGWIN_WIN32 0)  # Hide Cygwin error

# Minimum version required for FetchContent
cmake_minimum_required(VERSION 3.11)

# Name of your project
project(sbu)
# Create a list of source files (easier to maintain)
# set(sources sbu.cpp)
# Set the name for the executable
# set(executable_target sbu)
# Use a loop like openbabel/tools/CMakeLists.txt

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()
add_subdirectory(test)
add_library(mofidtest
    STATIC
    obdetails.cpp
)

# tools that do not require external MOFid headers/sources:
set(tools
        sobgrep
        searchdb
        tsfm_smiles
   )
if (EMSCRIPTEN)
  set (tools searchdb)  # disable extraneous tools from JS build
endif (EMSCRIPTEN)
# tools that do require external headers/sources:
set(linked_tools sbu)
set(mofid_includes
        obdetails.cpp
        deconstructor.cpp
        framework.cpp
        periodic.cpp
        pseudo_atom.cpp
        topology.cpp
        virtual_mol.cpp
)

################################################################################

set(LOCAL_OB_DATADIR "${CMAKE_SOURCE_DIR}/../openbabel/data"
    CACHE PATH "Install dir for OB data, no customization needed for MOFs")
# Set up include file for the data directory
configure_file(${CMAKE_SOURCE_DIR}/config_sbu.h.cmake
  ${CMAKE_BINARY_DIR}/includes/config_sbu.h)
include_directories(${CMAKE_BINARY_DIR}/includes/)


# Set the path containing OpenBabel2Config.cmake, needed for find_package below.
find_path(OpenBabel2_DIR OpenBabel2Config.cmake PATHS
          ${OpenBabel2_DIR}
          "${CMAKE_SOURCE_DIR}/../openbabel/build"
          "/usr/lib/openbabel"
         "/usr/local/lib/openbabel")

#
# Find and setup OpenBabel2.
#
find_package(OpenBabel2 REQUIRED)
include_directories(${OpenBabel2_INCLUDE_DIRS})

# Dependencies
#find_package(LibXml2)

# Two build types are available: Release (DEFAULT), Debug
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g3 -fno-inline")
endif()

foreach(tool ${tools})
  # The executable
  add_executable(${tool} ${tool}.cpp)  # no ${mofid_includes}
  # Link against imported openbabel target
  #target_link_libraries(${executable_target} openbabel ${LIBXML2_LIBRARIES})
  target_link_libraries(${tool} openbabel)
  # install(TARGETS ${executable_target} DESTINATION bin)
endforeach(tool)
foreach(linked_tool ${linked_tools})
  add_executable(${linked_tool} ${linked_tool}.cpp ${mofid_includes})
  target_link_libraries(${linked_tool} openbabel)
endforeach(linked_tool)


# Set up DLL's for Cygwin.
# Windows cannot find the paths to DLL's unless $PATH is modified or
# the DLL files are in the same directory as the exe.
# This section will fix DLL issues on Cygwin but not Windows in general.
IF (CYGWIN)
  set (ob_dlls cyginchi-0.dll cygopenbabel-5.dll)
  FOREACH(ob_dll ${ob_dlls})
    add_custom_command(
      TARGET sbu POST_BUILD COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${OpenBabel2_DIR}/bin/${ob_dll}" ${CMAKE_BINARY_DIR}
    )
  ENDFOREACH()
ENDIF (CYGWIN)

# Set up Emscripten-specific linker flags
IF (EMSCRIPTEN)
  # We could try to use target_compile_options to append more compile flags, but it looks like cmake
  # combines them together instead of two separate flags per
  # https://gitlab.kitware.com/cmake/cmake/issues/15826
  # Instead of using a new enough cmake for the SHELL: syntax, let's
  # pass the arguments another way
  #target_compile_options(searchdb PUBLIC "-s EXPORTED_FUNCTIONS='[_runSearchc]'")  # appends more compile flags
  #target_compile_options(searchdb PUBLIC "-s EXPORTED_FUNCTIONS=\"['_runSearchc']\"")  # appends more compile flags
  # New idea per https://github.com/emscripten-core/emscripten/issues/4398
  set_target_properties(searchdb PROPERTIES LINK_FLAGS "-s EXPORTED_FUNCTIONS=\"['_runSearchc']\"")
  set_target_properties(sbu PROPERTIES LINK_FLAGS "-s EXPORTED_FUNCTIONS=\"['_analyzeMOFc']\"")
ENDIF (EMSCRIPTEN)
