#
# This script is based on the example used to create static executables
# linking to the static OpenBabel2 library.
#

set(CMAKE_LEGACY_CYGWIN_WIN32 0)  # Hide Cygwin error

# This line is required for cmake backwards compatibility.
cmake_minimum_required(VERSION 2.6)

# Name of your project
project(sbu)
# Create a list of source files (easier to maintain)
set(sources sbu.cpp)
# Set the name for the executable
set(executable_target sbu)

################################################################################

set(LOCAL_OB_DATADIR "${CMAKE_SOURCE_DIR}/ob_datadir"
    CACHE PATH "Install dir for OB data customized for MOFs")
# Set up include file for the data directory
configure_file(${CMAKE_SOURCE_DIR}/config_sbu.h.cmake
  ${CMAKE_BINARY_DIR}/includes/config_sbu.h)
include_directories(${CMAKE_BINARY_DIR}/includes/)


# Set the path containing OpenBabel2Config.cmake, needed for find_package below.
find_path(OpenBabel2_DIR OpenBabel2Config.cmake PATHS
          ${OpenBabel2_DIR}
          "../openbabel/build"
          "/usr/lib/openbabel"
         "/usr/local/lib/openbabel")

#
# Find and setup OpenBabel2.
#
find_package(OpenBabel2 REQUIRED)
include_directories(${OpenBabel2_INCLUDE_DIRS})

# Dependencies
#find_package(LibXml2)

# The executable
add_executable(${executable_target} ${sources})
# Link against imported openbabel target
#target_link_libraries(${executable_target} openbabel ${LIBXML2_LIBRARIES})
target_link_libraries(${executable_target} openbabel)
# install(TARGETS ${executable_target} DESTINATION bin)
