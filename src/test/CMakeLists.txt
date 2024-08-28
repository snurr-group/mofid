project(test)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG b514bdc898e2951020cbdca1304b75f5950d1f59 # Release-1.15.2
)
FetchContent_MakeAvailable(googletest)

# Set the path containing OpenBabel3Config.cmake, needed for find_package below.
find_path(OpenBabel3_DIR OpenBabel3Config.cmake PATHS
    ${OpenBabel3_DIR}
    "${CMAKE_SOURCE_DIR}/../openbabel/build"
    "/usr/lib/openbabel"
    "/usr/local/lib/openbabel"
)

#
# Find and setup OpenBabel3.
#
find_package(Threads REQUIRED)
find_package(OpenBabel3 REQUIRED)
include_directories(${OpenBabel3_INCLUDE_DIRS})

add_executable(
    alltests
    alltests.cpp
)
target_include_directories(
    alltests
    PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}/includes
)
target_link_libraries(
    alltests
    mofidtest
    openbabel
    gtest
)

include(GoogleTest)
gtest_discover_tests(alltests)
