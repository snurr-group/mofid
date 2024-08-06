project(test)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # Release-1.14.0
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
)
target_link_libraries(
    alltests
    mofidtest
    openbabel
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(alltests)
