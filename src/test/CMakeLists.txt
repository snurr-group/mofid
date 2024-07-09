project(test)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # Release-1.14.0
)
FetchContent_MakeAvailable(googletest)

# Set the path containing OpenBabel2Config.cmake, needed for find_package below.
find_path(OpenBabel2_DIR OpenBabel2Config.cmake PATHS
    ${OpenBabel2_DIR}
    "${CMAKE_SOURCE_DIR}/../openbabel/build"
    "/usr/lib/openbabel"
    "/usr/local/lib/openbabel"
)

#
# Find and setup OpenBabel2.
#
find_package(OpenBabel2 REQUIRED)
include_directories(${OpenBabel2_INCLUDE_DIRS})

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
