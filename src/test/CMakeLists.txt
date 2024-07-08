project(test)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # Release-1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(
    example
    example.cpp
)
target_include_directories(
    example
    PRIVATE
    ${CMAKE_SOURCE_DIR}
)
target_link_libraries(
    example
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(example)
