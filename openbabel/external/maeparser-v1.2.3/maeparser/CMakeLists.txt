set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)

project(maeparser)

option(MAEPARSER_RIGOROUS_BUILD "Abort the build if the compiler issues any warnings" OFF )
option(MAEPARSER_BUILD_TESTS "Whether test executables should be built" ON)
option(MAEPARSER_BUILD_SHARED_LIBS "Build maeparser as a shared library (turn off for a static one)" ON)

if(MSVC)
    # C4251 disables warnings for export STL containers as arguments (returning a vector of things)
    add_definitions(/wd4251 /wd4275 /wd4996 /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DBOOST_ALL_NO_LIB)
endif(MSVC)

if(MAEPARSER_RIGOROUS_BUILD)
if(MSVC)
    add_definitions(/WX)
else(MSVC)
    add_definitions(-Wall -Wextra -Werror)
endif(MSVC)
endif(MAEPARSER_RIGOROUS_BUILD)

file(GLOB mae_sources "*.cpp")

if(MAEPARSER_BUILD_SHARED_LIBS)
    add_library(maeparser SHARED ${mae_sources})
    target_compile_definitions(maeparser PRIVATE "IN_MAEPARSER")
    target_compile_definitions(maeparser PRIVATE "BOOST_ALL_DYN_LINK")
    set_property(TARGET maeparser PROPERTY CXX_VISIBILITY_PRESET "hidden")
else(MAEPARSER_BUILD_SHARED_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
    if(NOT MSVC)
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()
    add_library(maeparser STATIC ${mae_sources})
    target_compile_definitions(maeparser PRIVATE "STATIC_MAEPARSER")
endif(MAEPARSER_BUILD_SHARED_LIBS)

find_package(Boost COMPONENTS iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(boost_libs ${Boost_LIBRARIES})

find_package(Boost COMPONENTS zlib QUIET)
if(Boost_ZLIB_FOUND)
    set(boost_libs ${boost_libs} ${Boost_LIBRARIES})
    message(STATUS "Using Boost zlib module for iostreams dependency.")
else(Boost_ZLIB_FOUND)
    find_package(ZLIB REQUIRED)
    set(boost_libs ${boost_libs} ${ZLIB_LIBRARIES})
    message(STATUS "Using zlib library for iostreams dependency.")
endif(Boost_ZLIB_FOUND)

target_link_libraries(maeparser ${boost_libs})

SET_TARGET_PROPERTIES (maeparser
    PROPERTIES
       VERSION 1.2.3
       SOVERSION 1
)

target_include_directories(maeparser PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)

install(TARGETS maeparser
    EXPORT maeparser-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin)

INSTALL(EXPORT maeparser-targets
    FILE ${PROJECT_NAME}-config.cmake
    DESTINATION lib/cmake)

file(GLOB mae_headers "*.hpp")
install(FILES ${mae_headers} DESTINATION include/maeparser)

# Tests
if (MAEPARSER_BUILD_TESTS)
    set(MEMORYCHECK_COMMAND_OPTIONS "${MEMORYCHECK_COMMAND_OPTIONS} --tool=memcheck \
            --time-stamp=yes --num-callers=20 --gen-suppressions=all --leak-check=full \
            --show-reachable=yes --trace-children=yes --error-exitcode=29")

    include(CTest)
    add_subdirectory(test)
endif()

