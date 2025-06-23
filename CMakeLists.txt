cmake_minimum_required(VERSION 3.0)
project(xmldig2csv)

set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fvar-tracking -fvar-tracking-assignments -Wall -Wextra -Wpedantic")

# LibXML2 for parsing the XMLdig file
find_package (LibXml2 REQUIRED)

# make sure these things only are applied on Linux systems
if (UNIX AND NOT APPLE)
    include_directories (${LIBXML2_INCLUDE_DIR})
endif()

include_directories (${PROJECT_SOURCE_DIR})
add_executable(xmldig2csv main.cpp)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(xmldig2csv PRIVATE DEBUG=1)
else()
    target_compile_definitions(xmldig2csv PRIVATE DEBUG=0)
endif()

target_link_libraries(xmldig2csv ${LIBXML2_LIBRARIES})
target_link_libraries(xmldig2csv "stdc++fs")
target_include_directories(xmldig2csv PRIVATE ${CMAKE_SOURCE_DIR}/third_party/argparse)

# add openMP to paralellize processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("-- OpenMP found")
    target_link_libraries(xmldig2csv OpenMP::OpenMP_CXX)
else()
    message("-- OpenMP NOT found, falling back to single-threaded mode")
endif()
