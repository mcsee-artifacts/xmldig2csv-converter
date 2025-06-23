cmake_minimum_required(VERSION 3.0)
project(xmldig2csv)

set(default_build_type "Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fvar-tracking -fvar-tracking-assignments")

# LibXML2 for parsing the XMLdig file
find_package (LibXml2 REQUIRED)

# make sure these things only are applied on Linux systems
if (UNIX AND NOT APPLE)
    find_package (PkgConfig REQUIRED)
    pkg_check_modules (LIBXMLXX REQUIRED IMPORTED_TARGET libxml++-2.6)
    include_directories (${LIBXML2_INCLUDE_DIR})
    include_directories ("/usr/include/libxml++2.6/")
endif()

include_directories (${PROJECT_SOURCE_DIR})
add_executable(xmldig2csv main.cpp)
target_link_libraries(xmldig2csv ${LIBXML2_LIBRARIES})
target_link_libraries(xmldig2csv "stdc++fs")

# add openMP to paralellize processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
    target_link_libraries(xmldig2csv OpenMP::OpenMP_CXX)
else()
    message("OpenMP NOT found, falling back to single-threaded mode.")
endif()
