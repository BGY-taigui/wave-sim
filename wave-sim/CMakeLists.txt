cmake_minimum_required(VERSION 3.12)
project(WaveSimProject)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

file(GLOB_RECURSE SOURCES 
    ${SRC_DIR}/*.cpp
    ${SRC_DIR}/classes/*.cpp
    )

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/classes
    )

add_executable(wave-sim ${SOURCES})

find_package(VTK REQUIRED)
target_link_libraries(wave-sim PRIVATE ${VTK_LIBRARIES})
target_include_directories(wave-sim PRIVATE ${VTK_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
include_directories(wave-sim PRIVATE ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(wave-sim PRIVATE ${ARMADILLO_LIBRARIES})