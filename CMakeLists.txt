cmake_minimum_required(VERSION 3.22)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELASE "-O2")
set(CMAKE_BUILD_PARALLEL_LEVEL 8)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

file(GLOB SOURCES
    src/*.cpp
)

project(Matrix)
include_directories(${PROJECT_SOURCE_DIR}/inc ${PROJECT_SOURCE_DIR}/src)
add_executable(Matrix main.cpp ${SOURCES})