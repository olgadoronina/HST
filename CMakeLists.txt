cmake_minimum_required(VERSION 3.6)
project(main)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
find_package(MPI REQUIRED)

include_directories(
        include
        /usr/include/mpi/)


set(SOURCE_FILES
        source/init.cpp
        source/io.cpp
        source/par.cpp
        source/utils.cpp
        source/visua.cpp
        source/xml_writer.cpp
        programs/main.cpp)

add_executable(main ${SOURCE_FILES})