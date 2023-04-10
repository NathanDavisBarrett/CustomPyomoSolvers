cmake_minimum_required(VERSION 3.9)

project(nathans_simplex_solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Threads REQUIRED)

add_subdirectory(pybind11)
pybind11_add_module(nathans_Simplex_solver_py NathansSimplexSolver.cpp)

target_link_properties(nathans_simplex_solver LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}")