cmake_minimum_required(VERSION 3.9)

project(nathans_simplex_solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS FALSE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Threads REQUIRED)

set_target_properties(nathans_simplex_solver PROPERTIES LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}")

add_subdirectory(pybind11)
pybind11_add_module(nathans_Simplex_solver_py NathansSimplexSolver.cpp)