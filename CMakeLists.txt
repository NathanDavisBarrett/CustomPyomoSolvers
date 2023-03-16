cmake_minimum_required(VERSION 3.4)
project(nathans_simplex_solver)
add_subdirectory(pybind11)
pybind11_add_module(nathans_Simplex_solver_py NathansSimplexSolver.cpp)

