# CustomPyomoSolvers

## Branch Overview
SimplexOnly: If you only want to use the Simplex Solver (no integer or nonlineaer capabilities).
MILPExtension: If you want to solve MILP problems (no nonlinear capabilities)

## How to use:
I'd recommend installing these solvers on a linux machine. I've never tested or designed them to run on anything else.

These solvers are written in c++ with direct access from Python via pybind11.

For a description of pybind11, check out their [github repo](https://github.com/pybind/pybind11).

1. **Make sure you have cmake installed**. For a description on how to install cmake see their [website](https://cmake.org/). Note that if you're using a managed cluster, the cluster managers likely already have cmake installed since it's pretty common. See if it's installed using `module avail` and `module load ____`.
2. **Clone this repo**. From the command line, use `cd ` to navigate the whereever you'd like this repo cloned to. Then use `git clone git@github.com:NathanDavisBarrett/CustomPyomoSolvers.git` to clone this repository.
3. **Create a build directory**. Create a build directory within the cloned repo using `cd CustomPyomoSolvers; mkdir build`.
4. **Build (compile) the c++ code**. Use the following commands to build the c++ code. `cd build; cmake ..; make`. This will create a python module file that can be imported by your original script. It will have a ".so" file type and will be located within the build directory.
5. **Copy the python/c++ module (the .so file from the previous step) and extra python scripts to whichever directory contains the python script you'd like to use these solvers in**. For example, if I wanted to use NathansSimplexSolver, I'd need to copy both build/WhateverTheSOFileIsNamed.so, NathansSimplexSolver.py, and CustomSolverResources.py to your working directory.
6. **Import the python module into your script**. The .py file (which will call the .so file and CustomSolverResources.py) is the file you'd like to import into your script. For example, if I wanted to use NathansSimplexSolver, I'd call `import NathansSimplexSolver` in my python script.
7. **Tell Pyomo to use the custom solver**. This can be done just like any other solver in pyomo in the pyomo.environ.SolverFacotry call, pass the name of the solver you imported as a string. For example, if I wanted to use NathansSimplexSolver, I'd call `pyomo.environ.SolverFactory('NathansSimplexSolver')`.
