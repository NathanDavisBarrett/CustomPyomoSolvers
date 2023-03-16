//from https://www.youtube.com/watch?v=_5T70cAXDJ0 Stopped at 2:08

#include<pybind11/pybind11.h>

namespace py = pybind11;

float DummyFunc(float arg1, float arg2) {
    return arg1 + arg2;
}

PYBIND11_MODULE(nathans_Simplex_solver_py, handle) {
    handle.doc() = "The C++ implementation of Nathan's Simplex Solver";
    handle.def("DummyFunc", &DummyFunc);
}