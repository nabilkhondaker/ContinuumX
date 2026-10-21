#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "continuumx/lattice_boltzmann.hpp"

namespace py = pybind11;

PYBIND11_MODULE(continuumx_core, m) {
    m.doc() = "ContinuumX C++ core: Lattice Boltzmann D2Q9 and coupling primitives";
    
    py::class_<continuumx::LatticeBoltzmannD2Q9>(m, "LatticeBoltzmannD2Q9")
        .def(py::init<int, int, double>(), py::arg("nx"), py::arg("ny"), py::arg("tau") = 0.6)
        .def("initialize", &continuumx::LatticeBoltzmannD2Q9::initialize,
             py::arg("rho0") = 1.0, py::arg("ux") = 0.0, py::arg("uy") = 0.0)
        .def("collide", &continuumx::LatticeBoltzmannD2Q9::collide)
        .def("stream", &continuumx::LatticeBoltzmannD2Q9::stream)
        .def("step", &continuumx::LatticeBoltzmannD2Q9::step)
        .def("apply_bounce_back", &continuumx::LatticeBoltzmannD2Q9::apply_bounce_back)
        .def("apply_inlet", &continuumx::LatticeBoltzmannD2Q9::apply_inlet)
        .def("apply_outlet", &continuumx::LatticeBoltzmannD2Q9::apply_outlet)
        .def("density", &continuumx::LatticeBoltzmannD2Q9::density)
        .def("velocity", &continuumx::LatticeBoltzmannD2Q9::velocity)
        .def("vorticity", &continuumx::LatticeBoltzmannD2Q9::vorticity)
        .def("add_force", &continuumx::LatticeBoltzmannD2Q9::add_force)
        .def("nx", &continuumx::LatticeBoltzmannD2Q9::nx)
        .def("ny", &continuumx::LatticeBoltzmannD2Q9::ny)
        .def("tau", &continuumx::LatticeBoltzmannD2Q9::tau)
        .def("omega", &continuumx::LatticeBoltzmannD2Q9::omega);
}
