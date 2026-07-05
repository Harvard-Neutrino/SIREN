#include <string>
#include <vector>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SIREN/io/HepMC3Writer.h"
#include "SIREN/dataclasses/InteractionTree.h"

namespace py = pybind11;
using namespace siren::io;

PYBIND11_MODULE(hepmc3, m) {
    py::class_<HepMC3Writer::Options>(m, "HepMC3WriterOptions")
        .def(py::init<>())
        .def_readwrite("siren_version", &HepMC3Writer::Options::siren_version)
        .def_readwrite("weight_names", &HepMC3Writer::Options::weight_names)
        .def_readwrite("provenance", &HepMC3Writer::Options::provenance);

    py::class_<HepMC3Writer>(m, "HepMC3Writer")
        .def(py::init<std::string const &, HepMC3Writer::Options const &>(),
             py::arg("filename"), py::arg("options") = HepMC3Writer::Options())
        .def("write", &HepMC3Writer::Write, py::arg("tree"), py::arg("event_number"))
        .def("close", &HepMC3Writer::Close);

    m.def("SaveInteractionTreesAsHepMC3", &SaveInteractionTreesAsHepMC3,
          py::arg("trees"), py::arg("filename"),
          py::arg("options") = HepMC3Writer::Options(),
          "Write a list of InteractionTrees to a HepMC3 Ascii file.");
}
