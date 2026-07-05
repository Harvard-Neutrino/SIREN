#include <string>
#include <vector>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SIREN/io/HepMC3Writer.h"
#include "SIREN/io/HepMC3Reader.h"
#include "SIREN/dataclasses/InteractionTree.h"

namespace py = pybind11;
using namespace siren::io;

PYBIND11_MODULE(hepmc3, m) {
    py::class_<HepMC3Writer::Options>(m, "HepMC3WriterOptions")
        .def(py::init<>())
        .def_readwrite("siren_version", &HepMC3Writer::Options::siren_version)
        .def_readwrite("weight_names", &HepMC3Writer::Options::weight_names)
        .def_readwrite("provenance", &HepMC3Writer::Options::provenance)
        // Weight provenance echoed as siren.weights_state; "computed"/"header" are
        // NuHepMC modes, "unweighted" suppresses all NuHepMC.* and siren.fatx.* keys.
        .def_readwrite("weights_state", &HepMC3Writer::Options::weights_state)
        // Extra non-PDG particle codes to declare (NuHepMC G.R.11), merged with
        // SIREN's built-in BSM set. Python type: Dict[int, Tuple[str, str]] mapping
        // code -> (name, description). A code equal to a built-in overrides it.
        .def_readwrite("additional_particle_numbers",
                       &HepMC3Writer::Options::additional_particle_numbers)
        // Run-level generation counts (metadata + FATX normalization); < 0 = unset.
        .def_readwrite("attempted_events", &HepMC3Writer::Options::attempted_events)
        .def_readwrite("accepted_events", &HepMC3Writer::Options::accepted_events)
        // The Injector's EventsToInject seed (pooled-weighting N_i); emitted as
        // siren.events_to_inject when >= 0. < 0 = unset.
        .def_readwrite("events_to_inject", &HepMC3Writer::Options::events_to_inject)
        // Flux-averaged cross section controls (NuHepMC E.C.4 / G.R.6).
        .def_readwrite("fatx_per_atom", &HepMC3Writer::Options::fatx_per_atom)
        .def_readwrite("fatx_partition_by_primary",
                       &HepMC3Writer::Options::fatx_partition_by_primary)
        .def_readwrite("cross_section_unit", &HepMC3Writer::Options::cross_section_unit)
        .def_readwrite("target_scale", &HepMC3Writer::Options::target_scale)
        .def_readwrite("gzip", &HepMC3Writer::Options::gzip);

    py::class_<HepMC3Writer>(m, "HepMC3Writer")
        .def(py::init<std::string const &, HepMC3Writer::Options const &>(),
             py::arg("filename"), py::arg("options") = HepMC3Writer::Options())
        .def("write", &HepMC3Writer::Write, py::arg("tree"), py::arg("event_number"))
        .def("close", &HepMC3Writer::Close);

    m.def("SaveInteractionTreesAsHepMC3", &SaveInteractionTreesAsHepMC3,
          py::arg("trees"), py::arg("filename"),
          py::arg("options") = HepMC3Writer::Options(),
          "Write a list of InteractionTrees to a HepMC3 Ascii file.");

    m.def("LoadInteractionTreesFromHepMC3", &LoadInteractionTreesFromHepMC3,
          py::arg("filename"), py::arg("strict") = true,
          "Read a HepMC3 Ascii file written by SIREN back into InteractionTrees. "
          "With strict=True (default) a missing/mistyped required siren.* attribute "
          "(helicity, primary initial position/time) raises; strict=False tolerates it.");
}
