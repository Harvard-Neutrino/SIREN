#include <set>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/DISFromSpline.h"
#include "../../public/SIREN/interactions/GammaDISFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_GammaDISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<GammaDISFromSpline, std::shared_ptr<GammaDISFromSpline>, DISFromSpline>
        gammadisfromspline(m, "GammaDISFromSpline");

    gammadisfromspline
        .def(init<>())
        .def(init<std::string, std::string,
                  std::vector<siren::dataclasses::ParticleType>,
                  std::vector<siren::dataclasses::ParticleType>,
                  std::string, std::string, std::string>(),
                arg("differential_xs_filename"),
                arg("total_xs_filename"),
                arg("primary_types"),
                arg("target_types"),
                arg("shower_model_path"),
                arg("pythia_events_path"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double,
                  std::vector<siren::dataclasses::ParticleType>,
                  std::vector<siren::dataclasses::ParticleType>,
                  std::string, std::string, std::string>(),
                arg("differential_xs_filename"),
                arg("total_xs_filename"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("shower_model_path"),
                arg("pythia_events_path"),
                arg("units") = std::string("cm"))
        .def(self == self)
        .def("Name", &GammaDISFromSpline::Name)
        .def("LoadShowerModel", &GammaDISFromSpline::LoadShowerModel)
        .def("LoadPythiaEvents", &GammaDISFromSpline::LoadPythiaEvents)
        .def("HasSpecies", &GammaDISFromSpline::HasSpecies)
        .def("YieldMean", &GammaDISFromSpline::YieldMean)
        .def("InteractionType", &GammaDISFromSpline::InteractionType)   // 1 = CC, 2 = NC
        .def("SetInteractionType", &GammaDISFromSpline::SetInteractionType)
        .def("SampleFinalState", &GammaDISFromSpline::SampleFinalState);
}
