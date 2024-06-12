#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/CharmDISFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_CharmDISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<CharmDISFromSpline, std::shared_ptr<CharmDISFromSpline>, CrossSection> charmdisfromspline(m, "CharmDISFromSpline");

    charmdisfromspline

        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(self == self)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&CharmDISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&CharmDISFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&CharmDISFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double, double>(&CharmDISFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&CharmDISFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&CharmDISFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&CharmDISFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&CharmDISFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&CharmDISFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&CharmDISFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&CharmDISFromSpline::FinalStateProbability);
}

