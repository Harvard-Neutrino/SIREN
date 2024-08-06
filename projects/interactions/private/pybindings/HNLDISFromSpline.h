#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLDISFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLDISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLDISFromSpline, std::shared_ptr<HNLDISFromSpline>, CrossSection> hnldisfromspline(m, "HNLDISFromSpline");

    hnldisfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, double, std::vector<double>, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::vector<char>, std::vector<char>, double, std::vector<double>, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, double, std::vector<double>, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, double, std::vector<double>, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, double, std::vector<double>, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, double, std::vector<double>, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("hnl_mass"),
                arg("dipole_coupling"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
                .def(self == self)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&HNLDISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType>(&HNLDISFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDISFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::Particle::ParticleType, double, double, double, double>(&HNLDISFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&HNLDISFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&HNLDISFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&HNLDISFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&HNLDISFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&HNLDISFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&HNLDISFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&HNLDISFromSpline::FinalStateProbability);
}

