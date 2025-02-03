#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/QuarkDISFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_QuarkDISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<QuarkDISFromSpline, std::shared_ptr<QuarkDISFromSpline>, CrossSection> quarkdisfromspline(m, "QuarkDISFromSpline");

    quarkdisfromspline

        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("quark_type"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::vector<char>, std::vector<char>, int, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_data"),
                arg("differential_xs_data"),
                arg("interaction"),
                arg("quark_type"),
                arg("target_mass"),
                arg("minimum_Q2"),
                arg("primary_types"),
                arg("target_types"),
                arg("units") = std::string("cm"))
        .def(init<std::string, std::string, int, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("quark_type"),
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
        .def(init<std::string, std::string, int, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>, std::string>(),
                arg("total_xs_filename"),
                arg("differential_xs_filename"),
                arg("interaction"),
                arg("quark_type"),
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
        .def("SetInteractionType",&QuarkDISFromSpline::SetInteractionType)
        .def("SetQuarkType",&QuarkDISFromSpline::SetQuarkType)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&QuarkDISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&QuarkDISFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&QuarkDISFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double, double>(&QuarkDISFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&QuarkDISFromSpline::InteractionThreshold)
        .def("FragmentationFraction",&QuarkDISFromSpline::FragmentationFraction)
        .def("GetPossibleTargets",&QuarkDISFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&QuarkDISFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&QuarkDISFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&QuarkDISFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&QuarkDISFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&QuarkDISFromSpline::FinalStateProbability);
}

