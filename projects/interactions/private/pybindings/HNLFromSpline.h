#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLFromSpline, std::shared_ptr<HNLFromSpline>, CrossSection> disfromspline(m, "HNLFromSpline");

    disfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::set<siren::dataclasses::ParticleType>, std::set<siren::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>>())
        .def(self == self)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType>(&HNLFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&HNLFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&HNLFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&HNLFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&HNLFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&HNLFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&HNLFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&HNLFromSpline::FinalStateProbability);
}

