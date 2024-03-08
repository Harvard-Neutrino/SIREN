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
    using namespace SI::interactions;

    class_<HNLFromSpline, std::shared_ptr<HNLFromSpline>, CrossSection> disfromspline(m, "HNLFromSpline");

    disfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<SI::dataclasses::ParticleType>, std::set<SI::dataclasses::ParticleType>>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<SI::dataclasses::ParticleType>, std::vector<SI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::set<SI::dataclasses::ParticleType>, std::set<SI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::set<SI::dataclasses::ParticleType>, std::set<SI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::vector<SI::dataclasses::ParticleType>, std::vector<SI::dataclasses::ParticleType>>())
        .def(init<std::string, std::string, std::vector<SI::dataclasses::ParticleType>, std::vector<SI::dataclasses::ParticleType>>())
        .def(self == self)
        .def("TotalCrossSection",overload_cast<SI::dataclasses::InteractionRecord const &>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<SI::dataclasses::ParticleType, double>(&HNLFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<SI::dataclasses::ParticleType, double, SI::dataclasses::ParticleType>(&HNLFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<SI::dataclasses::InteractionRecord const &>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double>(&HNLFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&HNLFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&HNLFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&HNLFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&HNLFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&HNLFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&HNLFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&HNLFromSpline::FinalStateProbability);
}

