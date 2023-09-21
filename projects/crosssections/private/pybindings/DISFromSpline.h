#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/DISFromSpline.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_DISFromSpline(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    class_<DISFromSpline, std::shared_ptr<DISFromSpline>, CrossSection> disfromspline(m, "DISFromSpline");

    disfromspline
        .def(init<>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def(init<std::string, std::string, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def(init<std::string, std::string, int, double, double, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
        .def(init<std::string, std::string, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
        .def(self == self)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double>(&DISFromSpline::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType>(&DISFromSpline::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DISFromSpline::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<double, double, double, double>(&DISFromSpline::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DISFromSpline::InteractionThreshold)
        .def("GetPossibleTargets",&DISFromSpline::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DISFromSpline::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DISFromSpline::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DISFromSpline::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DISFromSpline::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&DISFromSpline::FinalStateProbability);
}

