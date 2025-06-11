#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLDecay, std::shared_ptr<HNLDecay>, Decay> HNLDecay(m, "HNLDecay");

    HNLDecay
        .def(init<double, std::vector<double>, HNLDecay::ChiralNature>())
        .def(init<double, std::vector<double>, HNLDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(init<double, double, HNLDecay::ChiralNature>())
        .def(init<double, double, HNLDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(self == self)
        .def("GetHNLMass",&HNLDecay::GetHNLMass)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&HNLDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&HNLDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&HNLDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&HNLDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&HNLDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability",&HNLDecay::FinalStateProbability)
        ;

    enum_<HNLDecay::ChiralNature>(HNLDecay, "ChiralNature")
        .value("Dirac",HNLDecay::ChiralNature::Dirac)
        .value("Majorana",HNLDecay::ChiralNature::Majorana)
        .export_values()
        ;
}
