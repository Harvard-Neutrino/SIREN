#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLTwoBodyDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLTwoBodyDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLTwoBodyDecay, std::shared_ptr<HNLTwoBodyDecay>, Decay> HNLTwoBodyDecay(m, "HNLTwoBodyDecay");

    HNLTwoBodyDecay
        .def(init<double, std::vector<double>, HNLTwoBodyDecay::ChiralNature>())
        .def(init<double, std::vector<double>, HNLTwoBodyDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(init<double, double, HNLTwoBodyDecay::ChiralNature>())
        .def(init<double, double, HNLTwoBodyDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(self == self)
        .def("GetHNLMass",&HNLTwoBodyDecay::GetHNLMass)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLTwoBodyDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&HNLTwoBodyDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&HNLTwoBodyDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&HNLTwoBodyDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&HNLTwoBodyDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&HNLTwoBodyDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability",&HNLTwoBodyDecay::FinalStateProbability)
        ;

    enum_<HNLTwoBodyDecay::ChiralNature>(HNLTwoBodyDecay, "ChiralNature")
        .value("Dirac",HNLTwoBodyDecay::ChiralNature::Dirac)
        .value("Majorana",HNLTwoBodyDecay::ChiralNature::Majorana)
        .export_values()
        ;
}
