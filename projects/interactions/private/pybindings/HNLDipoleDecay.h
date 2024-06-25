#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLDipoleDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLDipoleDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLDipoleDecay, std::shared_ptr<HNLDipoleDecay>, Decay> HNLDipoleDecay(m, "HNLDipoleDecay");

    HNLDipoleDecay
        .def(init<double, std::vector<double>, HNLDipoleDecay::ChiralNature>())
        .def(init<double, std::vector<double>, HNLDipoleDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(init<double, double, HNLDipoleDecay::ChiralNature>())
        .def(init<double, double, HNLDipoleDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(self == self)
        .def("GetHNLMass",&HNLDipoleDecay::GetHNLMass)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDipoleDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&HNLDipoleDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&HNLDipoleDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&HNLDipoleDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&HNLDipoleDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&HNLDipoleDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability",&HNLDipoleDecay::FinalStateProbability)
        ;

    enum_<HNLDipoleDecay::ChiralNature>(HNLDipoleDecay, "ChiralNature")
        .value("Dirac",HNLDipoleDecay::ChiralNature::Dirac)
        .value("Majorana",HNLDipoleDecay::ChiralNature::Majorana)
        .export_values()
        ;
}
