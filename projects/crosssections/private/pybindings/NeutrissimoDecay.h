#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/NeutrissimoDecay.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_NeutrissimoDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    class_<NeutrissimoDecay, std::shared_ptr<NeutrissimoDecay>, Decay> neutrissimodecay(m, "NeutrissimoDecay");

    neutrissimodecay
        .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature>())
        .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature, std::set<LI::dataclasses::Particle::ParticleType> const &>())
        .def(init<double, double, NeutrissimoDecay::ChiralNature>())
        .def(init<double, double, NeutrissimoDecay::ChiralNature, std::set<LI::dataclasses::Particle::ParticleType> const &>())
        .def(self == self)
        .def("GetHNLMass",&NeutrissimoDecay::GetHNLMass)
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&NeutrissimoDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&NeutrissimoDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&NeutrissimoDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&NeutrissimoDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&NeutrissimoDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&NeutrissimoDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability",&NeutrissimoDecay::FinalStateProbability)
        ;

    enum_<NeutrissimoDecay::ChiralNature>(neutrissimodecay, "ChiralNature")
        .value("Dirac",NeutrissimoDecay::ChiralNature::Dirac)
        .value("Majorana",NeutrissimoDecay::ChiralNature::Majorana)
        .export_values()
        ;
}
