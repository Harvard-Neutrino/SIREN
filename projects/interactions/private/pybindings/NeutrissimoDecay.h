#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/NeutrissimoDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_NeutrissimoDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<NeutrissimoDecay, std::shared_ptr<NeutrissimoDecay>, Decay> neutrissimodecay(m, "NeutrissimoDecay");

    neutrissimodecay
        .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature>())
        .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(init<double, double, NeutrissimoDecay::ChiralNature>())
        .def(init<double, double, NeutrissimoDecay::ChiralNature, std::set<siren::dataclasses::ParticleType> const &>())
        .def(self == self)
        .def("GetHNLMass",&NeutrissimoDecay::GetHNLMass)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&NeutrissimoDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&NeutrissimoDecay::TotalDecayWidth, const_))
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
