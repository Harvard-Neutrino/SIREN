#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/Decay.h"
#include "../../public/SIREN/interactions/ElectroweakDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_ElectroweakDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<ElectroweakDecay, std::shared_ptr<ElectroweakDecay>, Decay>(m, "ElectroweakDecay")
        .def(init<std::set<siren::dataclasses::ParticleType> const &>())
        .def(self == self)
        .def("TotalDecayWidthAllFinalStates", overload_cast<siren::dataclasses::InteractionRecord const &>(&ElectroweakDecay::TotalDecayWidthAllFinalStates, const_))
        .def("TotalDecayWidth", overload_cast<siren::dataclasses::ParticleType>(&ElectroweakDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth", (double (siren::interactions::ElectroweakDecay::*)(siren::dataclasses::InteractionRecord const &) const)&ElectroweakDecay::TotalDecayWidth)
        .def("DifferentialDecayWidth", &ElectroweakDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures", &ElectroweakDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent", &ElectroweakDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability", &ElectroweakDecay::FinalStateProbability);
}
