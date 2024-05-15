#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "../../public/SIREN/interactions/Decay.h"
#include "../../public/SIREN/interactions/DarkNewsDecay.h"

void register_DarkNewsDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<DarkNewsDecay, std::shared_ptr<DarkNewsDecay>, Decay, siren::interactions::pyDarkNewsDecay> DarkNewsDecay(m, "DarkNewsDecay");

    DarkNewsDecay
        .def(init<>())
        .def("__eq__", [](const siren::interactions::DarkNewsDecay &self, const siren::interactions::DarkNewsDecay &other){ return self == other; })
        .def("equal", &siren::interactions::DarkNewsDecay::equal)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&DarkNewsDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&DarkNewsDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&DarkNewsDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&DarkNewsDecay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&DarkNewsDecay::DensityVariables)
        .def("FinalStateProbability",&DarkNewsDecay::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsDecay::SampleFinalState)
        .def("SampleRecordFromDarkNews",&DarkNewsDecay::SampleRecordFromDarkNews)
        ;

    RegisterTrampolinePickleMethods(DarkNewsDecay,pyDarkNewsDecay)
}
