#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/Interaction.h"
#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/Decay.h"
#include "../../public/SIREN/interactions/pyDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_Decay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<Decay, std::shared_ptr<Decay>, pyDecay, Interaction>(m, "Decay")
        .def(init<>())
        .def("TotalDecayLengthAllFinalStates",&Decay::TotalDecayLengthAllFinalStates)
        .def("TotalDecayLength",&Decay::TotalDecayLength) //;

        .def("TotalDecayWidthAllFinalStates",overload_cast<siren::dataclasses::InteractionRecord const &>(&Decay::TotalDecayWidthAllFinalStates, const_))
        .def("TotalDecayWidthAllFinalStates",overload_cast<siren::dataclasses::ParticleType const &>(&Decay::TotalDecayWidthAllFinalStates, const_))
        .def("TotalDecayWidth",&Decay::TotalDecayWidth)
        .def("DifferentialDecayWidth",&Decay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&Decay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&Decay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&Decay::DensityVariables)
        .def("FinalStateProbability",&Decay::FinalStateProbability)
        .def("SampleFinalState", (void (Decay::*)(siren::dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const)(&Decay::SampleFinalState))
        ;

}
