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
        .def("__eq__", [](const Decay &self, const Decay &other){ return self == other; })
        .def("equal", &Decay::equal)
        .def("TotalDecayWidthAllFinalStates", &Decay::TotalDecayWidthAllFinalStates)
        .def("TotalDecayWidth", (double (Decay::*)(siren::dataclasses::ParticleType) const)(&Decay::TotalDecayWidth))
        .def("TotalDecayWidth", (double (Decay::*)(siren::dataclasses::InteractionRecord const &) const)(&Decay::TotalDecayWidth))
        .def("TotalDecayLengthAllFinalStates", &Decay::TotalDecayLengthAllFinalStates)
        .def("TotalDecayLength", &Decay::TotalDecayLength)
        .def("DifferentialDecayWidth", &Decay::DifferentialDecayWidth)
        .def("SampleFinalState", &Decay::SampleFinalState)
        .def("SampleDecayTime", &Decay::SampleDecayTime)
        .def("SecondaryMasses", &Decay::SecondaryMasses)
        .def("SecondaryHelicities", &Decay::SecondaryHelicities)
        .def("GetPossibleSignatures", &Decay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent", &Decay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability", &Decay::FinalStateProbability)
        .def("DensityVariables", &Decay::DensityVariables)
        .def("Topology", &Decay::Topology)
        .def("Measure", &Decay::Measure)
        .def("Convention", &Decay::Convention);

}
