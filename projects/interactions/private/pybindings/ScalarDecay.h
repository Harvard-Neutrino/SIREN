#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/ScalarDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_ScalarDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<ScalarDecay, std::shared_ptr<ScalarDecay>, Decay> scalardecay(m, "ScalarDecay");

    scalardecay
        .def(init<double, double)
        .def(self == self)
        .def("GetMass",&ScalarDecay::GetMass)
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::InteractionRecord const &>(&ScalarDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<siren::dataclasses::ParticleType>(&ScalarDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&ScalarDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&ScalarDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&ScalarDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&ScalarDecay::GetPossibleSignaturesFromParent)
        .def("FinalStateProbability",&ScalarDecay::FinalStateProbability)
        ;
}
