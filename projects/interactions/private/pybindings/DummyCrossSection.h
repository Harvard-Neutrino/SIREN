#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/DummyCrossSection.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_DummyCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<DummyCrossSection, std::shared_ptr<DummyCrossSection>, CrossSection> dummycrosssection(m, "DummyCrossSection");

    dummycrosssection
        .def(init<>())
        .def("_equal", &DummyCrossSection::equal)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DummyCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::ParticleType, double, LI::dataclasses::ParticleType>(&DummyCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DummyCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DummyCrossSection::InteractionThreshold)
        .def("GetPossibleTargets",&DummyCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DummyCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DummyCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DummyCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DummyCrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&DummyCrossSection::FinalStateProbability)
        ;
}

