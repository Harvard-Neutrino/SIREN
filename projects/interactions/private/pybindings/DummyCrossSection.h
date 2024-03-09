#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/DummyCrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_DummyCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<DummyCrossSection, std::shared_ptr<DummyCrossSection>, CrossSection> dummycrosssection(m, "DummyCrossSection");

    dummycrosssection
        .def(init<>())
        .def("_equal", &DummyCrossSection::equal)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&DummyCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType>(&DummyCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&DummyCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DummyCrossSection::InteractionThreshold)
        .def("GetPossibleTargets",&DummyCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DummyCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DummyCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DummyCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DummyCrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&DummyCrossSection::FinalStateProbability)
        ;
}

