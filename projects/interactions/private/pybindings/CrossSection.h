#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/pyCrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_CrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<CrossSection, std::shared_ptr<CrossSection>, pyCrossSection>(m, "CrossSection")
        .def(init<>())
        .def("__eq__", [](const CrossSection &self, const CrossSection &other){ return self == other; })
        .def("equal", &CrossSection::equal)
        .def("TotalCrossSection", (double (CrossSection::*)(siren::dataclasses::InteractionRecord const &) const)(&CrossSection::TotalCrossSection))
        .def("TotalCrossSectionAllFinalStates", (double (CrossSection::*)(siren::dataclasses::InteractionRecord const &) const)(&CrossSection::TotalCrossSectionAllFinalStates))
        .def("TotalCrossSection", (double (CrossSection::*)(siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType) const)(&CrossSection::TotalCrossSection))
        .def("DifferentialCrossSection", &CrossSection::DifferentialCrossSection)
        .def("InteractionThreshold", &CrossSection::InteractionThreshold)
        .def("SampleFinalState", (void (CrossSection::*)(siren::dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const)(&CrossSection::SampleFinalState))
        .def("GetPossibleTargets", &CrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary", &CrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries", &CrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures", &CrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents", &CrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability", &CrossSection::FinalStateProbability)
        .def("DensityVariables", &CrossSection::DensityVariables)
        ;
}

