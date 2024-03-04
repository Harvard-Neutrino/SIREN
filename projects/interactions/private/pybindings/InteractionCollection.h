#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/InteractionCollection.h"
#include "../../public/LeptonInjector/interactions/Decay.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_InteractionCollection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<InteractionCollection, std::shared_ptr<InteractionCollection>>(m, "InteractionCollection")
        .def(init<>())
        .def(init<LI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>>())
        .def(init<LI::dataclasses::ParticleType, std::vector<std::shared_ptr<Decay>>>())
        .def(init<LI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>, std::vector<std::shared_ptr<Decay>>>())
        .def(self == self)
        .def("GetDecays",&InteractionCollection::GetDecays)
        .def("HasCrossSections",&InteractionCollection::HasCrossSections)
        .def("HasDecays",&InteractionCollection::HasDecays)
        .def("GetCrossSectionsForTarget",&InteractionCollection::GetCrossSectionsForTarget)
        .def("GetCrossSectionsByTarget",&InteractionCollection::GetCrossSectionsByTarget)
        .def("TargetTypes",&InteractionCollection::TargetTypes)
        .def("TotalDecayWidth",&InteractionCollection::TotalDecayWidth)
        .def("TotalDecayLength",&InteractionCollection::TotalDecayLength)
        .def("MatchesPrimary",&InteractionCollection::MatchesPrimary)
        ;
}
