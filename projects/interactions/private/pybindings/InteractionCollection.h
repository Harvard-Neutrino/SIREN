#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/InteractionCollection.h"
#include "../../public/SIREN/interactions/Decay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_InteractionCollection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace SI::interactions;

    class_<InteractionCollection, std::shared_ptr<InteractionCollection>>(m, "InteractionCollection")
        .def(init<>())
        .def(init<SI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>>())
        .def(init<SI::dataclasses::ParticleType, std::vector<std::shared_ptr<Decay>>>())
        .def(init<SI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>, std::vector<std::shared_ptr<Decay>>>())
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
