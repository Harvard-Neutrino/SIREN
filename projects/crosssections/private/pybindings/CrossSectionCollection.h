#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/CrossSectionCollection.h"
#include "../../public/LeptonInjector/crosssections/Decay.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_CrossSectionCollection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    class_<CrossSectionCollection, std::shared_ptr<CrossSectionCollection>>(m, "CrossSectionCollection")
        .def(init<>())
        .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>())
        .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<Decay>>>())
        .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>, std::vector<std::shared_ptr<Decay>>>())
        .def(self == self)
        .def("GetDecays",&CrossSectionCollection::GetDecays)
        .def("HasCrossSections",&CrossSectionCollection::HasCrossSections)
        .def("HasDecays",&CrossSectionCollection::HasDecays)
        .def("GetCrossSectionsForTarget",&CrossSectionCollection::GetCrossSectionsForTarget)
        .def("GetCrossSectionsByTarget",&CrossSectionCollection::GetCrossSectionsByTarget)
        .def("TargetTypes",&CrossSectionCollection::TargetTypes)
        .def("TotalDecayWidth",&CrossSectionCollection::TotalDecayWidth)
        .def("TotalDecayLength",&CrossSectionCollection::TotalDecayLength)
        .def("MatchesPrimary",&CrossSectionCollection::MatchesPrimary)
        ;
}
