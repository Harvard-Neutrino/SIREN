
#include <vector>

#include "../../public/LeptonInjector/crosssections/CrossSectionCollection.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(CrossSectionCollection,m) {
  using namespace LI::crosssections;

  class_<CrossSectionCollection>(m, "CrossSectionCollection")
    .def(init<>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<Decay>>>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>, std::vector<std::shared_ptr<Decay>>>())
    .def("TargetTypes",&CrossSectionCollection::TargetTypes)
    .def("TotalDecayWidth",&CrossSectionCollection::TotalDecayWidth)
    .def("TotalDecayLength",&CrossSectionCollection::TotalDecayLength);
}
