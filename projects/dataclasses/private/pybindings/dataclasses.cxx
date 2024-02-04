
#include <vector>
#include <set>

#include "../../public/LeptonInjector/dataclasses/Particle.h"
#include "../../public/LeptonInjector/dataclasses/ParticleID.h"
#include "../../public/LeptonInjector/dataclasses/ParticleType.h"
#include "../../public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../public/LeptonInjector/dataclasses/InteractionTree.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



using namespace pybind11;

PYBIND11_MODULE(dataclasses,m) {
  using namespace LI::dataclasses;

  class_<Particle, std::shared_ptr<Particle>> particle(m, "Particle");

  particle.def(init<>())
          .def(init<ParticleType>())
          .def_readwrite("type",&Particle::type)
          .def_readwrite("energy",&Particle::energy)
          .def("GetMass",&Particle::GetMass)
          .def("HasMass",&Particle::HasMass)
          .def("GetTypeString",&Particle::GetTypeString);

  enum_<ParticleType>(particle, "ParticleType")
#define X(a) .value( #a , ParticleType:: a )
#include "../../public/LeptonInjector/dataclasses/ParticleTypes.def"
#undef X
          .export_values();

  class_<InteractionSignature, std::shared_ptr<InteractionSignature>>(m, "InteractionSignature")
          .def(init<>())
          .def_readwrite("primary_type",&InteractionSignature::primary_type)
          .def_readwrite("target_type",&InteractionSignature::target_type)
          .def_readwrite("secondary_types",&InteractionSignature::secondary_types);

  class_<InteractionRecord, std::shared_ptr<InteractionRecord>>(m, "InteractionRecord")
          .def(init<>())
          .def_readwrite("signature",&InteractionRecord::signature)
          .def_readwrite("primary_mass",&InteractionRecord::primary_mass)
          .def_readwrite("primary_momentum",&InteractionRecord::primary_momentum)
          .def_readwrite("primary_helicity",&InteractionRecord::primary_helicity)
          .def_readwrite("target_mass",&InteractionRecord::target_mass)
          .def_readwrite("target_momentum",&InteractionRecord::target_momentum)
          .def_readwrite("target_helicity",&InteractionRecord::target_helicity)
          .def_readwrite("interaction_vertex",&InteractionRecord::interaction_vertex)
          .def_readwrite("secondary_masses",&InteractionRecord::secondary_masses)
          .def_readwrite("secondary_momenta",&InteractionRecord::secondary_momenta)
          .def_readwrite("secondary_helicity",&InteractionRecord::secondary_helicity)
          .def_readwrite("interaction_parameters",&InteractionRecord::interaction_parameters);

  class_<InteractionTreeDatum, std::shared_ptr<InteractionTreeDatum>>(m, "InteractionTreeDatum")
          .def(init<InteractionRecord&>())
          .def_readwrite("record",&InteractionTreeDatum::record)
          .def_readwrite("parent",&InteractionTreeDatum::parent)
          .def_readwrite("daughters",&InteractionTreeDatum::daughters)
          .def("depth",&InteractionTreeDatum::depth);

  class_<InteractionTree, std::shared_ptr<InteractionTree>>(m, "InteractionTree")
          .def(init<>())
          .def_readwrite("tree",&InteractionTree::tree)
          .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionTreeDatum&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry))
          .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionRecord&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry));

}
