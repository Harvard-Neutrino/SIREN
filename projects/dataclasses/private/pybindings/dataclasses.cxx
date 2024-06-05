
#include <vector>
#include <set>

#include "../../public/SIREN/dataclasses/Particle.h"
#include "../../public/SIREN/dataclasses/ParticleID.h"
#include "../../public/SIREN/dataclasses/ParticleType.h"
#include "../../public/SIREN/dataclasses/InteractionSignature.h"
#include "../../public/SIREN/dataclasses/InteractionRecord.h"
#include "../../public/SIREN/dataclasses/InteractionTree.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace pybind11;

PYBIND11_MODULE(dataclasses,m) {
  using namespace siren::dataclasses;

  class_<Particle, std::shared_ptr<Particle>> particle(m, "Particle");

  particle.def(init<>())
          .def(init<Particle const &>())
          .def(init<ParticleID, ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
          .def(init<ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
          .def("__str__", [](Particle const & p) { std::stringstream ss; ss << p; return ss.str(); })
          .def_readwrite("id",&Particle::id)
          .def_readwrite("type",&Particle::type)
          .def_readwrite("mass",&Particle::mass)
          .def_readwrite("momentum",&Particle::momentum)
          .def_readwrite("position",&Particle::position)
          .def_readwrite("length",&Particle::length)
          .def_readwrite("helicity",&Particle::helicity)
          .def("GenerateID",&Particle::GenerateID);

    enum_<ParticleType>(particle, "ParticleType", arithmetic())
#define X(a, b) .value( #a , ParticleType:: a )
#include "../../public/SIREN/dataclasses/ParticleTypes.def"
#undef X
        .export_values();

    class_<InteractionSignature, std::shared_ptr<InteractionSignature>>(m, "InteractionSignature")
        .def(init<>())
        .def("__str__", [](InteractionSignature const & p) { std::stringstream ss; ss << p; return ss.str(); })
        .def("__repr__", [](InteractionSignature const & s) {
            std::stringstream ss;
            ss << "InteractionSignature( ";
            ss << s.primary_type << " ";
            if(s.primary_type == ParticleType::unknown or s.target_type != ParticleType::unknown) {
                ss << s.target_type << " ";
            }
            ss << "-> ";
            for(auto const & secondary : s.secondary_types) {
                ss << secondary << " ";
            }
            ss << ")";
            return ss.str();
        })
        .def_readwrite("primary_type",&InteractionSignature::primary_type)
        .def_readwrite("target_type",&InteractionSignature::target_type)
        .def_readwrite("secondary_types",&InteractionSignature::secondary_types);

    class_<PrimaryDistributionRecord, std::shared_ptr<PrimaryDistributionRecord>>(m, "PrimaryDistributionRecord")
        .def(init<ParticleType>())
        .def_property_readonly("id",
            [](siren::dataclasses::PrimaryDistributionRecord const & pdr) {siren::dataclasses::ParticleID id = pdr.id; return id;})
        .def_property_readonly("type",
            [](siren::dataclasses::PrimaryDistributionRecord const & pdr) {siren::dataclasses::ParticleType pt = pdr.type; return pt;})
        .def("GetParticle", &PrimaryDistributionRecord::GetParticle)
        .def("SetParticle", &PrimaryDistributionRecord::SetParticle)
        .def_property("mass", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetMass)), &PrimaryDistributionRecord::SetMass)
        .def_property("energy", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetEnergy)), &PrimaryDistributionRecord::SetEnergy)
        .def_property("kinetic_energy", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetKineticEnergy)), &PrimaryDistributionRecord::SetKineticEnergy)
        .def_property("direction", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetDirection)), &PrimaryDistributionRecord::SetDirection)
        .def_property("three_momentum", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetThreeMomentum)), &PrimaryDistributionRecord::SetThreeMomentum)
        .def_property("four_momentum", ((std::array<double, 4> (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetFourMomentum)), &PrimaryDistributionRecord::SetFourMomentum)
        .def_property("length", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetLength)), &PrimaryDistributionRecord::SetLength)
        .def_property("initial_position", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetInitialPosition)), &PrimaryDistributionRecord::SetInitialPosition)
        .def_property("interaction_vertex", ((std::array<double, 3> const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetInteractionVertex)), &PrimaryDistributionRecord::SetInteractionVertex)
        .def_property("helicity", ((double const & (PrimaryDistributionRecord::*)())(&PrimaryDistributionRecord::GetHelicity)), &PrimaryDistributionRecord::SetHelicity)
        .def("Finalize", &PrimaryDistributionRecord::Finalize);

    class_<SecondaryParticleRecord, std::shared_ptr<SecondaryParticleRecord>>(m, "SecondaryParticleRecord")
        .def(init<InteractionRecord const &, size_t>())
        .def_property_readonly("id",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {siren::dataclasses::ParticleID id = spr.id; return id;})
        .def_property_readonly("type",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {siren::dataclasses::ParticleType pt = spr.type; return pt;})
        .def_property_readonly("initial_position",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {std::array<double, 3> ip = spr.initial_position; return ip;})
        .def("GetParticle", &SecondaryParticleRecord::GetParticle)
        .def("SetParticle", &SecondaryParticleRecord::SetParticle)
        .def_property("mass", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetMass)), &SecondaryParticleRecord::SetMass)
        .def_property("energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetEnergy)), &SecondaryParticleRecord::SetEnergy)
        .def_property("kinetic_energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetKineticEnergy)), &SecondaryParticleRecord::SetKineticEnergy)
        .def_property("direction", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetDirection)), &SecondaryParticleRecord::SetDirection)
        .def_property("three_momentum", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetThreeMomentum)), &SecondaryParticleRecord::SetThreeMomentum)
        .def_property("four_momentum", ((std::array<double, 4> (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetFourMomentum)), &SecondaryParticleRecord::SetFourMomentum)
        .def_property("helicity", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetHelicity)), &SecondaryParticleRecord::SetHelicity)
        .def("Finalize", &SecondaryParticleRecord::Finalize);

    class_<CrossSectionDistributionRecord, std::shared_ptr<CrossSectionDistributionRecord>>(m, "CrossSectionDistributionRecord")
        .def(init<InteractionRecord const &>())
        .def_property_readonly("record",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::InteractionRecord ir = cdr.record; return ir;})
        .def_property_readonly("signature",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::InteractionSignature is = cdr.signature; return is;})
        .def_property_readonly("primary_id",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::ParticleID id = cdr.primary_id; return id;})
        .def_property_readonly("primary_type",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::ParticleType pt = cdr.primary_type; return pt;})
        .def_property_readonly("primary_initial_position",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {std::array<double, 3> ip = cdr.primary_initial_position; return ip;})
        .def_property_readonly("primary_mass",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {double m = cdr.primary_mass; return m;})
        .def_property_readonly("primary_momentum",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {std::array<double, 4> p = cdr.primary_momentum; return p;})
        .def_property_readonly("primary_helicity",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {double h = cdr.primary_helicity; return h;})
        .def_property_readonly("interaction_vertex",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {std::array<double, 3> iv = cdr.interaction_vertex; return iv;})
        .def_property_readonly("target_id",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::ParticleID id = cdr.target_id; return id;})
        .def_property_readonly("target_type",
            [](siren::dataclasses::CrossSectionDistributionRecord const & cdr) {siren::dataclasses::ParticleType pt = cdr.target_type; return pt;})
        .def_property("target_mass", ((double const & (siren::dataclasses::CrossSectionDistributionRecord::*)() const)(&siren::dataclasses::CrossSectionDistributionRecord::GetTargetMass)), &siren::dataclasses::CrossSectionDistributionRecord::SetTargetMass)
        .def_property("target_helicity", ((double const & (siren::dataclasses::CrossSectionDistributionRecord::*)() const)(&siren::dataclasses::CrossSectionDistributionRecord::GetTargetHelicity)), &siren::dataclasses::CrossSectionDistributionRecord::SetTargetHelicity)
        .def_property("interaction_parameters", ((std::map<std::string, double> const & (siren::dataclasses::CrossSectionDistributionRecord::*)())(&siren::dataclasses::CrossSectionDistributionRecord::GetInteractionParameters)), &siren::dataclasses::CrossSectionDistributionRecord::SetInteractionParameters)
        .def("GetSecondaryParticleRecord",
                [](siren::dataclasses::CrossSectionDistributionRecord & cdr, size_t i) -> siren::dataclasses::SecondaryParticleRecord & {return cdr.GetSecondaryParticleRecord(i);},
                return_value_policy::reference_internal)
        .def("GetSecondaryParticleRecords",
                [](siren::dataclasses::CrossSectionDistributionRecord & cdr) -> std::vector<siren::dataclasses::SecondaryParticleRecord> & {return cdr.GetSecondaryParticleRecords();},
                return_value_policy::reference_internal)
        .def("Finalize", &CrossSectionDistributionRecord::Finalize);


  class_<InteractionRecord, std::shared_ptr<InteractionRecord>>(m, "InteractionRecord")
          .def(init<>())
          .def("__str__", [](InteractionRecord const & r) { std::stringstream ss; ss << r; return ss.str(); })
          .def_readwrite("signature",&InteractionRecord::signature)
          .def_readwrite("primary_mass",&InteractionRecord::primary_mass)
          .def_readwrite("primary_momentum",&InteractionRecord::primary_momentum)
          .def_readwrite("primary_helicity",&InteractionRecord::primary_helicity)
          .def_readwrite("target_mass",&InteractionRecord::target_mass)
          .def_readwrite("target_helicity",&InteractionRecord::target_helicity)
          .def_readwrite("interaction_vertex",&InteractionRecord::interaction_vertex)
          .def_readwrite("secondary_masses",&InteractionRecord::secondary_masses)
          .def_readwrite("secondary_momenta",&InteractionRecord::secondary_momenta)
          .def_readwrite("secondary_helicities",&InteractionRecord::secondary_helicities)
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

  m.def("SaveInteractionTrees",&SaveInteractionTrees);
  m.def("LoadInteractionTrees",&LoadInteractionTrees,pybind11::return_value_policy::reference);

}
