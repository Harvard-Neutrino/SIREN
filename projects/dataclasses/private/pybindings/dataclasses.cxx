
#include <vector>
#include <set>

#include "../../public/SIREN/dataclasses/Particle.h"
#include "../../public/SIREN/dataclasses/ParticleID.h"
#include "../../public/SIREN/dataclasses/ParticleType.h"
#include "../../public/SIREN/dataclasses/InteractionSignature.h"
#include "../../public/SIREN/dataclasses/InteractionRecord.h"
#include "../../public/SIREN/dataclasses/InteractionTree.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>


PYBIND11_MODULE(dataclasses, m) {
    namespace py = pybind11;
    using namespace siren::dataclasses;

    // Create a Python class binding for siren::dataclasses::ParticleID
    py::class_<siren::dataclasses::ParticleID>(m, "ParticleID")
        .def(py::init<>())
        .def(py::init<uint64_t, int32_t>(), py::arg("major"), py::arg("minor"))
        .def(py::init<const siren::dataclasses::ParticleID &>(), py::arg("other"))
        .def(py::init<siren::dataclasses::ParticleID>(), py::arg("other"))
        .def_property_readonly("major_id", &siren::dataclasses::ParticleID::GetMajorID)
        .def_property_readonly("minor_id", &siren::dataclasses::ParticleID::GetMinorID)
        .def("__bool__", &siren::dataclasses::ParticleID::operator bool)
        .def("__repr__", [](siren::dataclasses::ParticleID const & id) { return to_repr(id); })
        .def("__str__", [](siren::dataclasses::ParticleID const & id) { return to_str(id); })
        .def("is_set", &siren::dataclasses::ParticleID::IsSet)
        .def("set", &siren::dataclasses::ParticleID::SetID, py::arg("major"), py::arg("minor"))
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def_static("generate_id", &siren::dataclasses::ParticleID::GenerateID)
        ;

    py::class_<Particle, std::shared_ptr<Particle>>(m, "Particle")
        .def(py::init<>())
        .def(py::init<Particle const &>())
        .def(py::init<ParticleID, ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
        .def(py::init<ParticleType, double, std::array<double, 4>, std::array<double, 3>, double, double>())
        .def("__str__", [](Particle const & p) { std::stringstream ss; ss << p; return ss.str(); })
        .def_readwrite("id",&Particle::id)
        .def_readwrite("type",&Particle::type)
        .def_readwrite("mass",&Particle::mass)
        .def_readwrite("momentum",&Particle::momentum)
        .def_readwrite("position",&Particle::position)
        .def_readwrite("length",&Particle::length)
        .def_readwrite("helicity",&Particle::helicity)
        .def("generate_id",&Particle::GenerateID)
        ;

    py::enum_<ParticleType>(m, "ParticleType", py::arithmetic())
#define X(a, b) .value( #a , ParticleType:: a )
#include "../../public/SIREN/dataclasses/ParticleTypes.def"
#undef X
        .export_values();

    py::class_<InteractionSignature, std::shared_ptr<InteractionSignature>>(m, "InteractionSignature")
        .def(py::init<>())
        .def("__str__", [](InteractionSignature const & s) { return to_str(s); })
        .def("__repr__", [](InteractionSignature const & s) { return to_repr(s); })
        .def_readwrite("primary_type",&InteractionSignature::primary_type)
        .def_readwrite("target_type",&InteractionSignature::target_type)
        .def_readwrite("secondary_types",&InteractionSignature::secondary_types)
        ;

    py::class_<PrimaryDistributionRecord, std::shared_ptr<PrimaryDistributionRecord>>(m, "PrimaryDistributionRecord")
        .def(py::init<ParticleType>())
        .def_property_readonly("id",
            [](siren::dataclasses::PrimaryDistributionRecord const & pdr) {siren::dataclasses::ParticleID id = pdr.id; return id;})
        .def_property_readonly("type",
            [](siren::dataclasses::PrimaryDistributionRecord const & pdr) {siren::dataclasses::ParticleType pt = pdr.type; return pt;})
        //.def("GetParticle", &PrimaryDistributionRecord::GetParticle)
        //.def("SetParticle", &PrimaryDistributionRecord::SetParticle)
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
        .def("finalize", &PrimaryDistributionRecord::Finalize);

    py::class_<SecondaryParticleRecord, std::shared_ptr<SecondaryParticleRecord>>(m, "SecondaryParticleRecord")
        .def(py::init<InteractionRecord const &, size_t>())
        .def_property_readonly("id",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {siren::dataclasses::ParticleID id = spr.id; return id;})
        .def_property_readonly("type",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {siren::dataclasses::ParticleType pt = spr.type; return pt;})
        .def_property_readonly("initial_position",
            [](siren::dataclasses::SecondaryParticleRecord const & spr) {std::array<double, 3> ip = spr.initial_position; return ip;})
        //.def("GetParticle", &SecondaryParticleRecord::GetParticle)
        //.def("SetParticle", &SecondaryParticleRecord::SetParticle)
        .def_property("mass", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetMass)), &SecondaryParticleRecord::SetMass)
        .def_property("energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetEnergy)), &SecondaryParticleRecord::SetEnergy)
        .def_property("kinetic_energy", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetKineticEnergy)), &SecondaryParticleRecord::SetKineticEnergy)
        .def_property("direction", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetDirection)), &SecondaryParticleRecord::SetDirection)
        .def_property("three_momentum", ((std::array<double, 3> const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetThreeMomentum)), &SecondaryParticleRecord::SetThreeMomentum)
        .def_property("four_momentum", ((std::array<double, 4> (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetFourMomentum)), &SecondaryParticleRecord::SetFourMomentum)
        .def_property("helicity", ((double const & (SecondaryParticleRecord::*)())(&SecondaryParticleRecord::GetHelicity)), &SecondaryParticleRecord::SetHelicity)
        .def("finalize", &SecondaryParticleRecord::Finalize)
        ;

    py::class_<CrossSectionDistributionRecord, std::shared_ptr<CrossSectionDistributionRecord>>(m, "CrossSectionDistributionRecord")
        .def(py::init<InteractionRecord const &>())
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
        .def_property_readonly("secondary_particle_records",
                [](siren::dataclasses::CrossSectionDistributionRecord & cdr) -> std::vector<siren::dataclasses::SecondaryParticleRecord> & {return cdr.GetSecondaryParticleRecords();},
            py::return_value_policy::reference_internal)
        .def("get_econdary_particle_record",
                [](siren::dataclasses::CrossSectionDistributionRecord & cdr, size_t i) -> siren::dataclasses::SecondaryParticleRecord & {return cdr.GetSecondaryParticleRecord(i);},
                py::return_value_policy::reference_internal)
        .def("get_econdary_particle_records",
                [](siren::dataclasses::CrossSectionDistributionRecord & cdr) -> std::vector<siren::dataclasses::SecondaryParticleRecord> & {return cdr.GetSecondaryParticleRecords();},
                py::return_value_policy::reference_internal)
        .def("finalize", &CrossSectionDistributionRecord::Finalize)
        ;


    py::class_<InteractionRecord, std::shared_ptr<InteractionRecord>>(m, "InteractionRecord")
        .def(py::init<>())
        .def("__str__", [](InteractionRecord const & r) { return to_str(r); })
        .def("__repr__", [](InteractionRecord const & r) { return to_repr(r); })
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
        .def_readwrite("interaction_parameters",&InteractionRecord::interaction_parameters)
        ;

    py::class_<InteractionTreeDatum, std::shared_ptr<InteractionTreeDatum>>(m, "InteractionTreeDatum")
        .def(py::init<InteractionRecord&>())
        .def_readwrite("record",&InteractionTreeDatum::record)
        .def_readwrite("parent",&InteractionTreeDatum::parent)
        .def_readwrite("daughters",&InteractionTreeDatum::daughters)
        .def("depth",&InteractionTreeDatum::depth)
        ;

    py::class_<InteractionTree, std::shared_ptr<InteractionTree>>(m, "InteractionTree")
        .def(py::init<>())
        .def_readwrite("tree",&InteractionTree::tree)
        .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionTreeDatum&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry))
        .def("add_entry",static_cast<std::shared_ptr<InteractionTreeDatum> (InteractionTree::*)(InteractionRecord&,std::shared_ptr<InteractionTreeDatum>)>(&InteractionTree::add_entry))
        ;

    m.def("SaveInteractionTrees",&SaveInteractionTrees);
    m.def("LoadInteractionTrees",&LoadInteractionTrees, py::return_value_policy::reference);

}
