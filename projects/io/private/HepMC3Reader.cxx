#include "SIREN/io/HepMC3Reader.h"

#include <map>
#include <array>
#include <string>
#include <stdexcept>

#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"

#ifdef SIREN_HAS_HEPMC3

#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/ReaderAscii.h>

#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace io {

namespace {

namespace C = siren::utilities::Constants;
using siren::dataclasses::ParticleType;

ParticleType pdg_to_type(int pid) { return static_cast<ParticleType>(pid); }

// HepMC3 FourVector is (px, py, pz, E); SIREN stores [E, px, py, pz].
std::array<double, 4> siren_momentum(HepMC3::FourVector const & p) {
    return {p.e(), p.px(), p.py(), p.pz()};
}

template<typename Ptr>
double double_attr(Ptr const & obj, std::string const & name, double fallback = 0.0) {
    std::shared_ptr<HepMC3::DoubleAttribute> a = obj->template attribute<HepMC3::DoubleAttribute>(name);
    return a ? a->value() : fallback;
}

// The primary is the incoming particle that is not a target. Targets carry
// status 20 (primary interaction) or 22 (deeper cascade interaction).
HepMC3::GenParticlePtr find_primary(HepMC3::GenVertexPtr const & vertex) {
    for(auto const & pin : vertex->particles_in()) {
        if(pin->status() != 20 && pin->status() != 22) return pin;
    }
    return nullptr;
}

siren::dataclasses::InteractionRecord VertexToRecord(HepMC3::GenVertexPtr const & vertex) {
    siren::dataclasses::InteractionRecord rec;

    // Vertex position/time: CM -> internal meters, c*t slot -> internal time.
    HepMC3::FourVector const pos = vertex->position();
    rec.interaction_vertex = {pos.x() * C::cm, pos.y() * C::cm, pos.z() * C::cm};
    rec.interaction_time = (pos.t() * C::cm) / C::c;

    HepMC3::GenParticlePtr const primary = find_primary(vertex);
    if(primary) {
        rec.signature.primary_type = pdg_to_type(primary->pid());
        rec.primary_momentum = siren_momentum(primary->momentum());
        rec.primary_mass = primary->generated_mass();
        rec.primary_helicity = double_attr(primary, "siren.helicity");
        rec.primary_initial_position = {
            double_attr(primary, "siren.primary_initial_position.x") * C::cm,
            double_attr(primary, "siren.primary_initial_position.y") * C::cm,
            double_attr(primary, "siren.primary_initial_position.z") * C::cm};
        rec.primary_initial_time = double_attr(primary, "siren.primary_initial_time") * C::second;
    }

    // Target (status-20 primary or status-22 deeper-interaction incoming), if
    // present; its absence marks a decay vertex.
    HepMC3::GenParticlePtr target;
    for(auto const & pin : vertex->particles_in()) {
        if(pin->status() == 20 || pin->status() == 22) { target = pin; break; }
    }
    if(target) {
        rec.signature.target_type = pdg_to_type(target->pid());
        rec.target_mass = target->generated_mass();
        rec.target_helicity = double_attr(target, "siren.helicity");
    } else {
        rec.signature.target_type = ParticleType::Decay;
    }

    // Outgoing secondaries. ParticleIDs are not stored (not reweighting-relevant)
    // and are left default; tree topology is recovered from the vertex graph.
    for(auto const & out : vertex->particles_out()) {
        rec.signature.secondary_types.push_back(pdg_to_type(out->pid()));
        rec.secondary_momenta.push_back(siren_momentum(out->momentum()));
        rec.secondary_masses.push_back(out->generated_mass());
        rec.secondary_helicities.push_back(double_attr(out, "siren.helicity"));
        rec.secondary_times.push_back(double_attr(out, "siren.time") * C::second);
        rec.secondary_ids.push_back(siren::dataclasses::ParticleID());
    }

    // interaction_parameters from siren.param.<key> attributes on the vertex.
    std::string const prefix = "siren.param.";
    for(std::string const & name : vertex->attribute_names()) {
        if(name.size() > prefix.size() && name.compare(0, prefix.size(), prefix) == 0) {
            rec.interaction_parameters[name.substr(prefix.size())] = double_attr(vertex, name);
        }
    }

    return rec;
}

siren::dataclasses::InteractionTree GenEventToTree(HepMC3::GenEvent & evt) {
    siren::dataclasses::InteractionTree tree;
    tree.header.event_number = static_cast<std::uint64_t>(evt.event_number());
    tree.header.weights = evt.weights();

    std::vector<HepMC3::GenVertexPtr> const vertices = evt.vertices();

    // Reconstruct each vertex's record and its parent vertex id (via the shared
    // incoming primary's production vertex; a beam primary has none -> root).
    std::map<int, siren::dataclasses::InteractionRecord> record_by_vid;
    std::map<int, int> parent_vid; // child vid -> parent vid (0 == no parent)
    for(auto const & vertex : vertices) {
        record_by_vid[vertex->id()] = VertexToRecord(vertex);
        HepMC3::GenParticlePtr const primary = find_primary(vertex);
        HepMC3::GenVertexPtr const parent = primary ? primary->production_vertex() : nullptr;
        parent_vid[vertex->id()] = parent ? parent->id() : 0;
    }

    // Insert in topological order (a node only after its parent), independent of
    // the order vertices() reports, iterating the original order each round so
    // siblings keep a stable order.
    std::map<int, std::shared_ptr<siren::dataclasses::InteractionTreeDatum>> datum_by_vid;
    std::size_t placed = 0;
    bool progress = true;
    while(placed < vertices.size() && progress) {
        progress = false;
        for(auto const & vertex : vertices) {
            int const vid = vertex->id();
            if(datum_by_vid.count(vid)) continue;
            int const pvid = parent_vid[vid];
            bool const is_root = (pvid == 0) || (record_by_vid.find(pvid) == record_by_vid.end());
            if(is_root) {
                datum_by_vid[vid] = tree.add_entry(record_by_vid[vid]);
                ++placed; progress = true;
            } else if(datum_by_vid.count(pvid)) {
                datum_by_vid[vid] = tree.add_entry(record_by_vid[vid], datum_by_vid[pvid]);
                ++placed; progress = true;
            }
        }
    }
    // Defensive: place any node left over by a malformed graph as a root.
    for(auto const & vertex : vertices) {
        if(!datum_by_vid.count(vertex->id())) {
            datum_by_vid[vertex->id()] = tree.add_entry(record_by_vid[vertex->id()]);
        }
    }
    return tree;
}

} // namespace

std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>>
LoadInteractionTreesFromHepMC3(std::string const & filename) {
    HepMC3::ReaderAscii reader(filename);
    if(reader.failed()) {
        throw std::runtime_error("HepMC3Reader: could not open '" + filename + "' for reading");
    }
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>> trees;
    while(true) {
        HepMC3::GenEvent evt;
        reader.read_event(evt);
        if(reader.failed()) break;
        trees.push_back(std::make_shared<siren::dataclasses::InteractionTree>(GenEventToTree(evt)));
    }
    reader.close();
    return trees;
}

} // namespace io
} // namespace siren

#else // SIREN_HAS_HEPMC3 not defined

namespace siren {
namespace io {

std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>>
LoadInteractionTreesFromHepMC3(std::string const &) {
    throw std::runtime_error("SIREN was built without HepMC3 support (SIREN_WITH_HEPMC3=OFF or HepMC3 not found)");
}

} // namespace io
} // namespace siren

#endif // SIREN_HAS_HEPMC3
