#include "SIREN/io/HepMC3Reader.h"

#include <map>
#include <array>
#include <string>
#include <memory>
#include <fstream>
#include <stdexcept>

#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"

#ifdef SIREN_HAS_HEPMC3

#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/Units.h>
#ifdef SIREN_HEPMC3_HAS_COMPRESSION
#include <HepMC3/ReaderGZ.h>
#endif

#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace io {

namespace {

namespace C = siren::utilities::Constants;
using siren::dataclasses::ParticleType;

ParticleType pdg_to_type(int pid) { return static_cast<ParticleType>(pid); }

// Detect gzip by content (magic bytes 0x1f 0x8b) rather than by filename, so a
// gzip file without a .gz suffix still round-trips and a missing/plaintext file
// falls through to the Ascii reader (which reports the open/parse error).
bool looks_gzip(std::string const & filename) {
    std::ifstream f(filename, std::ios::binary);
    if(!f) return false;
    char b0 = 0, b1 = 0;
    f.get(b0);
    f.get(b1);
    return static_cast<unsigned char>(b0) == 0x1f && static_cast<unsigned char>(b1) == 0x8b;
}

// HepMC3 FourVector is (px, py, pz, E); SIREN stores [E, px, py, pz].
std::array<double, 4> siren_momentum(HepMC3::FourVector const & p) {
    return {p.e(), p.px(), p.py(), p.pz()};
}

// Optional attribute: returns the value if present and typed as a DoubleAttribute,
// otherwise the fallback. Used for the enumerate-if-present namespaces
// (siren.param.*, siren.time) and for every attribute in non-strict mode.
template<typename Ptr>
double double_attr(Ptr const & obj, std::string const & name, double fallback = 0.0) {
    std::shared_ptr<HepMC3::DoubleAttribute> a = obj->template attribute<HepMC3::DoubleAttribute>(name);
    return a ? a->value() : fallback;
}

// Required attribute: a SIREN file the writer produced always carries this key as
// a DoubleAttribute, so a missing or mistyped one signals a corrupt/foreign file.
// In strict mode that throws with the attribute name, the object, and the names of
// the attributes actually present; in lax mode it degrades to the 0.0 fallback.
template<typename Ptr>
double required_double_attr(Ptr const & obj, std::string const & name,
                            std::string const & object, bool strict) {
    std::shared_ptr<HepMC3::DoubleAttribute> a = obj->template attribute<HepMC3::DoubleAttribute>(name);
    if(a) return a->value();
    if(!strict) return 0.0;
    std::string present;
    for(std::string const & n : obj->attribute_names()) {
        if(!present.empty()) present += ", ";
        present += n;
    }
    throw std::runtime_error(
        "HepMC3Reader: required attribute '" + name + "' missing or not a DoubleAttribute on "
        + object + "; attributes present: [" + present + "]");
}

// The primary is the incoming particle that is not a target. Targets carry
// status 20 (primary interaction) or 22 (deeper cascade interaction).
HepMC3::GenParticlePtr find_primary(HepMC3::GenVertexPtr const & vertex) {
    for(auto const & pin : vertex->particles_in()) {
        if(pin->status() != 20 && pin->status() != 22) return pin;
    }
    return nullptr;
}

// Assigns fresh ParticleIDs to reconstructed records so a re-exported tree does
// not collapse topology (HepMC3Writer only joins particles across vertices on a
// set ParticleID). HepMC3 shares one GenParticlePtr between a parent vertex's
// outgoing secondary and the daughter vertex's incoming primary, so keying this
// map on the raw pointer identity gives the same ParticleID to both sides of
// that join, exactly mirroring the writer's shared-particle contract.
siren::dataclasses::ParticleID GetOrAssignId(
        HepMC3::GenParticlePtr const & particle,
        std::map<HepMC3::GenParticlePtr, siren::dataclasses::ParticleID> & ids) {
    auto it = ids.find(particle);
    if(it != ids.end()) return it->second;
    siren::dataclasses::ParticleID const id = siren::dataclasses::ParticleID::GenerateID();
    ids[particle] = id;
    return id;
}

siren::dataclasses::InteractionRecord VertexToRecord(
        HepMC3::GenVertexPtr const & vertex, bool strict,
        std::map<HepMC3::GenParticlePtr, siren::dataclasses::ParticleID> & particle_ids) {
    siren::dataclasses::InteractionRecord rec;

    // Vertex position/time: CM -> internal meters, c*t slot -> internal time.
    HepMC3::FourVector const pos = vertex->position();
    rec.interaction_vertex = {pos.x() * C::cm, pos.y() * C::cm, pos.z() * C::cm};
    rec.interaction_time = (pos.t() * C::cm) / C::c;

    HepMC3::GenParticlePtr const primary = find_primary(vertex);
    if(primary) {
        rec.primary_id = GetOrAssignId(primary, particle_ids);
        rec.signature.primary_type = pdg_to_type(primary->pid());
        rec.primary_momentum = siren_momentum(primary->momentum());
        rec.primary_mass = primary->generated_mass();
        // Unconditionally emitted by the writer -> required.
        rec.primary_helicity = required_double_attr(primary, "siren.helicity", "primary particle", strict);
        rec.primary_initial_position = {
            required_double_attr(primary, "siren.primary_initial_position.x", "primary particle", strict) * C::cm,
            required_double_attr(primary, "siren.primary_initial_position.y", "primary particle", strict) * C::cm,
            required_double_attr(primary, "siren.primary_initial_position.z", "primary particle", strict) * C::cm};
        rec.primary_initial_time =
            required_double_attr(primary, "siren.primary_initial_time", "primary particle", strict) * C::second;
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
        rec.target_helicity = required_double_attr(target, "siren.helicity", "target particle", strict);
    } else {
        // No target particle: either a genuine decay, or an unknown-target
        // interaction (the writer omits the target particle for both). The
        // siren.target_type_unknown vertex flag disambiguates the latter.
        std::shared_ptr<HepMC3::IntAttribute> const unknown_flag =
            vertex->attribute<HepMC3::IntAttribute>("siren.target_type_unknown");
        rec.signature.target_type = (unknown_flag && unknown_flag->value() != 0)
            ? ParticleType::unknown : ParticleType::Decay;
    }

    // Outgoing secondaries. Each gets a fresh ParticleID (or the id already
    // assigned to it as another vertex's primary, if this event was itself
    // re-exported) so a daughter vertex's incoming primary -- the same
    // GenParticlePtr -- resolves to the identical id and re-export can rejoin
    // the topology.
    for(auto const & out : vertex->particles_out()) {
        rec.signature.secondary_types.push_back(pdg_to_type(out->pid()));
        rec.secondary_momenta.push_back(siren_momentum(out->momentum()));
        rec.secondary_masses.push_back(out->generated_mass());
        // Helicity is unconditional (required); siren.time is only emitted when the
        // record carried a secondary time, so it stays enumerate-if-present.
        rec.secondary_helicities.push_back(required_double_attr(out, "siren.helicity", "secondary particle", strict));
        rec.secondary_times.push_back(double_attr(out, "siren.time") * C::second);
        rec.secondary_ids.push_back(GetOrAssignId(out, particle_ids));
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

siren::dataclasses::InteractionTree GenEventToTree(HepMC3::GenEvent & evt, bool strict) {
    siren::dataclasses::InteractionTree tree;
    // The GenEvent stores the event number narrowed to int; the writer preserves
    // the full 64-bit identity in a per-event siren.event_number ULongAttribute
    // when the value overflowed int, so prefer it when present.
    tree.header.event_number = static_cast<std::uint64_t>(evt.event_number());
    {
        std::shared_ptr<HepMC3::ULongAttribute> a =
            evt.attribute<HepMC3::ULongAttribute>("siren.event_number");
        if(a) tree.header.event_number = static_cast<std::uint64_t>(a->value());
    }
    tree.header.weights = evt.weights();

    std::vector<HepMC3::GenVertexPtr> const vertices = evt.vertices();

    // Reconstruct each vertex's record and its parent vertex id (via the shared
    // incoming primary's production vertex; a beam primary has none -> root).
    // particle_ids is shared across every vertex in this event so the same
    // GenParticlePtr always resolves to the same ParticleID.
    std::map<int, siren::dataclasses::InteractionRecord> record_by_vid;
    std::map<int, int> parent_vid; // child vid -> parent vid (0 == no parent)
    std::map<HepMC3::GenParticlePtr, siren::dataclasses::ParticleID> particle_ids;
    for(auto const & vertex : vertices) {
        record_by_vid[vertex->id()] = VertexToRecord(vertex, strict, particle_ids);
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
LoadInteractionTreesFromHepMC3(std::string const & filename, bool strict) {
    // Auto-detect gzip by content (magic bytes); deduce_reader is avoided to stay off
    // HepMC3's optional plugin/root/protobuf paths.
    std::unique_ptr<HepMC3::Reader> reader;
    if(looks_gzip(filename)) {
#ifdef SIREN_HEPMC3_HAS_COMPRESSION
        reader = std::make_unique<HepMC3::ReaderGZ<HepMC3::ReaderAscii>>(filename);
#else
        throw std::runtime_error("HepMC3Reader: '" + filename + "' looks gzip-compressed "
                                 "but this SIREN build's HepMC3 has no compression support");
#endif
    } else {
        reader = std::make_unique<HepMC3::ReaderAscii>(filename);
    }
    if(reader->failed()) {
        throw std::runtime_error("HepMC3Reader: could not open '" + filename + "' for reading");
    }
    {
        std::ifstream probe(filename, std::ios::binary | std::ios::ate);
        if(probe && probe.tellg() == 0) {
            throw std::runtime_error("HepMC3Reader: '" + filename + "' is empty");
        }
    }
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>> trees;
    while(true) {
        HepMC3::GenEvent evt;
        bool const read_ok = reader->read_event(evt);
        // ReaderAscii/ReaderGZ return false only on a genuine mid-file parse
        // failure (and set the stream's bad bit); a clean end-of-file instead
        // returns true with an event carrying no vertices, with the stream's eof
        // bit set. Distinguishing the two means a truncated/corrupt file is
        // reported instead of silently yielding a partial tree list.
        if(!read_ok) {
            throw std::runtime_error("HepMC3Reader: failed to parse an event in '" + filename + "'");
        }
        if(reader->failed()) {
            if(!evt.vertices().empty()) {
                // Successfully parsed the final event; the eof bit merely
                // reflects that no further event follows. Fall through and
                // process it normally, then stop on the next iteration.
            } else {
                break; // clean end of file
            }
        }
        HepMC3::Units::MomentumUnit const mom_unit = evt.momentum_unit();
        HepMC3::Units::LengthUnit const len_unit = evt.length_unit();
        if(mom_unit != HepMC3::Units::GEV || len_unit != HepMC3::Units::CM) {
            evt.set_units(HepMC3::Units::GEV, HepMC3::Units::CM);
        }
        auto tree = std::make_shared<siren::dataclasses::InteractionTree>(GenEventToTree(evt, strict));
        // Restore run-level weight provenance: the writer records how the per-event
        // CV weights were produced as a GenRunInfo siren.weights_state
        // StringAttribute; echo it into each tree header's provenance map when
        // present so a re-weighting/pooling pass can honor the original policy.
        std::shared_ptr<HepMC3::GenRunInfo> const ri = evt.run_info();
        if(ri) {
            std::shared_ptr<HepMC3::StringAttribute> ws =
                ri->attribute<HepMC3::StringAttribute>("siren.weights_state");
            if(ws) tree->header.provenance["siren.weights_state"] = ws->value();
        }
        trees.push_back(std::move(tree));
    }
    reader->close();
    return trees;
}

} // namespace io
} // namespace siren

#else // SIREN_HAS_HEPMC3 not defined

namespace siren {
namespace io {

std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>>
LoadInteractionTreesFromHepMC3(std::string const &, bool) {
    throw std::runtime_error("SIREN was built without HepMC3 support (SIREN_WITH_HEPMC3=OFF or HepMC3 not found)");
}

} // namespace io
} // namespace siren

#endif // SIREN_HAS_HEPMC3
