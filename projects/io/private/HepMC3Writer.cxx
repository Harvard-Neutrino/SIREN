#include "SIREN/io/HepMC3Writer.h"

#include <cmath>
#include <stdexcept>

#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"

#ifdef SIREN_HAS_HEPMC3

#include <HepMC3/Units.h>
#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/WriterAscii.h>

#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace io {

namespace {

namespace C = siren::utilities::Constants;
using siren::dataclasses::ParticleType;

int pdg(ParticleType t) { return static_cast<int32_t>(t); }

// SIREN four-momentum is [E, px, py, pz]; HepMC3 FourVector is (px, py, pz, E).
HepMC3::FourVector Momentum(std::array<double, 4> const & p) {
    return HepMC3::FourVector(p[1], p[2], p[3], p[0]);
}

// SIREN positions are internal meters (Constants::m == 1); GenEvent uses CM.
// The time slot carries c*t re-expressed in the same CM length unit so the
// 4-vector is Minkowski-consistent.
HepMC3::FourVector VertexPosition(std::array<double, 3> const & x, double t_internal) {
    double const s = 1.0 / C::cm; // internal-meter -> CM
    return HepMC3::FourVector(x[0] * s, x[1] * s, x[2] * s, (C::c * t_internal) * s);
}

std::vector<double> PositionCM(std::array<double, 3> const & x) {
    double const s = 1.0 / C::cm;
    return {x[0] * s, x[1] * s, x[2] * s};
}

double SecondsFromInternal(double t_internal) { return t_internal / C::second; }

int VertexStatus(siren::dataclasses::InteractionTreeDatum const & datum) {
    if(datum.is_root()) return 1;                                        // primary vertex
    if(datum.record.signature.target_type == ParticleType::Decay) return 22; // decay
    return 21;                                                           // secondary interaction
}

std::shared_ptr<HepMC3::DoubleAttribute> D(double v) {
    return std::make_shared<HepMC3::DoubleAttribute>(v);
}

} // namespace

struct HepMC3Writer::Impl {
    Impl(std::string const & filename, Options const & options)
        : options_(options) {
        run_info_ = std::make_shared<HepMC3::GenRunInfo>();

        HepMC3::GenRunInfo::ToolInfo tool;
        tool.name = "SIREN";
        tool.version = options_.siren_version;
        tool.description = "SIREN BSM neutrino/dark-sector event generator";
        run_info_->tools().push_back(tool);

        if(options_.weight_names.empty()) options_.weight_names = {"CV"};
        run_info_->set_weight_names(options_.weight_names);

        // NuHepMC version signalling (G.R.1-R.3).
        run_info_->add_attribute("NuHepMC.Version.Major", std::make_shared<HepMC3::IntAttribute>(1));
        run_info_->add_attribute("NuHepMC.Version.Minor", std::make_shared<HepMC3::IntAttribute>(0));
        run_info_->add_attribute("NuHepMC.Version.Patch", std::make_shared<HepMC3::IntAttribute>(0));

        for(auto const & kv : options_.provenance) {
            run_info_->add_attribute("siren." + kv.first,
                                     std::make_shared<HepMC3::StringAttribute>(kv.second));
        }

        writer_ = std::make_unique<HepMC3::WriterAscii>(filename, run_info_);
        if(writer_->failed()) {
            throw std::runtime_error("HepMC3Writer: could not open '" + filename + "' for writing");
        }
    }

    HepMC3::GenEvent TreeToGenEvent(siren::dataclasses::InteractionTree const & tree,
                                    int event_number) const {
        HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::CM);
        evt.set_run_info(run_info_);
        evt.set_event_number(event_number);

        // Central-value weight (and any additional named weights) from the header.
        std::size_t const nweights = run_info_->weight_names().size();
        std::vector<double> weights = tree.header.weights;
        if(weights.size() != nweights) {
            double const cv = weights.empty() ? 1.0 : weights.front();
            weights.assign(nweights, cv);
        }
        evt.weights() = weights;

        // Shared particles: a parent vertex's outgoing secondary is the same
        // physical particle as the daughter vertex's incoming primary; the
        // records' ParticleID linkage makes this an exact join (no momentum
        // matching).
        std::map<siren::dataclasses::ParticleID, HepMC3::GenParticlePtr> particle_by_id;

        for(auto const & datum_ptr : tree.tree) {
            siren::dataclasses::InteractionTreeDatum const & datum = *datum_ptr;
            siren::dataclasses::InteractionRecord const & rec = datum.record;

            HepMC3::GenVertexPtr vertex =
                std::make_shared<HepMC3::GenVertex>(VertexPosition(rec.interaction_vertex, rec.interaction_time));
            vertex->set_status(VertexStatus(datum));

            // Incoming primary particle. A non-root datum's primary is the same
            // physical particle as its parent's outgoing secondary (already
            // created and attributed when the parent vertex was processed).
            HepMC3::GenParticlePtr primary;
            bool primary_is_new = false;
            if(!datum.is_root()) {
                auto it = particle_by_id.find(rec.primary_id);
                if(it != particle_by_id.end()) {
                    primary = it->second;
                    primary->set_status(2); // decayed/re-interacted, no longer final
                }
            }
            if(!primary) {
                primary = std::make_shared<HepMC3::GenParticle>(
                    Momentum(rec.primary_momentum), pdg(rec.signature.primary_type), 4);
                primary->set_generated_mass(rec.primary_mass);
                primary_is_new = true;
            }
            vertex->add_particle_in(primary);

            // Target particle (nucleus/nucleon at rest), skipped for decay vertices.
            HepMC3::GenParticlePtr target;
            ParticleType const target_type = rec.signature.target_type;
            if(target_type != ParticleType::Decay && target_type != ParticleType::unknown) {
                target = std::make_shared<HepMC3::GenParticle>(
                    HepMC3::FourVector(0, 0, 0, rec.target_mass), pdg(target_type), 20);
                target->set_generated_mass(rec.target_mass);
                vertex->add_particle_in(target);
            }

            // Outgoing secondaries.
            std::vector<std::pair<HepMC3::GenParticlePtr, std::size_t>> outgoing;
            std::size_t const n = rec.signature.secondary_types.size();
            for(std::size_t j = 0; j < n; ++j) {
                std::array<double, 4> const mom =
                    (j < rec.secondary_momenta.size()) ? rec.secondary_momenta[j]
                                                       : std::array<double, 4>{0, 0, 0, 0};
                HepMC3::GenParticlePtr out = std::make_shared<HepMC3::GenParticle>(
                    Momentum(mom), pdg(rec.signature.secondary_types[j]), 1);
                if(j < rec.secondary_masses.size()) out->set_generated_mass(rec.secondary_masses[j]);
                if(j < rec.secondary_ids.size()) particle_by_id[rec.secondary_ids[j]] = out;
                vertex->add_particle_out(out);
                outgoing.emplace_back(out, j);
            }

            // Registering the vertex assigns particle/vertex ids; attributes can
            // only be attached once the objects belong to the event.
            evt.add_vertex(vertex);

            if(primary_is_new) {
                primary->add_attribute("siren.helicity", D(rec.primary_helicity));
                // Written as three scalar attributes rather than a vector
                // attribute so the writer compiles against older HepMC3 (3.2.x)
                // that lacks VectorDoubleAttribute.
                std::vector<double> const pos = PositionCM(rec.primary_initial_position);
                primary->add_attribute("siren.primary_initial_position.x", D(pos[0]));
                primary->add_attribute("siren.primary_initial_position.y", D(pos[1]));
                primary->add_attribute("siren.primary_initial_position.z", D(pos[2]));
                primary->add_attribute("siren.primary_initial_time",
                    D(SecondsFromInternal(rec.primary_initial_time)));
            }
            if(target) target->add_attribute("siren.helicity", D(rec.target_helicity));
            for(auto const & item : outgoing) {
                std::size_t const j = item.second;
                if(j < rec.secondary_helicities.size())
                    item.first->add_attribute("siren.helicity", D(rec.secondary_helicities[j]));
                if(j < rec.secondary_times.size())
                    item.first->add_attribute("siren.time", D(SecondsFromInternal(rec.secondary_times[j])));
            }
            // Reweighting-critical opaque map: one DoubleAttribute per key.
            for(auto const & kv : rec.interaction_parameters) {
                vertex->add_attribute("siren.param." + kv.first, D(kv.second));
            }
        }

        return evt;
    }

    void Write(siren::dataclasses::InteractionTree const & tree, int event_number) {
        HepMC3::GenEvent evt = TreeToGenEvent(tree, event_number);
        writer_->write_event(evt);
    }

    void Close() {
        if(writer_) {
            writer_->close();
            writer_.reset();
        }
    }

    Options options_;
    std::shared_ptr<HepMC3::GenRunInfo> run_info_;
    std::unique_ptr<HepMC3::WriterAscii> writer_;
};

HepMC3Writer::HepMC3Writer(std::string const & filename, Options const & options)
    : impl_(std::make_unique<Impl>(filename, options)) {}

HepMC3Writer::~HepMC3Writer() {
    if(impl_) impl_->Close();
}

void HepMC3Writer::Write(siren::dataclasses::InteractionTree const & tree, int event_number) {
    impl_->Write(tree, event_number);
}

void HepMC3Writer::Close() { impl_->Close(); }

} // namespace io
} // namespace siren

#else // SIREN_HAS_HEPMC3 not defined

namespace siren {
namespace io {

struct HepMC3Writer::Impl {};

static void ThrowUnsupported() {
    throw std::runtime_error("SIREN was built without HepMC3 support (SIREN_WITH_HEPMC3=OFF or HepMC3 not found)");
}

HepMC3Writer::HepMC3Writer(std::string const &, Options const &) { ThrowUnsupported(); }
HepMC3Writer::~HepMC3Writer() = default;
void HepMC3Writer::Write(siren::dataclasses::InteractionTree const &, int) { ThrowUnsupported(); }
void HepMC3Writer::Close() {}

} // namespace io
} // namespace siren

#endif // SIREN_HAS_HEPMC3

namespace siren {
namespace io {

// Delegating convenience constructor; works in both the HepMC3 and stub builds.
HepMC3Writer::HepMC3Writer(std::string const & filename)
    : HepMC3Writer(filename, Options()) {}

void SaveInteractionTreesAsHepMC3(
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>> const & trees,
    std::string const & filename,
    HepMC3Writer::Options const & options) {
    HepMC3Writer writer(filename, options);
    int index = 0;
    for(auto const & tree : trees) {
        if(!tree) { ++index; continue; }
        int const event_number =
            (tree->header.event_number != 0) ? static_cast<int>(tree->header.event_number) : index;
        writer.Write(*tree, event_number);
        ++index;
    }
    writer.Close();
}

} // namespace io
} // namespace siren
