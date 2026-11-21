#include "SIREN/io/HepMC3Writer.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <set>
#include <stdexcept>

#include "SIREN/dataclasses/InteractionTree.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/ParticleType.h"

// Signature-derived process identity (NuHepMC G.R.8). No HepMC3 types, so these
// live outside the SIREN_HAS_HEPMC3 guard for use by the pre-scan in both builds.
namespace siren {
namespace io {
namespace {

int pdg(siren::dataclasses::ParticleType t) { return static_cast<int32_t>(t); }

std::string ProcessKey(siren::dataclasses::InteractionSignature const & s) {
    std::string k = std::to_string(pdg(s.primary_type)) + "|"
                  + std::to_string(pdg(s.target_type)) + ">";
    for(auto const & st : s.secondary_types) k += std::to_string(pdg(st)) + ",";
    return k;
}

std::string ProcessName(siren::dataclasses::InteractionSignature const & s) {
    std::string n = std::to_string(pdg(s.primary_type)) + "+"
                  + std::to_string(pdg(s.target_type)) + " -> ";
    for(std::size_t j = 0; j < s.secondary_types.size(); ++j) {
        if(j) n += ",";
        n += std::to_string(pdg(s.secondary_types[j]));
    }
    return n;
}

} // namespace
} // namespace io
} // namespace siren

#ifdef SIREN_HAS_HEPMC3

#include <HepMC3/Units.h>
#include <HepMC3/Attribute.h>
#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAscii.h>
#ifdef SIREN_HEPMC3_HAS_COMPRESSION
#include <HepMC3/WriterGZ.h>
#endif

#include "SIREN/dataclasses/ParticleID.h"
#include "SIREN/dataclasses/ParticleType.h"
#include "SIREN/utilities/Constants.h"

namespace siren {
namespace io {

namespace {

namespace C = siren::utilities::Constants;
using siren::dataclasses::ParticleType;

// 1 GeV^-2 = (hbar c)^2 = 0.3894 mb = 3.894e8 pb (PDG).
constexpr double kGeVm2_to_pb = 3.894e8;

// SIREN four-momentum is [E, px, py, pz]; HepMC3 FourVector is (px, py, pz, E).
HepMC3::FourVector Momentum(std::array<double, 4> const & p) {
    return HepMC3::FourVector(p[1], p[2], p[3], p[0]);
}

// SIREN positions are internal meters (Constants::m == 1); GenEvent uses CM.
// Time slot carries c*t in the same CM length unit (Minkowski-consistent).
HepMC3::FourVector VertexPosition(std::array<double, 3> const & x, double t_internal) {
    double const s = 1.0 / C::cm; // internal-meter -> CM
    return HepMC3::FourVector(x[0] * s, x[1] * s, x[2] * s, (C::c * t_internal) * s);
}

std::vector<double> PositionCM(std::array<double, 3> const & x) {
    double const s = 1.0 / C::cm;
    return {x[0] * s, x[1] * s, x[2] * s};
}

double SecondsFromInternal(double t_internal) { return t_internal / C::second; }

// NuHepMC V.R.1 allows only one status-1 primary vertex; other root nodes are
// demoted to status 21 (secondary interaction).
int VertexStatus(siren::dataclasses::InteractionTreeDatum const & datum, bool is_primary_root) {
    if(is_primary_root) return 1;                                        // primary vertex
    if(datum.record.signature.target_type == ParticleType::Decay) return 22; // decay
    return 21;                                                           // secondary interaction
}

std::shared_ptr<HepMC3::DoubleAttribute> D(double v) {
    return std::make_shared<HepMC3::DoubleAttribute>(v);
}

// Write a set ParticleID as two SIREN-namespaced attributes. major_id is
// uint64 (ULongAttribute) and minor_id is int32 (IntAttribute). An unset id
// writes nothing
void WriteParticleID(HepMC3::GenParticlePtr const & particle,
                     siren::dataclasses::ParticleID const & id) {
    if(!id.IsSet()) return;
    particle->add_attribute("siren.id.major",
        std::make_shared<HepMC3::ULongAttribute>(static_cast<unsigned long>(id.GetMajorID())));
    particle->add_attribute("siren.id.minor",
        std::make_shared<HepMC3::IntAttribute>(id.GetMinorID()));
}

// Declare a NuHepMC ID registry on the run info: a VectorIntAttribute id-list
// under list_key, plus a Name and Description StringAttribute per id under
// info_stub + "[<id>]". Strict for attribute spelling and types.
void DeclareIdRegistry(HepMC3::GenRunInfo & ri,
                       std::string const & list_key,
                       std::string const & info_stub,
                       std::vector<int> const & ids,
                       std::map<int, std::pair<std::string, std::string>> const & info) {
    ri.add_attribute(list_key, std::make_shared<HepMC3::VectorIntAttribute>(ids));
    for(int id : ids) {
        auto it = info.find(id);
        std::string const name = (it != info.end()) ? it->second.first : std::string();
        std::string const desc = (it != info.end()) ? it->second.second : std::string();
        std::string const suffix = "[" + std::to_string(id) + "]";
        ri.add_attribute(info_stub + suffix + ".Name",
                         std::make_shared<HepMC3::StringAttribute>(name));
        ri.add_attribute(info_stub + suffix + ".Description",
                         std::make_shared<HepMC3::StringAttribute>(desc));
    }
}

// NuHepMC G.R.11 additional particle numbers. Reference writer/reader disagree
// on the per-code stub (AdditionalParticleNumber vs AdditionalParticleInfo), so
// Name is emitted under both and Description (always written, empty allowed)
// under Info.
void DeclareAdditionalParticleNumbers(
        HepMC3::GenRunInfo & ri,
        std::map<int, std::pair<std::string, std::string>> const & particles) {
    std::vector<int> codes;
    codes.reserve(particles.size());
    for(auto const & kv : particles) codes.push_back(kv.first);
    ri.add_attribute("NuHepMC.AdditionalParticleNumbers",
                     std::make_shared<HepMC3::VectorIntAttribute>(codes));
    for(auto const & kv : particles) {
        std::string const suffix = "[" + std::to_string(kv.first) + "]";
        ri.add_attribute("NuHepMC.AdditionalParticleNumber" + suffix + ".Name",
                         std::make_shared<HepMC3::StringAttribute>(kv.second.first));
        ri.add_attribute("NuHepMC.AdditionalParticleInfo" + suffix + ".Name",
                         std::make_shared<HepMC3::StringAttribute>(kv.second.first));
        ri.add_attribute("NuHepMC.AdditionalParticleInfo" + suffix + ".Description",
                         std::make_shared<HepMC3::StringAttribute>(kv.second.second));
    }
}

// SIREN's non-PDG particle codes (BSM neutrinos, HNLs, dark-sector states) that
// reuse the 59xx block; declared so external readers can interpret them.
std::map<int, std::pair<std::string, std::string>> SirenAdditionalParticleNumbers() {
    return {
        {5901,  {"HPrime",     "SIREN dark scalar H'"}},
        {5902,  {"PhiPrime",   "SIREN dark scalar phi'"}},
        {5903,  {"ALP",        "SIREN axion-like particle"}},
        {5910,  {"NuLight",    "SIREN light neutrino (flavor-summed)"}},
        {-5910, {"NuLightBar", "SIREN light antineutrino (flavor-summed)"}},
        {5911,  {"Nu1",        "SIREN neutrino mass eigenstate 1"}},
        {5912,  {"Nu2",        "SIREN neutrino mass eigenstate 2"}},
        {5913,  {"Nu3",        "SIREN neutrino mass eigenstate 3"}},
        {5914,  {"N4",         "SIREN heavy neutral lepton N4"}},
        {-5914, {"N4Bar",      "SIREN heavy neutral lepton N4 (antiparticle)"}},
        {5915,  {"N5",         "SIREN heavy neutral lepton N5"}},
        {-5915, {"N5Bar",      "SIREN heavy neutral lepton N5 (antiparticle)"}},
        {5916,  {"N6",         "SIREN heavy neutral lepton N6"}},
        {-5916, {"N6Bar",      "SIREN heavy neutral lepton N6 (antiparticle)"}},
        {5921,  {"ZPrime",     "SIREN dark vector Z'"}},
    };
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

        // "unweighted" turns off every NuHepMC.* key (per-event CV weight is a
        // 1.0 placeholder); output is then a plain HepMC3 file with siren.* only.
        run_info_->add_attribute("siren.weights_state",
            std::make_shared<HepMC3::StringAttribute>(options_.weights_state));
        bool const nuhepmc_mode = (options_.weights_state != "unweighted");

        // Version's mandatory backing (G.R.5 FATX + G.R.6 units + G.R.8 process
        // registry) needs at least one accepted event. This guard is reused by
        // the FATX block and the G.R.8 registry below.
        bool const have_events = options_.accepted_events > 0;
        bool const nuhepmc_conformant = nuhepmc_mode && have_events;

        // NuHepMC version signalling (G.R.1-R.3).
        if(nuhepmc_conformant) {
            run_info_->add_attribute("NuHepMC.Version.Major", std::make_shared<HepMC3::IntAttribute>(1));
            run_info_->add_attribute("NuHepMC.Version.Minor", std::make_shared<HepMC3::IntAttribute>(0));
            run_info_->add_attribute("NuHepMC.Version.Patch", std::make_shared<HepMC3::IntAttribute>(0));
        }

        for(auto const & kv : options_.provenance) {
            run_info_->add_attribute("siren." + kv.first,
                                     std::make_shared<HepMC3::StringAttribute>(kv.second));
        }

        if(nuhepmc_mode) {
            // NuHepMC vertex status registry (G.R.9). 1 is the standard primary code;
            // 21/22 are in the generator-dependent band and so must be declared.
            DeclareIdRegistry(*run_info_, "NuHepMC.VertexStatusIDs", "NuHepMC.VertexStatusInfo",
                {1, 21, 22},
                {{1,  {"PrimaryInteraction", "Primary interaction vertex where the beam particle interacts"}},
                 {21, {"SecondaryInteraction", "Interaction vertex of a secondary (cascade) particle"}},
                 {22, {"Decay", "Vertex where a particle decays"}}});

            // NuHepMC particle status registry (G.R.10). 1/2/4/20 are standard codes;
            // 22 (deeper-interaction target) is generator-dependent and must be declared.
            DeclareIdRegistry(*run_info_, "NuHepMC.ParticleStatusIDs", "NuHepMC.ParticleStatusInfo",
                {1, 2, 4, 20, 22},
                {{1,  {"FinalState", "Undecayed final-state particle"}},
                 {2,  {"DecayedOrReinteracted", "Particle that decayed or re-interacted downstream"}},
                 {4,  {"IncomingBeam", "Incoming beam particle of the primary interaction"}},
                 {20, {"Target", "Target of the primary interaction"}},
                 {22, {"NonTargetInteractionTarget", "Target of a deeper (non-primary) interaction in the cascade"}}});

            // NuHepMC additional particle numbers (G.R.11): SIREN's non-PDG codes,
            // plus any the caller supplies.
            std::map<int, std::pair<std::string, std::string>> particles = SirenAdditionalParticleNumbers();
            for(auto const & kv : options_.additional_particle_numbers) particles[kv.first] = kv.second;
            DeclareAdditionalParticleNumbers(*run_info_, particles);
        }

        // Run-level generation counts (also the FATX normalization). DoubleAttribute
        // avoids int overflow for attempts > INT_MAX; exact up to 2^53.
        if(options_.attempted_events >= 0)
            run_info_->add_attribute("siren.attempted_events",
                D(static_cast<double>(options_.attempted_events)));
        if(options_.accepted_events >= 0)
            run_info_->add_attribute("siren.accepted_events",
                D(static_cast<double>(options_.accepted_events)));

        // Pooled-weighting seed N_i (Injector EventsToInject). Metadata only.
        if(options_.events_to_inject >= 0)
            run_info_->add_attribute("siren.events_to_inject",
                D(static_cast<double>(options_.events_to_inject)));

        // NuHepMC process ID registry (G.R.8): one id per distinct root interaction
        // signature, assigned in the generator ("Other", >= 700) band by the pre-scan.
        if(nuhepmc_conformant && !options_.process_names.empty()) {
            std::vector<int> ids;
            std::map<int, std::pair<std::string, std::string>> info;
            for(auto const & kv : options_.process_names) {
                ids.push_back(kv.first);
                info[kv.first] = {kv.second, "SIREN interaction " + kv.second};
            }
            DeclareIdRegistry(*run_info_, "NuHepMC.ProcessIDs", "NuHepMC.ProcessInfo", ids, info);
        }

        // Flux-averaged total cross section, at run level (NuHepMC G.C.2, not
        // E.C.4: a single run-wide constant, not a per-event running estimate).
        // Each CV weight already carries 1/EventsToInject, so sum(CV) is the
        // unbiased estimator of the integral; only emitted when nuhepmc_conformant.
        if(nuhepmc_conformant) {
            run_info_->add_attribute("siren.fatx.weight_sum", D(options_.fatx_weight_sum));
            double const fatx = kGeVm2_to_pb * options_.fatx_weight_sum;
            run_info_->add_attribute("siren.fatx.value", D(fatx));
            run_info_->add_attribute("NuHepMC.Units.CrossSection.Unit",
                std::make_shared<HepMC3::StringAttribute>(options_.cross_section_unit));
            // TargetScale is only in-spec as "PerAtom" or "PerNucleon" (G.R.6);
            // omitted when fatx_per_atom is false.
            if(options_.fatx_per_atom) {
                run_info_->add_attribute("NuHepMC.Units.CrossSection.TargetScale",
                    std::make_shared<HepMC3::StringAttribute>(options_.target_scale));
            }
            run_info_->add_attribute("NuHepMC.FluxAveragedTotalCrossSection", D(fatx));

            // Optional per-primary breakdown (opt-in, needs >1 primary). Each
            // primary carries the normalized siren.fatx.<pdg> value plus its raw
            // weight_sum, so partitioned files can be losslessly pooled.
            if(options_.fatx_partition_by_primary
               && options_.fatx_weight_sum_by_primary.size() > 1) {
                run_info_->add_attribute("siren.fatx_partitioned",
                    std::make_shared<HepMC3::IntAttribute>(1));
                for(auto const & kv : options_.fatx_weight_sum_by_primary) {
                    std::string const stub = "siren.fatx." + std::to_string(kv.first);
                    run_info_->add_attribute(stub, D(kGeVm2_to_pb * kv.second));
                    run_info_->add_attribute(stub + ".weight_sum", D(kv.second));
                    auto const cit = options_.accepted_by_primary.find(kv.first);
                    long long const acc = (cit != options_.accepted_by_primary.end())
                                              ? cit->second : 0;
                    run_info_->add_attribute(stub + ".accepted",
                        D(static_cast<double>(acc)));
                }
            }
        }

        // Conventions adhered to (G.R.4). Only "G.C.2" (FATX block above) is
        // declared. E.C.5 is NOT: VectorDoubleAttribute's fixed 6-decimal format
        // truncates sub-microsecond lab time to exactly zero (see lab_pos below).
        if(nuhepmc_conformant) {
            std::vector<std::string> conventions;
            conventions.push_back("G.C.2");
            run_info_->add_attribute("NuHepMC.Conventions",
                std::make_shared<HepMC3::VectorStringAttribute>(conventions));
        }

        if(options_.gzip) {
#ifdef SIREN_HEPMC3_HAS_COMPRESSION
            writer_ = std::make_unique<HepMC3::WriterGZ<HepMC3::WriterAscii>>(filename, run_info_);
#else
            throw std::runtime_error("HepMC3Writer: gzip output requested but this SIREN "
                                     "build's HepMC3 has no compression support");
#endif
        } else {
            writer_ = std::make_unique<HepMC3::WriterAscii>(filename, run_info_);
        }
        if(writer_->failed()) {
            throw std::runtime_error("HepMC3Writer: could not open '" + filename + "' for writing");
        }
    }

    HepMC3::GenEvent TreeToGenEvent(siren::dataclasses::InteractionTree const & tree,
                                    std::uint64_t event_number) const {
        HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::CM);
        evt.set_run_info(run_info_);
        // Same gate as the run-level block: per-event E.R.3/E.R.5 attributes need
        // a backing NuHepMC.Version / process registry.
        bool const nuhepmc_conformant =
            (options_.weights_state != "unweighted") && (options_.accepted_events > 0);
        // HepMC3 stores the event number as int; values above INT_MAX are narrowed
        // here but preserved losslessly as a per-event siren.event_number
        // ULongAttribute (64-bit on SIREN's LP64 targets; would need revisiting
        // on LLP64/Windows).
        evt.set_event_number(static_cast<int>(event_number));
        if(event_number > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
            evt.add_attribute("siren.event_number",
                std::make_shared<HepMC3::ULongAttribute>(
                    static_cast<unsigned long>(event_number)));

        // NuHepMC E.R.3 signal_process_id from the root interaction signature
        // (plain IntAttribute, no NuHepMC. prefix). Ids come from the pre-scan.
        if(nuhepmc_conformant && !tree.tree.empty() && !options_.process_ids.empty()) {
            auto const it = options_.process_ids.find(ProcessKey(tree.tree.front()->record.signature));
            if(it != options_.process_ids.end())
                evt.add_attribute("signal_process_id",
                    std::make_shared<HepMC3::IntAttribute>(it->second));
        }

        // Central-value weight (and any additional named weights) from the header.
        std::size_t const nweights = run_info_->weight_names().size();
        std::vector<double> weights = tree.header.weights;
        if(weights.size() != nweights) {
            double const cv = weights.empty() ? 1.0 : weights.front();
            weights.assign(nweights, cv);
        }
        evt.weights() = weights;

        // A parent vertex's outgoing secondary and the daughter's incoming
        // primary are the same physical particle, joined via ParticleID. Keyed
        // only on set ids -- unset ParticleID compares equal to every other
        // unset id, so unset ids get their own distinct GenParticle instead.
        std::map<siren::dataclasses::ParticleID, HepMC3::GenParticlePtr> particle_by_id;

        // E.R.7 requires exactly one status-4 beam and one status-20 target per
        // event. A tree can hold more than one root (e.g. HepMC3Reader
        // reconstructing a foreign file); primary_root_seen ensures only the
        // first root gets the primary triple and later roots are demoted.
        bool primary_root_seen = false;

        for(auto const & datum_ptr : tree.tree) {
            siren::dataclasses::InteractionTreeDatum const & datum = *datum_ptr;
            siren::dataclasses::InteractionRecord const & rec = datum.record;

            bool const is_primary_root = datum.is_root() && !primary_root_seen;
            if(datum.is_root()) primary_root_seen = true;

            HepMC3::GenVertexPtr vertex =
                std::make_shared<HepMC3::GenVertex>(VertexPosition(rec.interaction_vertex, rec.interaction_time));
            vertex->set_status(VertexStatus(datum, is_primary_root));

            // Incoming primary particle. A non-root datum's primary is the same
            // physical particle as its parent's outgoing secondary. Join only
            // when the id is set, else unset ids would collapse onto one entry.
            HepMC3::GenParticlePtr primary;
            bool primary_is_new = false;
            if(!datum.is_root() && rec.primary_id.IsSet()) {
                auto it = particle_by_id.find(rec.primary_id);
                if(it != particle_by_id.end()) {
                    primary = it->second;
                    primary->set_status(2); // decayed/re-interacted, no longer final
                }
            }
            if(!primary) {
                // Status 4 (incoming beam) is reserved for the one primary root;
                // any other root's fresh incoming particle uses 2.
                int const primary_status = is_primary_root ? 4 : 2;
                primary = std::make_shared<HepMC3::GenParticle>(
                    Momentum(rec.primary_momentum), pdg(rec.signature.primary_type), primary_status);
                primary->set_generated_mass(rec.primary_mass);
                primary_is_new = true;
            }
            vertex->add_particle_in(primary);

            // Target particle (nucleus/nucleon at rest), skipped for decay vertices.
            HepMC3::GenParticlePtr target;
            ParticleType const target_type = rec.signature.target_type;
            if(target_type != ParticleType::Decay && target_type != ParticleType::unknown) {
                // The primary root's target is status 20 (NuHepMC target); every
                // other target is 22, so exactly one status-20 exists (E.R.7).
                int const target_status = is_primary_root ? 20 : 22;
                target = std::make_shared<HepMC3::GenParticle>(
                    HepMC3::FourVector(0, 0, 0, rec.target_mass), pdg(target_type), target_status);
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
                if(j < rec.secondary_ids.size() && rec.secondary_ids[j].IsSet())
                    particle_by_id[rec.secondary_ids[j]] = out;
                vertex->add_particle_out(out);
                outgoing.emplace_back(out, j);
            }

            // Registering the vertex assigns particle/vertex ids; attributes can
            // only be attached once the objects belong to the event.
            evt.add_vertex(vertex);

            if(primary_is_new) {
                // A reused (non-root) primary already got these as the parent's
                // outgoing secondary.
                primary->add_attribute("siren.helicity", D(rec.primary_helicity));
                WriteParticleID(primary, rec.primary_id);
            }
            // primary_initial_position/time are per-record (a daughter's initial
            // position is its parent's interaction vertex, not the root's), so
            // written for every vertex's primary. Three scalar attributes avoid
            // VectorDoubleAttribute, absent from HepMC3 3.2.x.
            {
                std::vector<double> const pos = PositionCM(rec.primary_initial_position);
                primary->add_attribute("siren.primary_initial_position.x", D(pos[0]));
                primary->add_attribute("siren.primary_initial_position.y", D(pos[1]));
                primary->add_attribute("siren.primary_initial_position.z", D(pos[2]));
                primary->add_attribute("siren.primary_initial_time",
                    D(SecondsFromInternal(rec.primary_initial_time)));
            }
            if(target) {
                target->add_attribute("siren.helicity", D(rec.target_helicity));
                WriteParticleID(target, rec.target_id);
            }
            for(auto const & item : outgoing) {
                std::size_t const j = item.second;
                // siren.helicity is required on every particle (strict reader
                // mode); missing entries default to 0.0.
                double const helicity = (j < rec.secondary_helicities.size())
                    ? rec.secondary_helicities[j] : 0.0;
                item.first->add_attribute("siren.helicity", D(helicity));
                if(j < rec.secondary_times.size())
                    item.first->add_attribute("siren.time", D(SecondsFromInternal(rec.secondary_times[j])));
                if(j < rec.secondary_ids.size())
                    WriteParticleID(item.first, rec.secondary_ids[j]);
            }
            // Each interaction_parameters entry {key -> double} becomes a scalar
            // DoubleAttribute "siren.param." + key, raw internal value (no unit
            // conversion). Reader strips the prefix to rebuild the map. Keys are
            // flat ASCII with no embedded dot (e.g. energy, bjorken_x, bjorken_y).
            for(auto const & kv : rec.interaction_parameters) {
                vertex->add_attribute("siren.param." + kv.first, D(kv.second));
            }
            // Decay and unknown-target interactions both omit the target
            // particle above; this flag disambiguates them for the reader.
            if(target_type == ParticleType::unknown) {
                vertex->add_attribute("siren.target_type_unknown",
                    std::make_shared<HepMC3::IntAttribute>(1));
            }
        }

        // Per-event lab position (E.R.5): the primary interaction vertex in CM.
        // Only 3 spatial entries, no 4th time entry (E.C.5): VectorDoubleAttribute's
        // fixed 6-decimal format would truncate sub-microsecond lab time to a false
        // zero. Lossless lab time is the ct slot of this vertex's GenVertex::position().
        if(nuhepmc_conformant && !tree.tree.empty()) {
            siren::dataclasses::InteractionRecord const & root = tree.tree.front()->record;
            std::vector<double> const lab = PositionCM(root.interaction_vertex);
            evt.add_attribute("lab_pos", std::make_shared<HepMC3::VectorDoubleAttribute>(lab));
        }

        return evt;
    }

    void Write(siren::dataclasses::InteractionTree const & tree, std::uint64_t event_number) {
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
    std::unique_ptr<HepMC3::Writer> writer_;
};

HepMC3Writer::HepMC3Writer(std::string const & filename, Options const & options)
    : impl_(std::make_unique<Impl>(filename, options)) {}

HepMC3Writer::~HepMC3Writer() {
    if(impl_) impl_->Close();
}

void HepMC3Writer::Write(siren::dataclasses::InteractionTree const & tree, std::uint64_t event_number) {
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
void HepMC3Writer::Write(siren::dataclasses::InteractionTree const &, std::uint64_t) { ThrowUnsupported(); }
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
    // Pre-scan every tree so run-level metadata (process registry, FATX weight
    // sums, accepted count) is known before GenRunInfo is finalized at construction.
    HepMC3Writer::Options opts = options;
    long long accepted = 0;
    // First pass: accumulate FATX sums/counts and collect one representative
    // signature per distinct root process key (dedup by key, keep the first).
    std::map<std::string, siren::dataclasses::InteractionSignature> sig_by_key;
    for(auto const & tree : trees) {
        if(!tree || tree->tree.empty()) continue;
        siren::dataclasses::InteractionSignature const & sig =
            tree->tree.front()->record.signature;
        std::string const key = ProcessKey(sig);
        sig_by_key.emplace(key, sig);
        // CV weight is slot 0, matching what TreeToGenEvent writes into evt.weights().
        std::vector<double> const & wv = tree->header.weights;
        double const w = wv.empty() ? 1.0 : wv.front();
        opts.fatx_weight_sum += w;
        opts.fatx_weight_sum_by_primary[pdg(sig.primary_type)] += w;
        opts.accepted_by_primary[pdg(sig.primary_type)] += 1;
        ++accepted;
    }
    if(opts.accepted_events < 0) opts.accepted_events = accepted;

    // Second pass: assign process ids in a deterministic signature order so the
    // same physics gets the same id across files, independent of tree order.
    std::vector<siren::dataclasses::InteractionSignature> sigs;
    sigs.reserve(sig_by_key.size());
    for(auto const & kv : sig_by_key) sigs.push_back(kv.second);
    std::sort(sigs.begin(), sigs.end());
    int next_process_id = 700; // NuHepMC generator ("Other") process-id band
    for(auto const & sig : sigs) {
        std::string const key = ProcessKey(sig);
        if(opts.process_ids.find(key) == opts.process_ids.end()) {
            int const id = next_process_id++;
            opts.process_ids[key] = id;
            opts.process_names[id] = ProcessName(sig);
        }
    }

    // Explicit (non-zero) header event numbers, so the running-index fallback
    // below never reassigns one already claimed explicitly.
    std::set<std::uint64_t> explicit_numbers;
    for(auto const & tree : trees) {
        if(tree && !tree->tree.empty() && tree->header.event_number != 0)
            explicit_numbers.insert(tree->header.event_number);
    }

    HepMC3Writer writer(filename, opts);
    std::uint64_t next_fallback = 0;
    for(auto const & tree : trees) {
        // Skip null and empty trees, matching the pre-scan above exactly.
        if(!tree || tree->tree.empty()) continue;
        std::uint64_t event_number;
        if(tree->header.event_number != 0) {
            event_number = tree->header.event_number;
        } else {
            while(explicit_numbers.count(next_fallback)) ++next_fallback;
            event_number = next_fallback++;
        }
        writer.Write(*tree, event_number);
    }
    writer.Close();
}

} // namespace io
} // namespace siren
