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

// Signature-derived process identity (NuHepMC G.R.8). These helpers use no HepMC3
// types, so they live outside the SIREN_HAS_HEPMC3 guard and are available to the
// pre-scan in SaveInteractionTreesAsHepMC3 in both the HepMC3 and stub builds.
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

// 1 GeV^-2 = (hbar c)^2 = 0.3894 mb = 3.894e8 pb (PDG). pdg() is defined in the
// always-compiled anonymous namespace above.
constexpr double kGeVm2_to_pb = 3.894e8;

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

// is_primary_root is true only for the first root datum encountered in a tree:
// NuHepMC V.R.1 permits exactly one vertex with the primary-vertex status code,
// so a forest (an InteractionTree with more than one root -- not produced by
// SIREN's own Injector, but reachable e.g. via HepMC3Reader reconstructing a
// foreign file with disconnected vertex chains) must not mark a second root's
// vertex as primary; it is treated like any other secondary interaction.
int VertexStatus(siren::dataclasses::InteractionTreeDatum const & datum, bool is_primary_root) {
    if(is_primary_root) return 1;                                        // primary vertex
    if(datum.record.signature.target_type == ParticleType::Decay) return 22; // decay
    return 21;                                                           // secondary interaction
}

std::shared_ptr<HepMC3::DoubleAttribute> D(double v) {
    return std::make_shared<HepMC3::DoubleAttribute>(v);
}

// Write a set ParticleID as two SIREN-namespaced attributes. major_id is
// uint64 (ULongAttribute; unsigned long is 64-bit on SIREN's LP64 targets, as
// with siren.event_number) and minor_id is int32 (IntAttribute). An unset id
// writes nothing, so the reader treats an absent siren.id.* as "no stored id"
// and only then falls back to generating one for topology.
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
// info_stub + "[<id>]". The list key is plural (...IDs); the per-id namespace is
// a distinct singular stub (...Info). Exact spelling and attribute types are
// load-bearing -- a strict reader throws on a key or type mismatch.
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

// NuHepMC G.R.11 additional particle numbers. The reference writer and reader
// disagree on the per-code stub (AdditionalParticleNumber vs AdditionalParticleInfo),
// so the Name is emitted under both and the Description under Info; a Description
// is always written (empty allowed) because the reference reader requires one.
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

        // Weight provenance is always echoed. "unweighted" turns off every NuHepMC.*
        // key: the per-event CV weight is only a 1.0 placeholder, so a NuHepMC reader
        // must not be told this is a conforming, rate-normalized file. In that mode the
        // output is a plain HepMC3 file carrying siren.* provenance only.
        run_info_->add_attribute("siren.weights_state",
            std::make_shared<HepMC3::StringAttribute>(options_.weights_state));
        bool const nuhepmc_mode = (options_.weights_state != "unweighted");

        // NuHepMC.Version's mandatory backing (G.R.5 FATX + G.R.6 units + G.R.8
        // process registry) can only be written once there is at least one real,
        // accepted event: with zero accepted events fatx_weight_sum/process_names
        // are both vacuous, so a Version declaration would outrun what the file
        // actually contains. Declare Version only when that backing is available;
        // see the FATX block and the G.R.8 registry below, which reuse this guard.
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

        // Run-level generation counts (metadata; also the FATX normalization). Emitted
        // as DoubleAttribute so large low-acceptance runs (attempts > INT_MAX) do not
        // overflow -- the exact integer is preserved up to 2^53.
        if(options_.attempted_events >= 0)
            run_info_->add_attribute("siren.attempted_events",
                D(static_cast<double>(options_.attempted_events)));
        if(options_.accepted_events >= 0)
            run_info_->add_attribute("siren.accepted_events",
                D(static_cast<double>(options_.accepted_events)));

        // The pooled-weighting seed N_i (Injector EventsToInject). Pure metadata:
        // never used in a normalization here, but persisted so a downstream pooler
        // can reconstruct each file's intended event budget.
        if(options_.events_to_inject >= 0)
            run_info_->add_attribute("siren.events_to_inject",
                D(static_cast<double>(options_.events_to_inject)));

        // NuHepMC process ID registry (G.R.8): one id per distinct root interaction
        // signature, assigned in the generator ("Other", >= 700) band by the pre-scan.
        // Gated on nuhepmc_conformant: an empty registry cannot back a Version claim.
        if(nuhepmc_conformant && !options_.process_names.empty()) {
            std::vector<int> ids;
            std::map<int, std::pair<std::string, std::string>> info;
            for(auto const & kv : options_.process_names) {
                ids.push_back(kv.first);
                info[kv.first] = {kv.second, "SIREN interaction " + kv.second};
            }
            DeclareIdRegistry(*run_info_, "NuHepMC.ProcessIDs", "NuHepMC.ProcessInfo", ids, info);
        }

        // Flux-averaged total cross section, at run level (NuHepMC G.C.2: the value
        // is a single run-wide constant written once, not a per-event running
        // estimate, so G.C.2 -- not E.C.4 -- is the applicable convention).
        //
        // Weight convention traced from siren::injection::Weighter::EventWeight
        // (projects/injection/private/Weighter.cxx): each accepted event's CV
        // weight is physical_probability(x) / (EventsToInject * p_gen(x)), i.e. it
        // already carries a 1/N factor where N = EventsToInject is the injector's
        // configured target event count for the run that produced it. Since events
        // are drawn x_i ~ p_gen, summing the CV weight over the N draws is already
        // the unbiased Monte Carlo estimator of the integral defining the
        // flux-averaged total cross section:
        //   E[sum_i w(x_i)] = N * E_pgen[physical(x) / (N p_gen(x))] = integral physical(x) dx
        // Rejected/failed injection attempts contribute physical_probability == 0
        // and are simply absent from the tree list rather than included as
        // zero-weight terms, which does not bias this sum. Dividing sum(CV) by
        // attempted_events or accepted_events AGAIN would double-normalize: the
        // reported value would shrink like 1/(run size) instead of converging to a
        // constant, so it is not an estimator of any cross section. The fix is to
        // NOT divide a second time -- fatx_weight_sum, summed once over the
        // accepted events, is already the (unnormalized-by-anything-further) GeV^-2
        // estimate; only the pb unit conversion is applied here.
        //
        // The rate normalization is only meaningful for weighted output and only
        // once at least one accepted event exists (nuhepmc_conformant), so under
        // "unweighted" or a zero-event run no siren.fatx.* or NuHepMC FATX keys are
        // emitted.
        if(nuhepmc_conformant) {
            run_info_->add_attribute("siren.fatx.weight_sum", D(options_.fatx_weight_sum));
            double const fatx = kGeVm2_to_pb * options_.fatx_weight_sum;
            // The reserved key + its units back the version declaration.
            run_info_->add_attribute("siren.fatx.value", D(fatx));
            run_info_->add_attribute("NuHepMC.Units.CrossSection.Unit",
                std::make_shared<HepMC3::StringAttribute>(options_.cross_section_unit));
            // NuHepMC.Units.CrossSection.TargetScale is only in-spec as "PerAtom" or
            // "PerNucleon" (G.R.6); SIREN's rate-weight default (fatx_per_atom ==
            // false) is neither, so the key is omitted rather than filled with an
            // out-of-spec placeholder -- a reader must treat an absent TargetScale
            // as "the per-atom/per-nucleon normalization is not asserted".
            if(options_.fatx_per_atom) {
                run_info_->add_attribute("NuHepMC.Units.CrossSection.TargetScale",
                    std::make_shared<HepMC3::StringAttribute>(options_.target_scale));
            }
            run_info_->add_attribute("NuHepMC.FluxAveragedTotalCrossSection", D(fatx));

            // Optional per-primary breakdown (opt-in, needs >1 primary). Each
            // primary carries the normalized siren.fatx.<pdg> value plus the raw
            // ingredient (the CV weight sum) that formed it -- so partitioned files
            // can be losslessly pooled: summing the raw weight_sum across files and
            // rescaling by kGeVm2_to_pb reproduces the combined estimate exactly.
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

        // Conventions adhered to (G.R.4 signalling), built up from whatever this
        // file actually backs so the declared set never outruns the content:
        //  - "G.C.2": the FluxAveragedTotalCrossSection block above, only when
        //    actually emitted (nuhepmc_conformant).
        // E.C.5 (lab time in the lab_pos vector) is NOT declared: HepMC3's
        // VectorDoubleAttribute serializes each entry via a fixed 6-decimal format,
        // which truncates any sub-microsecond lab time to exactly zero seconds (see
        // the lab_pos construction below), so the file cannot actually back that
        // convention's time component.
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
        // Same gate the run-level NuHepMC block uses: per-event NuHepMC attributes
        // (E.R.3 signal_process_id, E.R.5 lab_pos) must not be emitted in a file
        // that carries no NuHepMC.Version / process registry to back them.
        bool const nuhepmc_conformant =
            (options_.weights_state != "unweighted") && (options_.accepted_events > 0);
        // HepMC3's GenEvent stores the event number as int; SIREN event numbers
        // above INT_MAX cannot be represented in the format and are narrowed here.
        // When that narrowing is lossy the full 64-bit identity is preserved as a
        // per-event ULongAttribute the reader prefers over the narrowed field.
        // ULongAttribute (not the wider ULongLongAttribute HepMC3 also offers) is
        // used because the reader looks this attribute up by that exact type;
        // `unsigned long` is 64-bit on SIREN's POSIX/LP64 targets, so no event
        // number is actually truncated here, but this would need revisiting on an
        // LLP64 (Windows) build where unsigned long is only 32 bits.
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

        // Shared particles: a parent vertex's outgoing secondary is the same
        // physical particle as the daughter vertex's incoming primary; the
        // records' ParticleID linkage makes this an exact join (no momentum
        // matching), keyed only on set ids -- an unset ParticleID compares equal
        // to every other unset id, so unset ids are never inserted or looked up
        // here and each such record instead gets its own distinct GenParticle.
        std::map<siren::dataclasses::ParticleID, HepMC3::GenParticlePtr> particle_by_id;

        // E.R.7 requires exactly one status-4 beam and one status-20 target per
        // event. SIREN's own Injector only ever produces single-rooted trees, but
        // an InteractionTree is not structurally prevented from holding more than
        // one root (e.g. HepMC3Reader reconstructing a foreign file with
        // disconnected vertex chains). primary_root_seen tracks whether the one
        // allowed primary-vertex/beam/target triple has already been assigned, so
        // a second root's vertex/primary/target are demoted to the same
        // secondary-interaction status band as any other cascade vertex instead of
        // producing multiple status-4/status-20 particles in one event.
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
            // physical particle as its parent's outgoing secondary (already
            // created and attributed when the parent vertex was processed). The
            // join only fires when the id is actually set: an unset ParticleID
            // compares equal to every other unset id, so joining on it would
            // collapse unrelated particles onto a single map entry.
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
                // any other root's fresh incoming particle uses 2 (not a
                // NuHepMC-standard "beam", but a valid generic incoming-particle
                // status already declared final-state-adjacent in the particle
                // status registry) so E.R.7's beam-uniqueness invariant holds.
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
                // The one primary root's target is the NuHepMC target (status 20);
                // every other target (deeper cascade interactions, and any
                // non-primary root in a forest) is marked 22 so exactly one
                // status-20 target exists per event (E.R.7).
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
                // Helicity and ParticleID for a reused (non-root) primary were
                // already written on it as the parent's outgoing secondary.
                primary->add_attribute("siren.helicity", D(rec.primary_helicity));
                WriteParticleID(primary, rec.primary_id);
            }
            // primary_initial_position/time are per-record quantities -- a
            // daughter's initial position is its parent's interaction vertex, not
            // the root's -- so they are written for every vertex's primary,
            // including the shared particle reused from the parent, and are read
            // back by the secondary vertex-position distributions during
            // reweighting. Three scalar attributes avoid VectorDoubleAttribute,
            // which is absent from HepMC3 3.2.x.
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
                // siren.helicity is unconditional on every particle (the reader's
                // strict mode requires it); a record with fewer helicity entries
                // than secondaries defaults the missing ones to 0.0 rather than
                // producing a file its own strict reader would reject.
                double const helicity = (j < rec.secondary_helicities.size())
                    ? rec.secondary_helicities[j] : 0.0;
                item.first->add_attribute("siren.helicity", D(helicity));
                if(j < rec.secondary_times.size())
                    item.first->add_attribute("siren.time", D(SecondsFromInternal(rec.secondary_times[j])));
                if(j < rec.secondary_ids.size())
                    WriteParticleID(item.first, rec.secondary_ids[j]);
            }
            // siren.param.* namespace contract:
            // each interaction_parameters entry {key -> double} becomes exactly one
            // scalar DoubleAttribute on the vertex named "siren.param." + key, carrying
            // the raw internal value (no unit conversion). The reader rebuilds the map
            // by stripping the "siren.param." prefix. Keys are generator-defined flat
            // ASCII with no embedded dot (e.g. energy, bjorken_x, bjorken_y). This map
            // is reweighting-critical: the Weighter consumes it opaquely via
            // FinalStateProbability, so every key must survive.
            for(auto const & kv : rec.interaction_parameters) {
                vertex->add_attribute("siren.param." + kv.first, D(kv.second));
            }
            // A genuine decay (target_type == Decay) and an unknown-target
            // interaction both omit the target particle above, so the reader
            // cannot tell them apart from the vertex graph alone; this flag
            // disambiguates without touching NuHepMC vertex status.
            if(target_type == ParticleType::unknown) {
                vertex->add_attribute("siren.target_type_unknown",
                    std::make_shared<HepMC3::IntAttribute>(1));
            }
        }

        // Per-event lab position (E.R.5): the primary interaction vertex in CM.
        // lab_pos deliberately carries only the three spatial entries, not a
        // fourth time entry (which would be the E.C.5 convention): HepMC3's
        // VectorDoubleAttribute serializes every entry via a fixed 6-fractional-
        // digit format (std::to_string), so a realistic sub-microsecond lab time
        // would silently truncate to exactly 0.000000 s -- a false zero is worse
        // than an absent value. The authoritative, lossless lab time is the ct
        // slot of this same root vertex's GenVertex::position() (a plain scientific-
        // notation float field, not a std::to_string-serialized attribute); a
        // reader that needs the lab time should use that instead of lab_pos[3].
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
    // Pre-scan every tree once so all run-level metadata (process registry, FATX
    // weight sums, accepted count) is known before the writer's GenRunInfo is
    // finalized at construction -- this avoids relying on late run-info flushing.
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
        // CV weight is slot 0 (the 'CV' weight_names entry), matching what
        // TreeToGenEvent writes into evt.weights(); keep FATX consistent with it.
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

    // Explicit (non-zero) header event numbers, collected up front so the
    // running-index fallback below never reassigns a number another tree in this
    // same file already claims explicitly.
    std::set<std::uint64_t> explicit_numbers;
    for(auto const & tree : trees) {
        if(tree && !tree->tree.empty() && tree->header.event_number != 0)
            explicit_numbers.insert(tree->header.event_number);
    }

    HepMC3Writer writer(filename, opts);
    std::uint64_t next_fallback = 0;
    for(auto const & tree : trees) {
        // Skip null and empty trees, matching the pre-scan above exactly -- a
        // written GenEvent count that disagreed with the accepted/FATX count
        // would desync the file's event count from its own normalization metadata.
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
