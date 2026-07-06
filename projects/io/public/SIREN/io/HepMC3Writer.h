#pragma once
#ifndef SIREN_HepMC3Writer_H
#define SIREN_HepMC3Writer_H

#include <map>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>

namespace siren { namespace dataclasses { struct InteractionTree; } }

namespace siren {
namespace io {

// Writes SIREN InteractionTrees to a HepMC3 Ascii file (one GenEvent per tree),
// following NuHepMC conventions. The public interface carries no HepMC3 types;
// when SIREN is built without HepMC3 the methods throw at runtime.
class HepMC3Writer {
public:
    struct Options {
        std::string siren_version;                          // ToolInfo version (may be empty)
        std::vector<std::string> weight_names = {"CV"};     // GenRunInfo weight names
        std::map<std::string, std::string> provenance;      // extra GenRunInfo siren.* attributes

        // Echoed as GenRunInfo StringAttribute "siren.weights_state". Values:
        // "computed" (weighter-computed per-event CV), "header" (per-event weight
        // taken from the tree's header as-is), or "unweighted" (CV is a 1.0
        // placeholder). "unweighted" suppresses all NuHepMC.* attributes and
        // siren.fatx.* diagnostics.
        std::string weights_state = "header";
        // Extra non-PDG particle codes to declare (NuHepMC G.R.11), merged with
        // SIREN's built-in BSM set: code -> {name, description}.
        std::map<int, std::pair<std::string, std::string>> additional_particle_numbers;

        // Run-level generation counts, stored as siren.attempted_events /
        // siren.accepted_events. With < 1 accepted_events the writer emits no
        // NuHepMC.* keys at all. Neither count is a divisor in the FATX formula
        // below. A negative value means "not provided".
        long long attempted_events = -1;   // total sampled including rejected
        long long accepted_events  = -1;   // events saved (auto-filled from tree count if < 0)

        // The Injector's EventsToInject seed (pooled-weighting N_i target). Already
        // baked into each event's CV weight (Weighter::EventWeight), not used in any
        // normalization here. Emitted as siren.events_to_inject when >= 0. A
        // negative value means "not provided".
        long long events_to_inject = -1;

        // Flux-averaged total cross section (NuHepMC G.C.2, G.R.6 units). CV weight
        // (Weighter::EventWeight) already divides by EventsToInject, so
        // fatx_weight_sum summed over accepted events is already the unbiased
        // total cross section in GeV^-2; the writer only applies GeV^-2 -> pb and
        // does NOT also divide by attempted_events/accepted_events. Emitted when
        // weights_state != "unweighted" and >= 1 accepted event.
        // NuHepMC.Units.CrossSection.TargetScale is in-spec only as
        // "PerAtom"/"PerNucleon" (G.R.6): emitted via target_scale only when
        // fatx_per_atom is true, else omitted (value is an unnormalized rate
        // weight). Default false, since SIREN's per-event weight is a rate weight,
        // not a per-atom cross section.
        bool fatx_per_atom = false;
        bool fatx_partition_by_primary = false;  // emit per-primary siren.fatx.<pdg>
        std::string cross_section_unit = "pb";   // NuHepMC.Units.CrossSection.Unit
        std::string target_scale = "PerAtom";    // NuHepMC.Units.CrossSection.TargetScale, when fatx_per_atom

        // Gzip-compress the output (HepMC3 WriterGZ). Requires zlib support in the
        // HepMC3 build; throws at construction if unsupported. Output should carry
        // a .gz suffix (e.g. events.hepmc3.gz) for reader auto-detection.
        bool gzip = false;

        // Internal: process registry + FATX accumulators populated by the pre-scan in
        // SaveInteractionTreesAsHepMC3. Callers normally leave these empty.
        std::map<std::string, int> process_ids;      // signature key -> process id
        std::map<int, std::string> process_names;    // process id -> human name
        double fatx_weight_sum = 0.0;                // sum of CV weights over accepted
        std::map<int, double> fatx_weight_sum_by_primary;  // per primary PDG
        std::map<int, long long> accepted_by_primary;      // per primary PDG accepted count
    };

    HepMC3Writer(std::string const & filename, Options const & options);
    explicit HepMC3Writer(std::string const & filename);
    ~HepMC3Writer();

    HepMC3Writer(HepMC3Writer const &) = delete;
    HepMC3Writer & operator=(HepMC3Writer const &) = delete;

    // Append one tree as a GenEvent with the given event number.
    void Write(siren::dataclasses::InteractionTree const & tree, std::uint64_t event_number);
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Convenience: write every tree to filename (event numbers taken from each
// tree's header, or the running index when the header number is zero).
void SaveInteractionTreesAsHepMC3(
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>> const & trees,
    std::string const & filename,
    HepMC3Writer::Options const & options = HepMC3Writer::Options());

} // namespace io
} // namespace siren

#endif // SIREN_HepMC3Writer_H
