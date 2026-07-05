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
// following NuHepMC conventions. The public interface carries no HepMC3 types so
// that consumers build without HepMC3 headers; when SIREN is built without
// HepMC3 the methods throw at runtime.
class HepMC3Writer {
public:
    struct Options {
        std::string siren_version;                          // ToolInfo version (may be empty)
        std::vector<std::string> weight_names = {"CV"};     // GenRunInfo weight names
        std::map<std::string, std::string> provenance;      // extra GenRunInfo siren.* attributes

        // Weight provenance, always echoed as the GenRunInfo StringAttribute
        // "siren.weights_state". Allowed values: "computed" (per-event CV weights
        // were computed by a weighter before writing), "header" (per-event weights
        // come from each tree's header, trusted as-is), or "unweighted" (no
        // meaningful per-event weight; the CV slot is a 1.0 placeholder). Only the
        // first two are NuHepMC modes: when weights_state == "unweighted" the file
        // carries NO NuHepMC.* attributes at all (plain HepMC3 with siren.* only),
        // and no siren.fatx.* diagnostics (rate normalization is meaningless).
        std::string weights_state = "header";
        // Extra non-PDG particle codes to declare (NuHepMC G.R.11), merged with
        // SIREN's built-in BSM set: code -> {name, description}.
        std::map<int, std::pair<std::string, std::string>> additional_particle_numbers;

        // Run-level generation counts, stored as metadata (siren.attempted_events /
        // siren.accepted_events). accepted_events also gates whether the file has
        // any real content to back a NuHepMC.Version declaration: with < 1 accepted
        // event the writer stays silent on all NuHepMC.* keys. Neither count is a
        // divisor in the FATX formula below -- see fatx_per_atom. A negative value
        // means "not provided".
        long long attempted_events = -1;   // total sampled including rejected
        long long accepted_events  = -1;   // events saved (auto-filled from tree count if < 0)

        // The Injector's EventsToInject seed (the pooled-weighting N_i target). Not
        // used in any normalization here (it is already baked into each event's CV
        // weight by siren::injection::Weighter::EventWeight); emitted as
        // siren.events_to_inject when >= 0 so a downstream pooler can reconstruct
        // the intended per-file event budget. A negative value means "not provided".
        long long events_to_inject = -1;

        // Flux-averaged total cross section (NuHepMC G.C.2, a run-level constant;
        // G.R.6 units). SIREN's per-event CV weight (Weighter::EventWeight) already
        // divides by the injector's EventsToInject, so fatx_weight_sum -- summed
        // once over the accepted events -- is already the unbiased flux-averaged
        // total cross section in GeV^-2; the writer only applies the GeV^-2 -> pb
        // unit conversion and does NOT divide by attempted_events/accepted_events
        // (that would double-normalize and make the reported value shrink with run
        // size instead of converging). Emitted whenever weights_state != "unweighted"
        // and there is at least one accepted event. NuHepMC.Units.CrossSection.
        // TargetScale is only in-spec as "PerAtom"/"PerNucleon" (G.R.6), so it is
        // only emitted when fatx_per_atom asserts one of those via target_scale;
        // when fatx_per_atom is false the key is omitted rather than filled with an
        // out-of-spec placeholder, and a reader should treat the value as an
        // unnormalized rate weight. The default is false because SIREN's per-event
        // weight is a rate weight, not a per-atom cross section, so claiming a
        // per-atom TargetScale by default would be false.
        bool fatx_per_atom = false;
        bool fatx_partition_by_primary = false;  // emit per-primary siren.fatx.<pdg>
        std::string cross_section_unit = "pb";   // NuHepMC.Units.CrossSection.Unit
        std::string target_scale = "PerAtom";    // NuHepMC.Units.CrossSection.TargetScale, when fatx_per_atom

        // Gzip-compress the output (HepMC3 WriterGZ). Requires a HepMC3 build with
        // zlib support; throws at construction if unsupported. Output should carry a
        // .gz suffix so the reader auto-detects it (e.g. events.hepmc3.gz).
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
