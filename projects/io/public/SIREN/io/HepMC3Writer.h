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
        // Extra non-PDG particle codes to declare (NuHepMC G.R.11), merged with
        // SIREN's built-in BSM set: code -> {name, description}.
        std::map<int, std::pair<std::string, std::string>> additional_particle_numbers;

        // Run-level generation counts, stored as metadata to normalize
        // the flux-averaged cross section. A negative value means "not provided".
        long long attempted_events = -1;   // total sampled including rejected
        long long accepted_events  = -1;   // events saved (auto-filled from tree count if < 0)

        // Flux-averaged total cross section (NuHepMC E.C.4 / G.R.6). SIREN's per-event
        // weight is a *rate* weight, so sum(weights)/attempted is only a true per-atom
        // cross section when the physical flux is unit-normalized and the target
        // column-density normalization is divided out. The siren.fatx.* diagnostics are
        // always written; the reserved NuHepMC.FluxAveragedTotalCrossSection key is only
        // emitted when fatx_per_atom is set (opt-in that the value is a per-atom sigma).
        bool fatx_per_atom = false;
        bool fatx_partition_by_primary = false;  // emit per-primary siren.fatx.<pdg>
        std::string cross_section_unit = "pb";   // NuHepMC.Units.CrossSection.Unit
        std::string target_scale = "PerAtom";    // NuHepMC.Units.CrossSection.TargetScale
        int cv_weight_index = 0;                 // weight slot used as the CV rate weight

        // Internal: process registry + FATX accumulators populated by the pre-scan in
        // SaveInteractionTreesAsHepMC3. Callers normally leave these empty.
        std::map<std::string, int> process_ids;      // signature key -> process id
        std::map<int, std::string> process_names;    // process id -> human name
        double fatx_weight_sum = 0.0;                // sum of CV weights over accepted
        std::map<int, double> fatx_weight_sum_by_primary;  // per primary PDG
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
