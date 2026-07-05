#pragma once
#ifndef SIREN_HepMC3Writer_H
#define SIREN_HepMC3Writer_H

#include <map>
#include <string>
#include <memory>
#include <vector>
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
    };

    HepMC3Writer(std::string const & filename, Options const & options);
    explicit HepMC3Writer(std::string const & filename);
    ~HepMC3Writer();

    HepMC3Writer(HepMC3Writer const &) = delete;
    HepMC3Writer & operator=(HepMC3Writer const &) = delete;

    // Append one tree as a GenEvent with the given event number.
    void Write(siren::dataclasses::InteractionTree const & tree, int event_number);
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
