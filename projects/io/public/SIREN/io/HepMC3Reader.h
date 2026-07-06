#pragma once
#ifndef SIREN_HepMC3Reader_H
#define SIREN_HepMC3Reader_H

#include <string>
#include <memory>
#include <vector>

namespace siren { namespace dataclasses { struct InteractionTree; } }

namespace siren {
namespace io {

// Read a HepMC3 Ascii file written by SIREN back into InteractionTrees (one tree
// per GenEvent). Carries no HepMC3 types in its signature; throws at runtime when
// SIREN was built without HepMC3.
//
// strict (default true) requires siren.helicity on every particle and
// siren.primary_initial_position.{x,y,z} / siren.primary_initial_time on the
// primary to be present and correctly typed, else throws std::runtime_error.
// If false, missing attributes fall back to 0.0. siren.param.* / siren.time
// are always enumerate-if-present regardless of strict.
std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>>
LoadInteractionTreesFromHepMC3(std::string const & filename, bool strict = true);

} // namespace io
} // namespace siren

#endif // SIREN_HepMC3Reader_H
