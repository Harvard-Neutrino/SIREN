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
std::vector<std::shared_ptr<siren::dataclasses::InteractionTree>>
LoadInteractionTreesFromHepMC3(std::string const & filename);

} // namespace io
} // namespace siren

#endif // SIREN_HepMC3Reader_H
