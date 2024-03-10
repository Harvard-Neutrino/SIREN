#pragma once
#ifndef SIREN_WeightingUtils_H
#define SIREN_WeightingUtils_H

#include <memory>                 // for shared_ptr

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace injection {

double CrossSectionProbability(std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &);

} // namespace injection
} // namespace siren

#endif // SIREN_WeightingUtils_H
