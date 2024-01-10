#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H

#include <memory>                 // for shared_ptr

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }

namespace LI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::InteractionRecord const &);

} // namespace injection
} // namespace LI

#endif // LI_WeightingUtils_H
