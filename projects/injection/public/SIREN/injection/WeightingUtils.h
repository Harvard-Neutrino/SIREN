#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H

#include <memory>                 // for shared_ptr

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }

namespace SI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<SI::detector::DetectorModel const>, std::shared_ptr<SI::interactions::InteractionCollection const>, SI::dataclasses::InteractionRecord const &);

} // namespace injection
} // namespace SI

#endif // LI_WeightingUtils_H
