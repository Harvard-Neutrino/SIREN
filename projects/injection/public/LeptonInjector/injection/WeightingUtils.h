#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H

#include <memory>                 // for shared_ptr

namespace LI { namespace crosssections { class CrossSectionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class EarthModel; } }

namespace LI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord const &);

} // namespace injection
} // namespace LI

#endif // LI_WeightingUtils_H
