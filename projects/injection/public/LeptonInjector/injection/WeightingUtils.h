#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H

#include <memory>

#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"

namespace LI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord const &);

} // namespace injection
} // namespace LI

#endif // LI_WeightingUtils_H
