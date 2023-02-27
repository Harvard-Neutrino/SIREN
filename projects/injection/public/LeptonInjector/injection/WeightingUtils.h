#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H


namespace LI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::crosssections::InteractionRecord const &);

} // namespace injection
} // namespace LI

#endif // LI_WeightingUtils_H
