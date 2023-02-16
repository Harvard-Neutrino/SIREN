#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H


namespace LeptonInjector {

double CrossSectionProbability(std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::crosssections::InteractionRecord const &);

} // namespace LeptonInjector

#endif // LI_WeightingUtils_H
