#pragma once
#ifndef LI_WeightingUtils_H
#define LI_WeightingUtils_H


namespace LeptonInjector {

double CrossSectionProbability(std::shared_ptr<earthmodel::EarthModel const>, std::shared_ptr<CrossSectionCollection const>, InteractionRecord const &);

} // namespace LeptonInjector

#endif // LI_WeightingUtils_H
