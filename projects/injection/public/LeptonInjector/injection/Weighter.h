#pragma once
#ifndef LI_Weighter_H
#define LI_Weighter_H

#include <queue>
#include <memory> // adds shared pointer
#include <iostream>

#include <photospline/bspline.h>
#include <photospline/splinetable.h>
#include <photospline/cinter/splinetable.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/Coordinates.h"
#include "LeptonInjector/Distributions.h"
#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/WeightingUtils.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

namespace LeptonInjector {

class LeptonWeighter {
private:
    std::vector<std::shared_ptr<InjectorBase>> injectors;
    std::shared_ptr<earthmodel::EarthModel> earth_model;
    //TODO Think about whether we want to pass a CrossSection collection, or a vector of cross sections
    //TODO Think about what to do with multiple neutrino primary types. Do we want to support mutiple types across one CrossSectionCollection, across one InjectorBase, across one LeptonWeighter?
    std::shared_ptr<CrossSectionCollection> cross_sections;
    std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions;

    std::vector<std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>> unique_distributions;
    std::vector<unsigned int> common_gen_idxs;
    std::vector<unsigned int> common_phys_idxs;
    std::vector<std::vector<unsigned int>> distinct_gen_idxs_by_injector;
    std::vector<std::vector<unsigned int>> distinct_physical_idxs_by_injector;
    std::vector<std::tuple<std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>> unique_contexts;
    std::vector<unsigned int> context_idx_by_injector;
    double normalization;

    bool user_supplied_position_distribution = false;

    void Initialize();
public:
    //TODO Think about the relationship between interaction probability and the positional distribution. Check that the math works out
    //TODO Add versions of these functions that take precomputed intersections
    double InteractionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const;
    double InteractionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const;
    double UnnormalizedPositionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const;
    double UnnormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const;
    //TODO Add a function to check that we have the right match up of variables between generator and physical distribution
    //TODO Figure out a way to check that physical and generation probabilities match, and ignore those when weighting
    LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<CrossSectionCollection> cross_sections, std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions);
    double EventWeight(InteractionRecord const & record) const;
    double SimplifiedEventWeight(InteractionRecord const & record) const;
};

} //namespace LeptonInjector

CEREAL_CLASS_VERSION(LeptonInjector::LeptonWeighter, 0);

#endif // LI_Weighter_H
