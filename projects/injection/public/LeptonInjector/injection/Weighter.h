#pragma once
#ifndef LI_Weighter_H
#define LI_Weighter_H

#include <tuple>                  // for tuple
#include <memory>                 // for shared_ptr
#include <vector>                 // for vector
#include <utility>                // for pair

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace injection { class Injector; } }
namespace LI { namespace math { class Vector3D; } }

namespace LI {
namespace injection {

class LeptonWeighter {
private:
    std::vector<std::shared_ptr<Injector>> injectors;
    std::shared_ptr<LI::detector::DetectorModel> detector_model;
    //TODO Think about whether we want to pass a CrossSection collection, or a vector of cross sections
    //TODO Think about what to do with multiple neutrino primary types. Do we want to support mutiple types across one InteractionCollection, across one Injector, across one LeptonWeighter?
    std::shared_ptr<LI::interactions::InteractionCollection> interactions;
    std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>> physical_distributions;

    std::vector<std::tuple<std::shared_ptr<LI::distributions::WeightableDistribution>, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<LI::interactions::InteractionCollection>>> unique_distributions;
    std::vector<unsigned int> common_gen_idxs;
    std::vector<unsigned int> common_phys_idxs;
    std::vector<std::vector<unsigned int>> distinct_gen_idxs_by_injector;
    std::vector<std::vector<unsigned int>> distinct_physical_idxs_by_injector;
    std::vector<std::tuple<std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<LI::interactions::InteractionCollection>>> unique_contexts;
    std::vector<unsigned int> context_idx_by_injector;
    double normalization;

    bool user_supplied_position_distribution = false;

    void Initialize();
public:
    //TODO Think about the relationship between interaction probability and the positional distribution. Check that the math works out
    //TODO Add versions of these functions that take precomputed intersections
    double InteractionProbability(std::shared_ptr<Injector const> injector, LI::dataclasses::InteractionRecord const & record) const;
    double InteractionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds, LI::dataclasses::InteractionRecord const & record) const;
    double UnnormalizedPositionProbability(std::shared_ptr<Injector const> injector, LI::dataclasses::InteractionRecord const & record) const;
    double UnnormalizedPositionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds, LI::dataclasses::InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds, LI::dataclasses::InteractionRecord const & record) const;
    //TODO Add a function to check that we have the right match up of variables between generator and physical distribution
    //TODO Figure out a way to check that physical and generation probabilities match, and ignore those when weighting
    LeptonWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<LI::interactions::InteractionCollection> interactions, std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>> physical_distributions);
    double EventWeight(LI::dataclasses::InteractionRecord const & record) const;
    double SimplifiedEventWeight(LI::dataclasses::InteractionRecord const & record) const;
};

} //namespace injection
} //namespace LI

CEREAL_CLASS_VERSION(LI::injection::LeptonWeighter, 0);


#endif // LI_Weighter_H

