#pragma once
#ifndef SIREN_Weighter_TCC
#define SIREN_Weighter_TCC
#include "SIREN/injection/Weighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp
#include <initializer_list>                                       // for ini...
#include <iostream>                                               // for ope...
#include <set>                                                    // for set
#include <stdexcept>                                              // for out...

#include "SIREN/interactions/Decay.h"            // for Dec...
#include "SIREN/interactions/CrossSection.h"            // for Cro...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/dataclasses/InteractionSignature.h"      // for Int...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"           // for Inj...
#include "SIREN/geometry/Geometry.h"                     // for Geo...
#include "SIREN/injection/Injector.h"                // for Inj...
#include "SIREN/injection/Process.h"                     // for Phy...
#include "SIREN/injection/WeightingUtils.h"              // for Cro...
#include "SIREN/math/Vector3D.h"                         // for Vec...

#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "SIREN/injection/Injector.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/dataclasses/InteractionSignature.h"

namespace siren {
namespace injection {

double one_minus_exp_of_negative(double x) {
    if(x < 1e-1) {
        return std::exp(std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0);
    } else {
        return 1.0 - std::exp(-x);
    }
}

double log_one_minus_exp_of_negative(double x) {
    if(x < 1e-1) {
        return std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0;
    } else if(x > 3) {
        double ex = std::exp(-x);
        double ex2 = ex * ex;
        double ex3 = ex2 * ex;
        double ex4 = ex3 * ex;
        double ex5 = ex4 * ex;
        double ex6 = ex5 * ex;
        return -(ex + ex2 / 2.0 + ex3 / 3.0 + ex4 / 4.0 + ex5 / 5.0 + ex6 / 6.0);
    } else {
        return std::log(1.0 - std::exp(-x));
    }
}

//---------------
// class ProcessWeighter
//---------------

template<typename ProcessType>
void ProcessWeighter<ProcessType>::Initialize() {
    normalization = 1.0;
    for(auto physical_dist : phys_process->GetPhysicalDistributions()) {
        const siren::distributions::PhysicallyNormalizedDistribution* p = dynamic_cast<const siren::distributions::PhysicallyNormalizedDistribution*>(physical_dist.get());
        if(p) {
            if(p->IsNormalizationSet()) {
                normalization *= p->GetNormalization();
            }
        }
    }
    unique_gen_distributions = GetInjectionDistributions();
    unique_phys_distributions = phys_process->GetPhysicalDistributions();
    for(typename std::vector<std::shared_ptr<typename ProcessType::InjectionType>>::reverse_iterator gen_it = unique_gen_distributions.rbegin();
            gen_it != unique_gen_distributions.rend(); ++gen_it) {
        for(std::vector<std::shared_ptr<siren::distributions::WeightableDistribution>>::reverse_iterator phys_it = unique_phys_distributions.rbegin();
                phys_it != unique_phys_distributions.rend(); ++phys_it) {
            if((*gen_it) == (*phys_it)) {
                unique_gen_distributions.erase(std::next(gen_it).base());
                unique_phys_distributions.erase(std::next(phys_it).base());
                break;
            }
        }
    }
}

template<typename ProcessType>
double ProcessWeighter<ProcessType>::InteractionProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const {
    using siren::detector::DetectorPosition;
    using siren::detector::DetectorDirection;
    siren::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    siren::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    siren::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<siren::interactions::CrossSection>>> const & cross_sections_by_target = phys_process->GetInteractions()->GetCrossSectionsByTarget();
    std::vector<siren::dataclasses::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    double total_decay_length = phys_process->GetInteractions()->TotalDecayLength(record);

    siren::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = detector_model->GetTargetMass(target_xs.first);
        std::vector<std::shared_ptr<siren::interactions::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<siren::dataclasses::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(std::get<0>(bounds)), DetectorPosition(std::get<1>(bounds)), targets, total_cross_sections, total_decay_length);

    double interaction_probability;
    if(total_interaction_depth < 1e-6) {
        interaction_probability = total_interaction_depth;
    } else {
        interaction_probability = one_minus_exp_of_negative(total_interaction_depth);
    }
    return interaction_probability;
}

template<typename ProcessType>
double ProcessWeighter<ProcessType>::NormalizedPositionProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const {
    using siren::detector::DetectorPosition;
    using siren::detector::DetectorDirection;
    siren::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    siren::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    siren::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<siren::interactions::CrossSection>>> const & cross_sections_by_target = phys_process->GetInteractions()->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<siren::dataclasses::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    double total_decay_length = phys_process->GetInteractions()->TotalDecayLength(record);
    siren::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = detector_model->GetTargetMass(target_xs.first);
        std::vector<std::shared_ptr<siren::interactions::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<siren::dataclasses::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(std::get<0>(bounds)), DetectorPosition(std::get<1>(bounds)), targets, total_cross_sections, total_decay_length); // unitless
    double traversed_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(std::get<0>(bounds)), DetectorPosition(interaction_vertex), targets, total_cross_sections, total_decay_length);
    double interaction_density = detector_model->GetInteractionDensity(intersections, DetectorPosition(interaction_vertex), targets, total_cross_sections, total_decay_length); //units of m^-1

    double prob_density;
    // This is equivalent to equation 11 of the SIREN paper
    // Reach out to the authors if you disagree and we can send the derivation :)
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

template<typename ProcessType>
double ProcessWeighter<ProcessType>::PhysicalProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds,
        siren::dataclasses::InteractionRecord const & record ) const {

    double physical_probability = 1.0;
    double prob = InteractionProbability(bounds, record);
    physical_probability *= prob;

    prob = NormalizedPositionProbability(bounds, record);
    physical_probability *= prob;

    prob = siren::injection::CrossSectionProbability(detector_model, phys_process->GetInteractions(), record);
    physical_probability *= prob;

    for(auto physical_dist : unique_phys_distributions) {
        physical_probability *= physical_dist->GenerationProbability(detector_model, phys_process->GetInteractions(), record);
    }

    return normalization * physical_probability;
}

template<typename ProcessType>
double ProcessWeighter<ProcessType>::GenerationProbability(siren::dataclasses::InteractionTreeDatum const & datum ) const {
    double gen_probability = siren::injection::CrossSectionProbability(detector_model, inj_process->GetInteractions(), datum.record);

    for(auto gen_dist : unique_gen_distributions) {
        gen_probability *= gen_dist->GenerationProbability(detector_model, inj_process->GetInteractions(), datum.record);
    }
    return gen_probability;
}

template<typename ProcessType>
double ProcessWeighter<ProcessType>::EventWeight(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds,
        siren::dataclasses::InteractionTreeDatum const & datum) const {
    return PhysicalProbability(bounds,datum.record)/GenerationProbability(datum);
}

template<typename ProcessType>
ProcessWeighter<ProcessType>::ProcessWeighter(std::shared_ptr<siren::injection::PhysicalProcess> phys_process, std::shared_ptr<ProcessType> inj_process, std::shared_ptr<siren::detector::DetectorModel> detector_model)
    : phys_process(phys_process)
      , inj_process(inj_process)
      , detector_model(detector_model)
{
    Initialize();
}

//template<typename ProcessType>
//std::vector<std::shared_ptr<ProcessType::InjectionType>> const & ProcessWeighter<ProcessType>::GetInjectionDistributions() {
//    return unique_gen_distributions;
//}

template<>
std::vector<std::shared_ptr<siren::distributions::PrimaryInjectionDistribution>> const & PrimaryProcessWeighter::GetInjectionDistributions() {
    return inj_process->GetPrimaryInjectionDistributions();
}

template<>
std::vector<std::shared_ptr<siren::distributions::SecondaryInjectionDistribution>> const & SecondaryProcessWeighter::GetInjectionDistributions() {
    return inj_process->GetSecondaryInjectionDistributions();
}

} // namespace injection
} // namespace siren

#endif // SIREN_Weighter_TCC
