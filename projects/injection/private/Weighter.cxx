#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/injection/Weighter.h"

#include "LeptonInjector/injection/LeptonInjector.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include <rk/rk.hh>

namespace LeptonInjector {

//---------------
// class LeptonWeighter
//---------------

namespace {
    template <class InIt>
    typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
        typedef typename std::iterator_traits<InIt>::value_type real;
        real sum = real(0);
        real running_error = real(0);
        real temp;
        real difference;

        for (; begin != end; ++begin) {
            difference = *begin;
            difference -= running_error;
            temp = sum;
            temp += difference;
            running_error = temp;
            running_error -= sum;
            running_error -= difference;
            sum = std::move(temp);
        }
        return sum;
    }

    template<typename T>
    T accumulate(std::initializer_list<T> list) {
        return accumulate(list.begin(), list.end());
    }

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
}


double LeptonWeighter::InteractionProbability(std::shared_ptr<InjectorBase const> injector, LI::crosssections::InteractionRecord const & record) const {
    std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds = injector->InjectionBounds(record);
    return InteractionProbability(bounds, record);
}

double LeptonWeighter::InteractionProbability(std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds, LI::crosssections::InteractionRecord const & record) const {
    LI::geometry::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::geometry::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), earth_model->GetEarthCoordDirFromDetCoordDir(primary_direction));
    std::map<LI::utilities::Particle::ParticleType, std::vector<std::shared_ptr<LI::crosssections::CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<LI::utilities::Particle::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    LI::crosssections::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<LI::crosssections::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, bounds.second, targets, total_cross_sections);
    double interaction_probability;
    if(total_interaction_depth < 1e-6) {
        interaction_probability = total_interaction_depth;
    } else {
        interaction_probability = one_minus_exp_of_negative(total_interaction_depth);
    }
    return interaction_probability;
}

double LeptonWeighter::UnnormalizedPositionProbability(std::shared_ptr<InjectorBase const> injector, LI::crosssections::InteractionRecord const & record) const {
    std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds = injector->InjectionBounds(record);
    return UnnormalizedPositionProbability(bounds, record);
}

double LeptonWeighter::UnnormalizedPositionProbability(std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds, LI::crosssections::InteractionRecord const & record) const {
    LI::geometry::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::geometry::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), primary_direction);
    std::map<LI::utilities::Particle::ParticleType, std::vector<std::shared_ptr<LI::crosssections::CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LI::utilities::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    LI::crosssections::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<LI::crosssections::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, bounds.second, targets, total_cross_sections);
    double traversed_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections);
    double interaction_density = earth_model->GetInteractionDensity(intersections, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density;
    } else {
        prob_density = interaction_density * exp(-traversed_interaction_depth);
    }

    return prob_density;
}

double LeptonWeighter::NormalizedPositionProbability(std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds, LI::crosssections::InteractionRecord const & record) const {
    LI::geometry::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::geometry::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), primary_direction);
    std::map<LI::utilities::Particle::ParticleType, std::vector<std::shared_ptr<LI::crosssections::CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LI::utilities::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    LI::crosssections::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<LI::crosssections::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, bounds.second, targets, total_cross_sections);
    double traversed_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections);
    double interaction_density = earth_model->GetInteractionDensity(intersections, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

void LeptonWeighter::Initialize() {
    // Clear distributions
    unique_distributions.clear();
    common_gen_idxs.clear();
    common_phys_idxs.clear();
    distinct_gen_idxs_by_injector.clear();
    distinct_physical_idxs_by_injector.clear();
    unique_contexts.clear();
    context_idx_by_injector.clear();
    normalization = 1.0;

    // Weights are is given by
    //  w = (\sum_i (\prod_j p_gen^ij / p_phys^ij) )^-1
    // We first want to determine which pairs of p_gen^ij and p_phys^ij cancel
    // Secondly we want to determine which p_gen^j are common across all injectors {i}
    //  and similarly which p_phys^j are common across all injectors {i}
    // The calculation can then be simplified by not computing terms that cancel,
    //  pulling out common terms, and finding duplicate terms

    // To do this we will track unique terms in each ratio
    // Initially we assume all term are unique

    // Initialize the state for physical distributions
    // true ==> distribution does not cancel and is not common
    std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> physical_init_state;
    for(auto physical_dist : physical_distributions) {
        physical_init_state.push_back(std::make_pair(true, physical_dist));
        const PhysicallyNormalizedDistribution* p = dynamic_cast<const PhysicallyNormalizedDistribution*>(physical_dist.get());
        if(p) {
            if(p->IsNormalizationSet()) {
                normalization *= p->GetNormalization();
            }
        }
    }
    std::vector<std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>>> physical_distribution_state(injectors.size(), physical_init_state);
    assert(physical_distribution_state.size() == injectors.size());

    // Initialize the state for generation distributions
    // true ==> distribution does not cancel and is not common
    std::vector<std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>>> generation_distribution_state;
    generation_distribution_state.reserve(injectors.size());
    for(auto injector : injectors) {
        std::vector<std::shared_ptr<InjectionDistribution>> dists = injector->GetInjectionDistributions();
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> dist_state;
        dist_state.reserve(dists.size());
        for(auto dist : dists) {
            dist_state.push_back(std::make_pair(true, dist));
        }
        generation_distribution_state.push_back(dist_state);
    }
    assert(generation_distribution_state.size() == injectors.size());

    // Now we can try to identify term that cancel
    for(unsigned int i=0; i<injectors.size(); ++i) {
        // Consider each injector separately
        std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> & phys_dists = physical_distribution_state[i];
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists = generation_distribution_state[i];
        // Must check every pair of physical and injection distribution (unless already cancelled)
        for(unsigned int phys_idx=0; phys_idx<phys_dists.size(); ++phys_idx) {
            std::pair<bool, std::shared_ptr<WeightableDistribution>> & phys_dist = phys_dists[phys_idx];
            if(not phys_dist.first) // Skip if already cancelled
                continue;
            for(unsigned int gen_idx=0; gen_idx<gen_dists.size(); ++gen_idx) {
                std::pair<bool, std::shared_ptr<InjectionDistribution>> & gen_dist = gen_dists[gen_idx];
                if(not gen_dist.first) { // Skip if already cancelled
                    continue;
                }
                // Check if dists are equivalent
                // Must consider the EarthModel and CrossSectionCollection context in the comparison
                std::shared_ptr<WeightableDistribution> gen_dist_ptr(gen_dist.second);
                bool equivalent_dists =
                    phys_dist.second->AreEquivalent( // physical dist
                            earth_model, // physical context
                            cross_sections, // physical context
                            gen_dist_ptr, // generation dist
                            injectors[i]->GetEarthModel(), // generation context
                            injectors[i]->GetCrossSections()); // generation context
                if(not equivalent_dists) {
                    continue;
                }
                phys_dist.first = false;
                gen_dist.first = false;
                break; // This physical dist is cancelled out so we can skip additional comparisons
            }
        }
    }

    // With cancelled terms marked, we can now collect distributions that are common across all terms
    // The one exception to this is vertex position distributions
    // Physical vertex position distributions depend on the injection bounds and so cannot be common across terms

    // Physical distributions have the same EarthModel+CrossSection context so we do not need to compare them
    // We just need to check that these distributions have not been cancelled out for any terms
    std::vector<unsigned int> common_physical_dist_idxs;
    for(unsigned int phys_idx=0; phys_idx<physical_distributions.size(); ++phys_idx) {
        bool has_been_cancelled = false;
        for(unsigned int i=0; i<injectors.size() and not has_been_cancelled; ++i) {
            has_been_cancelled |= (not physical_distribution_state[i][phys_idx].first);
        }
        // Skip distributions that are cancelled out
        if(has_been_cancelled)
            continue;
        // Skip vertex position distributions and note that it is user-supplied
        if(dynamic_cast<const VertexPositionDistribution*>(physical_distributions[phys_idx].get())) {
            user_supplied_position_distribution = true;
            continue;
        }
        // Remove distribution from distinct distributions
        for(unsigned int i=0; i<injectors.size() and not has_been_cancelled; ++i) {
            physical_distribution_state[i][phys_idx].first = false;
        }
        // Add distriution to common distributions
        common_physical_dist_idxs.push_back(phys_idx);
    }

    std::vector<unsigned int> common_generation_dist_idxs;

    unsigned int i=0;
    std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists_0 = generation_distribution_state[i];
    for(unsigned int gen_idx_0=0; gen_idx_0<gen_dists_0.size(); ++gen_idx_0) {
        std::pair<bool, std::shared_ptr<InjectionDistribution>> & gen_dist_0 = gen_dists_0[gen_idx_0];
        if(not gen_dist_0.first)
            continue;
        bool is_common = true;
        std::vector<unsigned int> common_idxs(injectors.size(), 0);
        common_idxs[i] = gen_idx_0;
        for(unsigned int j=i+1; j<injectors.size(); ++j) {
            bool found_common = false;
            std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists_1 = generation_distribution_state[j];
            for(unsigned int gen_idx_1=0; gen_idx_1<gen_dists_1.size(); ++gen_idx_1) {
                std::pair<bool, std::shared_ptr<InjectionDistribution>> & gen_dist_1 = gen_dists_1[gen_idx_1];
                if(not gen_dist_1.first)
                    continue;
                bool equivalent_dists =
                    gen_dist_0.second->AreEquivalent( // gen dist 0
                            injectors[i]->GetEarthModel(), // gen dist 0 context
                            injectors[i]->GetCrossSections(), // gen dist 0 context
                            (std::shared_ptr<WeightableDistribution>)(gen_dist_1.second), // gen dist 1
                            injectors[j]->GetEarthModel(), // gen dist 1 context
                            injectors[j]->GetCrossSections()); // gen dist 1 context
                if(not equivalent_dists)
                    continue;
                found_common = true;
                common_idxs[j] = gen_idx_1;
                break; // We found a gen dist cancelled out so we can skip additional comparisons
            }
            if(not found_common) {
                // No matching distribution in this injector
                // Term is not common across injectors
                is_common = false;
                // We can stop checking other injectors for this term
                break;
            }
        }
        if(not is_common)
            continue;
        // Remove distribution from distinct distribution list
        for(unsigned int inj_idx=0; inj_idx<injectors.size(); ++inj_idx) {
            generation_distribution_state[inj_idx][common_idxs[inj_idx]].first = false;
        }
        // Add distribution to list of common distriubtions
        common_generation_dist_idxs.push_back(common_idxs[0]); // Use the position in the first injector as an ID
    }

    // Now we can collect all the unique distributions
    for(unsigned int gen_idx : common_generation_dist_idxs) {
        // These are common to all injectors, so we pull information from the first injector
        std::shared_ptr<WeightableDistribution> dist = generation_distribution_state[0][gen_idx].second;
        std::shared_ptr<LI::detector::EarthModel> dist_earth = injectors[0]->GetEarthModel();
        std::shared_ptr<LI::crosssections::CrossSectionCollection> dist_cross_sections = injectors[0]->GetCrossSections();
        std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>> p) -> bool {
            return std::get<0>(p)->AreEquivalent(std::get<1>(p), std::get<2>(p), dist, dist_earth, dist_cross_sections);
        };
        auto it = std::find_if(unique_distributions.begin(), unique_distributions.end(), predicate);
        if(it != unique_distributions.end()) {
            unsigned int index = std::distance(unique_distributions.begin(), it);
            common_gen_idxs.push_back(index);
        } else {
            unique_distributions.push_back(std::make_tuple(dist, injectors[0]->GetEarthModel(), injectors[0]->GetCrossSections()));
            common_gen_idxs.push_back(unique_distributions.size()-1);
        }
    }

    for(unsigned int phys_idx : common_physical_dist_idxs) {
        std::shared_ptr<WeightableDistribution> dist = physical_distributions[phys_idx];
        std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>> p) -> bool {
            return std::get<0>(p)->AreEquivalent(std::get<1>(p), std::get<2>(p), dist, earth_model, cross_sections);
        };
        auto it = std::find_if(unique_distributions.begin(), unique_distributions.end(), predicate);
        if(it != unique_distributions.end()) {
            unsigned int index = std::distance(unique_distributions.begin(), it);
            common_phys_idxs.push_back(index);
        } else {
            unique_distributions.push_back(std::make_tuple(dist, earth_model, cross_sections));
            common_phys_idxs.push_back(unique_distributions.size()-1);
        }
    }

    for(unsigned int injector_idx=0; injector_idx<injectors.size(); ++injector_idx) {
        std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> & phys_dists = physical_distribution_state[injector_idx];
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists = generation_distribution_state[injector_idx];

        std::vector<unsigned int> gen_idxs;
        std::vector<unsigned int> phys_idxs;
        for(unsigned int gen_idx=0; gen_idx<gen_dists.size(); ++gen_idx) {
            bool included = gen_dists[gen_idx].first;
            if(not included)
                continue;
            std::shared_ptr<WeightableDistribution> dist = gen_dists[gen_idx].second;
            // These are common to all injectors, so we pull information from the first injector
            std::shared_ptr<LI::detector::EarthModel> dist_earth = injectors[injector_idx]->GetEarthModel();
            std::shared_ptr<LI::crosssections::CrossSectionCollection> dist_cross_sections = injectors[injector_idx]->GetCrossSections();
            std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>> p) -> bool {
                return std::get<0>(p)->AreEquivalent(std::get<1>(p), std::get<2>(p), dist, dist_earth, dist_cross_sections);
            };
            auto it = std::find_if(unique_distributions.begin(), unique_distributions.end(), predicate);
            if(it != unique_distributions.end()) {
                unsigned int index = std::distance(unique_distributions.begin(), it);
                gen_idxs.push_back(index);
            } else {
                unique_distributions.push_back(std::make_tuple(dist, injectors[injector_idx]->GetEarthModel(), injectors[injector_idx]->GetCrossSections()));
                gen_idxs.push_back(unique_distributions.size()-1);
            }
        }

        for(unsigned int phys_idx=0; phys_idx<phys_dists.size(); ++phys_idx) {
            bool included = phys_dists[phys_idx].first;
            if(not included)
                continue;
            std::shared_ptr<WeightableDistribution> dist = phys_dists[phys_idx].second;
            std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>> p) -> bool {
                return std::get<0>(p)->AreEquivalent(std::get<1>(p), std::get<2>(p), dist, earth_model, cross_sections);
            };
            auto it = std::find_if(unique_distributions.begin(), unique_distributions.end(), predicate);
            if(it != unique_distributions.end()) {
                unsigned int index = std::distance(unique_distributions.begin(), it);
                phys_idxs.push_back(index);
            } else {
                unique_distributions.push_back(std::make_tuple(dist, earth_model, cross_sections));
                phys_idxs.push_back(unique_distributions.size()-1);
            }
        }
        distinct_gen_idxs_by_injector.push_back(gen_idxs);
        distinct_physical_idxs_by_injector.push_back(phys_idxs);
    }

    //TODO
    // Find unique contexts
    // std::vector<std::tuple<std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>>> unique_contexts;
    // std::vector<unsigned int> context_idx_by_injector;
}

LeptonWeighter::LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection> cross_sections, std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions)
    : injectors(injectors)
    , earth_model(earth_model)
    , cross_sections(cross_sections)
    , physical_distributions(physical_distributions)
{
    Initialize();
}

double LeptonWeighter::EventWeight(LI::crosssections::InteractionRecord const & record) const {
    // The weight is given by
    //  w = (\sum_i p_gen^i / p_phys^i)^-1

    // The generation probabilities are different between each injector.
    // Most of the physical probabilities are common between all injectors.
    // The physical interaction probability and physical position distribution
    //  depend on the position boundaries of the injection
    //  and thus are different for each injection.
    // Thus the weighting can be given by
    //  w = p_physCommon / (\sum_i p_gen^i / (p_physPos^i * p_physInt^i))

    // However, the normalization of the physical position distribution is identical to the interaction probability.
    // Thus, the two will cancel out and we are left with only the unnormalized position probability
    //  w = p_physCommon / (\sum_i p_gen^i / p_physPosNonNorm^i)

    // The ratio between physical and generation probabilities that differ between injectors
    std::vector<double> gen_over_phys;
    gen_over_phys.reserve(injectors.size());

    // From each injector we need the generation probability and the unnormalized position probability (interaction probability * position probability)
    for(auto injector : injectors) {
        double generation_probability = injector->GenerationProbability(record);
        std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds = injector->InjectionBounds(record);
        double physical_probability = 1.0;

        /*
        if(user_supplied_position_distribution) {
            // Need pos_prob * int_prob
            // pos_prob already supplied
            // just need int_prob
            physical_probability *= InteractionProbability((std::shared_ptr<InjectorBase const>)injector, record);
        } else {
            // Need pos_prob * int_prob
            // nothing is already supplied
            // need pos_prob and int_prob
            // pos_prob * int_prob == unnormalized pos_prob
            physical_probability *= UnnormalizedPositionProbability((std::shared_ptr<InjectorBase const>)injector, record);
        }
        */
        double prob = InteractionProbability(bounds, record);
        physical_probability *= prob;
        prob = NormalizedPositionProbability(bounds, record);
        physical_probability *= prob;
        prob = CrossSectionProbability(injector->GetEarthModel(), injector->GetCrossSections(), record);
        physical_probability *= prob;
        // Number of events is already in GenerationProbability
        // double num_events = injector->EventsToInject();
        gen_over_phys.push_back(generation_probability / physical_probability);
    }

    // The denominator is the sum over the ratios for each injector
    double injection_specific_factors = accumulate(gen_over_phys.begin(), gen_over_phys.end());

    // One physical probability density is computed for each distribution, independent of the injectors
    double common_physical_probability = 1.0;
    for(auto physical_distribution : physical_distributions) {
        double prob = physical_distribution->GenerationProbability(earth_model, cross_sections, record);
        common_physical_probability *= prob;
    }

    double weight = common_physical_probability / injection_specific_factors;
    return normalization * weight;
}

double LeptonWeighter::SimplifiedEventWeight(LI::crosssections::InteractionRecord const & record) const {
    std::vector<double> probs;
    probs.reserve(unique_distributions.size());
    for(unsigned int i=0; i<unique_distributions.size(); ++i) {
        std::tuple<
            std::shared_ptr<WeightableDistribution>,
            std::shared_ptr<LI::detector::EarthModel>,
            std::shared_ptr<LI::crosssections::CrossSectionCollection>
        > const & p = unique_distributions[i];
        probs.push_back(std::get<0>(p)->GenerationProbability(std::get<1>(p), std::get<2>(p), record));
    }

    double phys_over_gen = 1.0;
    for(unsigned int i=0; i<common_phys_idxs.size(); ++i) {
        phys_over_gen *= probs[common_phys_idxs[i]];
    }
    double prob = CrossSectionProbability(earth_model, cross_sections, record);
    phys_over_gen *= prob;
    for(unsigned int i=0; i<common_gen_idxs.size(); ++i) {
        phys_over_gen /= probs[common_gen_idxs[i]];
    }

    std::vector<double> gen_over_phys;
    gen_over_phys.reserve(injectors.size());
    for(unsigned int i=0; i<injectors.size(); ++i) {
        double prob = 1.0;
        prob *= injectors[i]->EventsToInject();
        for(unsigned int j=0; j<distinct_gen_idxs_by_injector[i].size(); ++j) {
            prob *= probs[distinct_gen_idxs_by_injector[i][j]];
        }
        double cross_section_probability = CrossSectionProbability(injectors[i]->GetEarthModel(), injectors[i]->GetCrossSections(), record);
        prob *= cross_section_probability;
        for(unsigned int j=0; j<distinct_physical_idxs_by_injector[i].size(); ++j) {
            prob /= probs[distinct_physical_idxs_by_injector[i][j]];
        }
        /*
        if(user_supplied_position_distribution) {
            // Need pos_prob * int_prob
            // pos_prob already supplied
            // just need int_prob
            double int_prob = InteractionProbability((std::shared_ptr<InjectorBase const>)injectors[i], record);
            prob /= int_prob;
        } else {
            // Need pos_prob * int_prob
            // nothing is already supplied
            // need pos_prob and int_prob
            // pos_prob * int_prob == unnormalized pos_prob
            double pos_prob = UnnormalizedPositionProbability((std::shared_ptr<InjectorBase const>)injectors[i], record);
            prob /= pos_prob;
        }*/
        std::pair<LI::geometry::Vector3D, LI::geometry::Vector3D> bounds = injectors[i]->InjectionBounds(record);
        double interaction_probability = InteractionProbability(bounds, record);
        double normalized_position_probability = NormalizedPositionProbability(bounds, record);
        prob /= interaction_probability;
        prob /= normalized_position_probability;

        // TODO
        // Use unique contexts to compute cross section probability
        gen_over_phys.push_back(prob);
    }

    double gen_over_phys_d = accumulate(gen_over_phys.begin(), gen_over_phys.end());
    double weight = phys_over_gen / gen_over_phys_d;
    return normalization * weight;
}

} // namespace LeptonInjector

