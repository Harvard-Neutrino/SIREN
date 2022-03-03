#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/Weighter.h"

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"

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
}


double LeptonWeighter::InteractionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const {
    std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
    return InteractionProbability(bounds, record);
}

double LeptonWeighter::InteractionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<LeptonInjector::Particle::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            total_xs += xs->TotalCrossSection(record);
        }
        total_cross_sections.push_back(total_xs);
    }
    std::vector<double> particle_counts = earth_model->GetParticleColumnDepth(intersections, bounds.first, bounds.second, targets);
    std::vector<double> interaction_probabilities;
    interaction_probabilities.reserve(particle_counts.size());
    for(unsigned int i=0; i<particle_counts.size(); ++i) {
        interaction_probabilities.push_back(particle_counts[i] * total_cross_sections[i]);
    }
    double exponent = accumulate(interaction_probabilities.begin(), interaction_probabilities.end());
    return one_minus_exp_of_negative(exponent);
}

double LeptonWeighter::UnnormalizedPositionProbability(std::shared_ptr<InjectorBase const> injector, InteractionRecord const & record) const {
    std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
    return UnnormalizedPositionProbability(bounds, record);
}

double LeptonWeighter::UnnormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);
    std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LeptonInjector::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            total_xs += xs->TotalCrossSection(record);
        }
        total_cross_sections.push_back(total_xs);
    }

    std::vector<double> particle_counts = earth_model->GetParticleColumnDepth(intersections, bounds.first, interaction_vertex, targets);
    std::vector<double> interaction_probabilities; interaction_probabilities.reserve(n_targets);
    for(unsigned int i=0; i<n_targets; ++i) {
        interaction_probabilities.push_back(particle_counts[i] * total_cross_sections[i]);
    }
    double exponent = accumulate(interaction_probabilities.begin(), interaction_probabilities.end());
    double probability_density = std::exp(-exponent);

    std::vector<double> particle_densities = earth_model->GetParticleDensity(intersections, interaction_vertex, targets.begin(), targets.end());
    std::vector<double> jacobians; jacobians.reserve(n_targets);
    for(unsigned int i=0; i<n_targets; ++ i) {
        jacobians.push_back(particle_densities[i] * total_cross_sections[i]);
    }
    double mass_density = earth_model->GetMassDensity(intersections, interaction_vertex);
    double jacobian = mass_density * accumulate(jacobians.begin(), jacobians.end()) * 100; // cm^-1 --> m^-1

    return probability_density * jacobian;
}

double LeptonWeighter::NormalizedPositionProbability(std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds, InteractionRecord const & record) const {
    double norm = InteractionProbability(bounds, record);
    return UnnormalizedPositionProbability(bounds, record) / norm;
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
    ///std::cerr << "##### Initializing physical distribution states;" << std::endl;
    std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> physical_init_state;
    for(auto physical_dist : physical_distributions) {
        physical_init_state.push_back(std::make_pair(true, physical_dist));
    }
    std::vector<std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>>> physical_distribution_state(injectors.size(), physical_init_state);
    assert(physical_distribution_state.size() == injectors.size());
    ///std::cerr << "##### Initialized physical distribution states;" << std::endl;

    // Initialize the state for generation distributions
    // true ==> distribution does not cancel and is not common
    ///std::cerr << "##### Initializing generation distribution states;" << std::endl;
    std::vector<std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>>> generation_distribution_state;
    generation_distribution_state.reserve(injectors.size());
    unsigned int __injector_idx = 0;
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
    ///std::cerr << "##### Initialized generation distribution states;" << std::endl;


    ///std::cerr << "##### Looking for terms that cancel;" << std::endl;
    // Now we can try to identify term that cancel
    for(unsigned int i=0; i<injectors.size(); ++i) {
        ///std::cerr << "Working with terms from injector " << i << std::endl;
        // Consider each injector separately
        std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> & phys_dists = physical_distribution_state[i];
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists = generation_distribution_state[i];
        // Must check every pair of physical and injection distribution (unless already cancelled)
        for(unsigned int phys_idx=0; phys_idx<phys_dists.size(); ++phys_idx) {
            ///std::cerr << "Working with physical dist "<< phys_idx << std::endl;
            std::pair<bool, std::shared_ptr<WeightableDistribution>> & phys_dist = phys_dists[phys_idx];
            ///if(phys_dist.second)
                ///std::cerr << "Named: " << phys_dist.second->Name() << std::endl;
            ///else
                ///std::cerr << "Named: NULL" << std::endl;
            if(not phys_dist.first) // Skip if already cancelled
                continue;
            for(unsigned int gen_idx=0; gen_idx<gen_dists.size(); ++gen_idx) {
                std::pair<bool, std::shared_ptr<InjectionDistribution>> & gen_dist = gen_dists[gen_idx];
                ///std::cerr << "Working with gen dist "<< phys_idx << std::endl;
                ///if(gen_dist.second) {
                    ///std::cerr << "Named: " << gen_dist.second->Name() << std::endl;
                ///} else {
                    ///std::cerr << "Named: NULL" << std::endl;
                ///}
                if(not gen_dist.first) { // Skip if already cancelled
                    ///std::cerr << gen_dist.second->Name() << " already cancelled" << std::endl;
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
                    ///std::cerr << phys_dist.second->Name() << " != " << gen_dist.second->Name() << std::endl;
                    continue;
                }
                ///std::cerr << phys_dist.second->Name() << " == " << gen_dist.second->Name() << std::endl;
                phys_dist.first = false;
                gen_dist.first = false;
                break; // This physical dist is cancelled out so we can skip additional comparisons
            }
        }
    }
    ///std::cerr << "##### Done looking for terms that cancel;" << std::endl;

    for(unsigned int i=0; i<injectors.size(); ++i) {
        ///std::cerr << "Injector " << i << std::endl;
        // Consider each injector separately
        std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> & phys_dists = physical_distribution_state[i];
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists = generation_distribution_state[i];
        ///std::cerr << "\tGeneration:" << std::endl;
        for(unsigned int j=0; j<gen_dists.size(); ++j) {
            ///std::cerr << "\t\t" << (gen_dists[j].first?"True":"False") << " " << gen_dists[j].second->Name() << std::endl;
        }
        ///std::cerr << "\tPhysical:" << std::endl;
        for(unsigned int j=0; j<phys_dists.size(); ++j) {
            ///std::cerr << "\t\t" << (gen_dists[j].first?"True":"False") << " " << phys_dists[j].second->Name() << std::endl;
        }
    }



    // With cancelled terms marked, we can now collect distributions that are common across all terms
    // The one exception to this is vertex position distributions
    // Physical vertex position distributions depend on the injection bounds and so cannot be common across terms

    // Physical distributions have the same EarthModel+CrossSection context so we do not need to compare them
    // We just need to check that these distributions have not been cancelled out for any terms
    ///std::cerr << "##### Looking for physical common terms that have not cancelled;" << std::endl;
    std::vector<unsigned int> common_physical_dist_idxs;
    for(unsigned int phys_idx=0; phys_idx<physical_distributions.size(); ++phys_idx) {
        ///std::cerr << "Looking at physical distribution " << phys_idx << std::endl;
        bool has_been_cancelled = false;
        ///std::cerr << "Checking if dist has been cancelled already" << std::endl;
        for(unsigned int i=0; i<injectors.size() and not has_been_cancelled; ++i) {
            has_been_cancelled |= (not physical_distribution_state[i][phys_idx].first);
        }
        // Skip distributions that are cancelled out
        if(has_been_cancelled)
            continue;
        // Skip vertex position distributions and note that it is user-supplied
        ///std::cerr << "Checking if distribution is a VertexPositionDistribution" << std::endl;
        if(dynamic_cast<const VertexPositionDistribution*>(physical_distributions[phys_idx].get())) {
            user_supplied_position_distribution = true;
            continue;
        }
        // Remove distribution from distinct distributions
        ///std::cerr << "Removing common distribution" << std::endl;
        for(unsigned int i=0; i<injectors.size() and not has_been_cancelled; ++i) {
            physical_distribution_state[i][phys_idx].first = false;
        }
        // Add distriution to common distributions
        common_physical_dist_idxs.push_back(phys_idx);
    }
    ///std::cerr << "##### Done looking for physical common terms that have not cancelled;" << std::endl;

    std::vector<unsigned int> common_generation_dist_idxs;

    unsigned int i=0;
    ///std::cerr << "##### Looking for generation common terms that have not cancelled;" << std::endl;
    ///std::cerr << "Looking at injector " << i << std::endl;
    std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists_0 = generation_distribution_state[i];
    for(unsigned int gen_idx_0=0; gen_idx_0<gen_dists_0.size(); ++gen_idx_0) {
        ///std::cerr << "Term 0: injector " << i << " gen " << gen_idx_0 << std::endl;
        std::pair<bool, std::shared_ptr<InjectionDistribution>> & gen_dist_0 = gen_dists_0[gen_idx_0];
        if(not gen_dist_0.first)
            continue;
        bool is_common = true;
        std::vector<unsigned int> common_idxs(injectors.size(), 0);
        common_idxs[i] = gen_idx_0;
        for(unsigned int j=i+1; j<injectors.size(); ++j) {
            ///std::cerr << "Looking at injector " << j << std::endl;
            bool found_common = false;
            std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists_1 = generation_distribution_state[j];
            for(unsigned int gen_idx_1=0; gen_idx_1<gen_dists_1.size(); ++gen_idx_1) {
                ///std::cerr << "Term 1: injector " << j << " gen " << gen_idx_1 << std::endl;
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
        ///std::cerr << "Removing common gen distribution" << std::endl;
        for(unsigned int inj_idx=0; inj_idx<injectors.size(); ++inj_idx) {
            generation_distribution_state[inj_idx][common_idxs[inj_idx]].first = false;
        }
        // Add distribution to list of common distriubtions
        common_generation_dist_idxs.push_back(common_idxs[0]); // Use the position in the first injector as an ID
    }
    ///std::cerr << "##### Done looking for generation common terms that have not cancelled;" << std::endl;

    // Now we can collect all the unique distributions
    ///std::cerr << "##### Gathering unique gen distributions;" << std::endl;
    for(unsigned int gen_idx : common_generation_dist_idxs) {
        // These are common to all injectors, so we pull information from the first injector
        std::shared_ptr<WeightableDistribution> dist = generation_distribution_state[0][gen_idx].second;
        std::shared_ptr<earthmodel::EarthModel> dist_earth = injectors[0]->GetEarthModel();
        std::shared_ptr<CrossSectionCollection> dist_cross_sections = injectors[0]->GetCrossSections();
        std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>> p) -> bool {
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
    ///std::cerr << "##### Done gathering unique gen distributions;" << std::endl;

    ///std::cerr << "##### Gathering unique physical distributions;" << std::endl;
    for(unsigned int phys_idx : common_physical_dist_idxs) {
        std::shared_ptr<WeightableDistribution> dist = physical_distributions[phys_idx];
        std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>> p) -> bool {
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
    ///std::cerr << "##### Done gathering unique physical distributions;" << std::endl;

    ///std::cerr << "##### Gathering unique distinct distributions;" << std::endl;
    for(unsigned int injector_idx=0; injector_idx<injectors.size(); ++injector_idx) {
        ///std::cerr << "Looking at injector " << injector_idx << std::endl;
        std::vector<std::pair<bool, std::shared_ptr<WeightableDistribution>>> & phys_dists = physical_distribution_state[injector_idx];
        std::vector<std::pair<bool, std::shared_ptr<InjectionDistribution>>> & gen_dists = generation_distribution_state[injector_idx];

        std::vector<unsigned int> gen_idxs;
        std::vector<unsigned int> phys_idxs;
        ///std::cerr << "Looking at gen dists" << std::endl;
        for(unsigned int gen_idx=0; gen_idx<gen_dists.size(); ++gen_idx) {
            bool included = gen_dists[gen_idx].first;
            if(not included)
                continue;
            std::shared_ptr<WeightableDistribution> dist = gen_dists[gen_idx].second;
            if(dist) {
                ///std::cerr << "Named: " << dist->Name() << std::endl;
            } else {
                ///std::cerr << "Named: NULL" << std::endl;
            }
            // These are common to all injectors, so we pull information from the first injector
            std::shared_ptr<earthmodel::EarthModel> dist_earth = injectors[injector_idx]->GetEarthModel();
            std::shared_ptr<CrossSectionCollection> dist_cross_sections = injectors[injector_idx]->GetCrossSections();
            std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>> p) -> bool {
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
        ///std::cerr << "Done looking at gen dists" << std::endl;

        ///std::cerr << "Looking at physical dists" << std::endl;
        for(unsigned int phys_idx=0; phys_idx<phys_dists.size(); ++phys_idx) {
            bool included = phys_dists[phys_idx].first;
            if(not included)
                continue;
            std::shared_ptr<WeightableDistribution> dist = phys_dists[phys_idx].second;
            if(dist) {
                ///std::cerr << "Named: " << dist->Name() << std::endl;
            } else {
                ///std::cerr << "Named: NULL" << std::endl;
            }
            std::function<bool(std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>)> predicate = [&] (std::tuple<std::shared_ptr<WeightableDistribution>, std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>> p) -> bool {
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
        ///std::cerr << "Done looking at physical dists" << std::endl;
    }
    ///std::cerr << "##### Done gathering unique distinct distributions;" << std::endl;

    ///std::cerr << "Common generation distributions:" << std::endl;
    for(unsigned int i=0; i<common_gen_idxs.size(); ++i) {
        ///std::cerr << "\t" << std::get<0>(unique_distributions[common_gen_idxs[i]])->Name() << std::endl;
    }
    ///std::cerr << "Common physical distributions:" << std::endl;
    for(unsigned int i=0; i<common_phys_idxs.size(); ++i) {
        ///std::cerr << "\t" << std::get<0>(unique_distributions[common_phys_idxs[i]])->Name() << std::endl;
    }
    for(unsigned int i=0; i<injectors.size(); ++i) {
        ///std::cerr << "Injector " << i << ":" << std::endl;
        ///std::cerr << "\t" << "Generation:" << std::endl;
        for(unsigned int j=0; j<distinct_gen_idxs_by_injector[i].size(); ++j) {
            ///std::cerr << "\t\t" << std::get<0>(unique_distributions[distinct_gen_idxs_by_injector[i][j]])->Name() << std::endl;
        }
        ///std::cerr << "\t" << "Physical:" << std::endl;
        for(unsigned int j=0; j<distinct_physical_idxs_by_injector[i].size(); ++j) {
            ///std::cerr << "\t\t" << std::get<0>(unique_distributions[distinct_physical_idxs_by_injector[i][j]])->Name() << std::endl;
        }
    }

    //TODO
    // Find unique contexts
    // std::vector<std::tuple<std::shared_ptr<earthmodel::EarthModel>, std::shared_ptr<CrossSectionCollection>>> unique_contexts;
    // std::vector<unsigned int> context_idx_by_injector;
}

LeptonWeighter::LeptonWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<earthmodel::EarthModel> earth_model, std::shared_ptr<CrossSectionCollection> cross_sections, std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions)
    : injectors(injectors)
    , earth_model(earth_model)
    , cross_sections(cross_sections)
    , physical_distributions(physical_distributions)
{
    Initialize();
}

double LeptonWeighter::EventWeight(InteractionRecord const & record) const {
    ///std::cerr << "Basic event weight" << std::endl;
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
        ///std::cerr << "\tNew Injector" << std::endl;
        double generation_probability = injector->GenerationProbability(record);
        ///std::cerr << "\t\tGenerationProbability: " << generation_probability << std::endl;
        std::pair<earthmodel::Vector3D, earthmodel::Vector3D> bounds = injector->InjectionBounds(record);
        double physical_probability = 1.0;
        if(user_supplied_position_distribution) {
            // Need pos_prob * int_prob
            // pos_prob already supplied
            // just need int_prob
            physical_probability *= InteractionProbability((std::shared_ptr<InjectorBase const>)injector, record);
            ///std::cerr << "\t\tInteractionProbability: " << physical_probability << std::endl;
        } else {
            // Need pos_prob * int_prob
            // nothing is already supplied
            // need pos_prob and int_prob
            // pos_prob * int_prob == unnormalized pos_prob
            physical_probability *= UnnormalizedPositionProbability((std::shared_ptr<InjectorBase const>)injector, record);
            ///std::cerr << "\t\tUnnormalizedPositionProbability: " << physical_probability << std::endl;
        }
        // Number of events is already in GenerationProbability
        // double num_events = injector->EventsToInject();
        // ///std::cerr << "\t\tNumEvents: " << num_events << std::endl;
        gen_over_phys.push_back(generation_probability / physical_probability);
    }

    // The denominator is the sum over the ratios for each injector
    double injection_specific_factors = accumulate(gen_over_phys.begin(), gen_over_phys.end());
    ///std::cerr << "\tInjectionSpecificFactors: " << injection_specific_factors << std::endl;

    // One physical probability density is computed for each distribution, independent of the injectors
    double common_physical_probability = 1.0;
    ///std::cerr << "\tPhysicalProbabilities" << std::endl;
    for(auto physical_distribution : physical_distributions) {
        double prob = physical_distribution->GenerationProbability(earth_model, cross_sections, record);
        ///std::cerr << "\t\t" << physical_distribution->Name() << ": " << prob << std::endl;
        common_physical_probability *= prob;
    }
    ///std::cerr << "\tCommonPhysicalProbability: " << common_physical_probability << std::endl;

    double weight = common_physical_probability / injection_specific_factors;
    ///std::cerr << "\tWeight: " << weight << std::endl;
    return weight;
}

double LeptonWeighter::SimplifiedEventWeight(InteractionRecord const & record) const {
    std::vector<double> probs;
    probs.reserve(unique_distributions.size());
    ///std::cerr << "Simplified Event Weight" << std::endl;
    for(unsigned int i=0; i<unique_distributions.size(); ++i) {
        std::tuple<
            std::shared_ptr<WeightableDistribution>,
            std::shared_ptr<earthmodel::EarthModel>,
            std::shared_ptr<CrossSectionCollection>
        > const & p = unique_distributions[i];
        probs.push_back(std::get<0>(p)->GenerationProbability(std::get<1>(p), std::get<2>(p), record));
    }

    ///std::cerr << "\tCommon physical probs" << std::endl;
    double phys_over_gen = 1.0;
    for(unsigned int i=0; i<common_phys_idxs.size(); ++i) {
        ///std::cerr << "\t\t" << std::get<0>(unique_distributions[common_phys_idxs[i]])->Name() << ": " << probs[common_phys_idxs[i]] << std::endl;
        phys_over_gen *= probs[common_phys_idxs[i]];
    }
    ///std::cerr << "\tCommon gen probs" << std::endl;
    for(unsigned int i=0; i<common_gen_idxs.size(); ++i) {
        ///std::cerr << "\t\t" << std::get<0>(unique_distributions[common_gen_idxs[i]])->Name() << ": " << probs[common_gen_idxs[i]] << std::endl;
        phys_over_gen /= probs[common_gen_idxs[i]];
    }

    std::vector<double> gen_over_phys;
    gen_over_phys.reserve(injectors.size());
    ///std::cerr << "\tInjector specific probs" << std::endl;
    for(unsigned int i=0; i<injectors.size(); ++i) {
        ///std::cerr << "\tNew Injector" << std::endl;
        double prob = 1.0;
        prob *= injectors[i]->EventsToInject();
        ///std::cerr << "\t\tNumEvents: " << prob << std::endl;
        for(unsigned int j=0; j<distinct_gen_idxs_by_injector[i].size(); ++j) {
            ///std::cerr << "\t\t" << std::get<0>(unique_distributions[distinct_gen_idxs_by_injector[i][j]])->Name() << ": " << probs[distinct_gen_idxs_by_injector[i][j]] << std::endl;
            prob *= probs[distinct_gen_idxs_by_injector[i][j]];
        }
        for(unsigned int j=0; j<distinct_physical_idxs_by_injector[i].size(); ++j) {
            ///std::cerr << "\t\t" << std::get<0>(unique_distributions[distinct_physical_idxs_by_injector[i][j]])->Name() << ": " << probs[distinct_physical_idxs_by_injector[i][j]] << std::endl;
            prob /= probs[distinct_physical_idxs_by_injector[i][j]];
        }
        if(user_supplied_position_distribution) {
            // Need pos_prob * int_prob
            // pos_prob already supplied
            // just need int_prob
            double int_prob = InteractionProbability((std::shared_ptr<InjectorBase const>)injectors[i], record);
            ///std::cerr << "\t\tInteractionProbability: " << int_prob << std::endl;
            prob /= int_prob;
        } else {
            // Need pos_prob * int_prob
            // nothing is already supplied
            // need pos_prob and int_prob
            // pos_prob * int_prob == unnormalized pos_prob
            double pos_prob = UnnormalizedPositionProbability((std::shared_ptr<InjectorBase const>)injectors[i], record);
            ///std::cerr << "\t\tUnnormalizedPositionProbability: " << pos_prob << std::endl;
            prob /= pos_prob;
        }

        // TODO
        // Use unique contexts to compute cross section probability
        double cross_section_prob = CrossSectionProbability(injectors[i]->GetEarthModel(), injectors[i]->GetCrossSections(), record);
        ///std::cerr << "\t\tCrossSectionProbability: " << cross_section_prob << std::endl;
        prob *= cross_section_prob;
        gen_over_phys.push_back(prob);
    }

    double gen_over_phys_d = accumulate(gen_over_phys.begin(), gen_over_phys.end());
    ///std::cerr << "\tPhysOverGen: " << phys_over_gen << std::endl;
    ///std::cerr << "\tGenOverPhys: " << gen_over_phys_d << std::endl;
    double weight = phys_over_gen / gen_over_phys_d;
    ///std::cerr << "\tWeight: " << weight << std::endl;
    return weight;
}

double LeptonWeighter::CrossSectionProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections,  InteractionRecord const & record) {
    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();
    std::set<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
    std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);

    double total_prob = 0.0;
    double selected_prob = 0.0;
    for(auto const target : available_targets) {
        if(possible_targets.find(target) != possible_targets.end()) {
            // Get target density
            double target_density = earth_model->GetParticleDensity(intersections, interaction_vertex, target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections->GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<InteractionSignature> signatures = cross_section->GetPossibleSignatures();
                for(auto const & signature : signatures) {
                    // Add total cross section times density to the total prob
                    double target_prob = target_density * cross_section->TotalCrossSection(record);
                    total_prob += target_prob;
                    // Add up total cross section times density times final state prob for matching signatures
                    if(signature == record.signature) {
                        selected_prob += target_prob * cross_section->FinalStateProbability(record);
                    }
                }
            }
        }
    }
    return selected_prob / total_prob;
}

} // namespace LeptonInjector

