#include "LeptonInjector/injection/TreeWeighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp
#include <initializer_list>                                       // for ini...
#include <iostream>                                               // for ope...
#include <set>                                                    // for set
#include <stdexcept>                                              // for out...
#include "LeptonInjector/interactions/CrossSection.h"            // for Cro...
#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/dataclasses/InteractionSignature.h"      // for Int...
#include "LeptonInjector/detector/DetectorModel.h"                   // for Ear...
#include "LeptonInjector/detector/Coordinates.h"
#include "LeptonInjector/distributions/Distributions.h"           // for Inj...
#include "LeptonInjector/geometry/Geometry.h"                     // for Geo...
#include "LeptonInjector/injection/Injector.h"                // for Inj...
#include "LeptonInjector/injection/Process.h"                     // for Phy...
#include "LeptonInjector/injection/WeightingUtils.h"              // for Cro...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...

#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/injection/Injector.h"

#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/InteractionCollection.h"

#include "LeptonInjector/dataclasses/InteractionSignature.h"

#include <rk/rk.hh>

namespace LI {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

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
} // namespace

//---------------
// class LeptonTreeWeighter
//---------------

void LeptonTreeWeighter::Initialize() {
    int i = 0;
    primary_process_weighters.reserve(injectors.size());
    secondary_process_weighter_maps.reserve(injectors.size());
    for(auto const & injector : injectors) {
        assert(primary_physical_process->MatchesHead(injector->GetPrimaryProcess()));
        primary_process_weighters.push_back(std::make_shared<LeptonProcessWeighter>(LeptonProcessWeighter(primary_physical_process, injector->GetPrimaryProcess(), detector_model)));
        std::map<LI::dataclasses::Particle::ParticleType,
            std::shared_ptr<LeptonProcessWeighter>
                > injector_sec_process_weighter_map;
        std::map<LI::dataclasses::Particle::ParticleType,
            std::shared_ptr<LI::injection::InjectionProcess>
                > injector_sec_process_map = injector->GetSecondaryProcessMap();
        for(auto const & sec_phys_process : secondary_physical_processes) {
            try{
                std::shared_ptr<LI::injection::InjectionProcess> sec_inj_process = injector_sec_process_map.at(sec_phys_process->GetPrimaryType());
                assert(sec_phys_process->MatchesHead(sec_inj_process)); // make sure cross section collection matches
                injector_sec_process_weighter_map[sec_phys_process->GetPrimaryType()] = std::make_shared<LeptonProcessWeighter>(LeptonProcessWeighter(sec_phys_process,sec_inj_process,detector_model));
            } catch(const std::out_of_range& oor) {
                std::cout << "Out of Range error: " << oor.what() << '\n';
                std::cout << "Initialization Incomplete: Particle " <<  sec_phys_process->GetPrimaryType() << " does not exist in injector\n";
                return;
            }
        }
        if(injector_sec_process_weighter_map.size() != injector_sec_process_map.size()) {
            std::cout << "Initialization Incomplete: No one-to-one mapping between injection and physical distributions for injector " << i << "\n";
            return;
        }
        secondary_process_weighter_maps.push_back(injector_sec_process_weighter_map);
    }
}

double LeptonTreeWeighter::EventWeight(LI::dataclasses::InteractionTree const & tree) const {
    // The weight is given by
    //
    // [sum_{injectors i}
    //  x prod_{tree datum d}
    //  x (prod_{generation dist j} p_gen^{idj})
    //  / (prod_{physical dist j} p_phys^{idj}) ] ^-1
    //
    double inv_weight = 0;
    for(unsigned int idx = 0; idx < injectors.size(); ++idx) {
        double physical_probability = 1.0;
        double generation_probability = injectors[idx]->EventsToInject();//GenerationProbability(tree);
        for(auto const & datum : tree.tree) {
            std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds;
            if(datum->depth()==0) {
                bounds = injectors[idx]->InjectionBounds(datum->record);
                physical_probability *= primary_process_weighters[idx]->PhysicalProbability(bounds, datum->record);
                generation_probability *= primary_process_weighters[idx]->GenerationProbability(*datum);
            }
            else {
                try {
                    bounds = injectors[idx]->InjectionBounds(*datum, datum->record.signature.primary_type);
                    double phys_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->PhysicalProbability(bounds, datum->record);
                    double gen_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->GenerationProbability(*datum);
                    physical_probability *= phys_prob;
                    generation_probability *= gen_prob;
                } catch(const std::out_of_range& oor) {
                    std::cout << "Out of Range error: " << oor.what() << '\n';
                    return 0;
                }
            }
        }
        inv_weight += generation_probability / physical_probability;
    }
    return 1./inv_weight;
}

LeptonTreeWeighter::LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<LI::injection::PhysicalProcess>> secondary_physical_processes)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(secondary_physical_processes)
{
    Initialize();
}

LeptonTreeWeighter::LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(std::vector<std::shared_ptr<LI::injection::PhysicalProcess>>())
{
    Initialize();
}

//---------------
// class LeptonProcessWeighter
//---------------

void LeptonProcessWeighter::Initialize() {
    normalization = 1.0;
    for(auto physical_dist : phys_process->GetPhysicalDistributions()) {
        const LI::distributions::PhysicallyNormalizedDistribution* p = dynamic_cast<const LI::distributions::PhysicallyNormalizedDistribution*>(physical_dist.get());
        if(p) {
            if(p->IsNormalizationSet()) {
                normalization *= p->GetNormalization();
            }
        }
    }
    unique_gen_distributions = inj_process->GetInjectionDistributions();
    unique_phys_distributions = phys_process->GetPhysicalDistributions();
    for(std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>>::reverse_iterator gen_it = unique_gen_distributions.rbegin();
            gen_it != unique_gen_distributions.rend(); ++gen_it) {
        for(std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>>::reverse_iterator phys_it = unique_phys_distributions.rbegin();
                phys_it != unique_phys_distributions.rend(); ++phys_it) {
            if((*gen_it) == (*phys_it)) {
                unique_gen_distributions.erase(std::next(gen_it).base());
                unique_phys_distributions.erase(std::next(phys_it).base());
                break;
            }
        }
    }
}

double LeptonProcessWeighter::InteractionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<LI::interactions::CrossSection>>> const & cross_sections_by_target = phys_process->GetInteractions()->GetCrossSectionsByTarget();
    std::vector<LI::dataclasses::Particle::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    double total_decay_length = phys_process->GetInteractions()->TotalDecayLength(record);

    LI::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = detector_model->GetTargetMass(target_xs.first);
        std::vector<std::shared_ptr<LI::interactions::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<LI::dataclasses::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(bounds.first), DetectorPosition(bounds.second), targets, total_cross_sections, total_decay_length);

    double interaction_probability;
    if(total_interaction_depth < 1e-6) {
        interaction_probability = total_interaction_depth;
    } else {
        interaction_probability = one_minus_exp_of_negative(total_interaction_depth);
    }
    return interaction_probability;
}

double LeptonProcessWeighter::NormalizedPositionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<LI::interactions::CrossSection>>> const & cross_sections_by_target = phys_process->GetInteractions()->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LI::dataclasses::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    double total_decay_length = phys_process->GetInteractions()->TotalDecayLength(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = detector_model->GetTargetMass(target_xs.first);
        std::vector<std::shared_ptr<LI::interactions::CrossSection>> const & xs_list = target_xs.second;
        double total_xs = 0.0;
        for(auto const & xs : xs_list) {
            std::vector<LI::dataclasses::InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
            for(auto const & signature : signatures) {
                fake_record.signature = signature;
                // Add total cross section
                total_xs += xs->TotalCrossSection(fake_record);
            }
        }
        total_cross_sections.push_back(total_xs);
    }

    double total_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(bounds.first), DetectorPosition(bounds.second), targets, total_cross_sections, total_decay_length); // unitless
    double traversed_interaction_depth = detector_model->GetInteractionDepthInCGS(intersections, DetectorPosition(bounds.first), DetectorPosition(interaction_vertex), targets, total_cross_sections, total_decay_length);
    double interaction_density = detector_model->GetInteractionDensity(intersections, DetectorPosition(interaction_vertex), targets, total_cross_sections, total_decay_length); //units of m^-1

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

double LeptonProcessWeighter::PhysicalProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds,
        LI::dataclasses::InteractionRecord const & record ) const {

    double physical_probability = 1.0;
    double prob = InteractionProbability(bounds, record);
    physical_probability *= prob;

    prob = NormalizedPositionProbability(bounds, record);
    physical_probability *= prob;

    prob = LI::injection::CrossSectionProbability(detector_model, phys_process->GetInteractions(), record);
    physical_probability *= prob;

    for(auto physical_dist : unique_phys_distributions) {
        physical_probability *= physical_dist->GenerationProbability(detector_model, phys_process->GetInteractions(), record);
    }

    return normalization * physical_probability;
}

double LeptonProcessWeighter::GenerationProbability(LI::dataclasses::InteractionTreeDatum const & datum ) const {
    double gen_probability = LI::injection::CrossSectionProbability(detector_model, phys_process->GetInteractions(), datum.record);

    for(auto gen_dist : unique_gen_distributions) {
        gen_probability *= gen_dist->GenerationProbability(detector_model, phys_process->GetInteractions(), datum);
    }
    return gen_probability;
}

double LeptonProcessWeighter::EventWeight(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds,
        LI::dataclasses::InteractionTreeDatum const & datum) const {
    return PhysicalProbability(bounds,datum.record)/GenerationProbability(datum);
}

LeptonProcessWeighter::LeptonProcessWeighter(std::shared_ptr<LI::injection::PhysicalProcess> phys_process, std::shared_ptr<LI::injection::InjectionProcess> inj_process, std::shared_ptr<LI::detector::DetectorModel> detector_model)
    : phys_process(phys_process)
      , inj_process(inj_process)
      , detector_model(detector_model)
{
    Initialize();
}

} // namespace injection
} // namespace LI
