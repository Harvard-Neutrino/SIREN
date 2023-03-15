#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/injection/TreeWeighter.h"

#include "LeptonInjector/injection/InjectorBase.h"

#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"

#include "LeptonInjector/dataclasses/InteractionSignature.h"

#include <rk/rk.hh>

namespace LI {
namespace injection {

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

//---------------
// class LeptonTreeWeighter
//---------------

void LeptonTreeWeighter::Initialize() {
  int i = 0;
  primary_process_weighters.reserve(injectors.size());
  secondary_process_weighter_maps.reserve(injectors.size());
  for(auto const & injector : injectors) {
    assert(primary_physical_process->MatchesHead(injector->GetPrimaryProcess()));
    primary_process_weighters.push_back(std::make_shared<LeptonProcessWeighter>(LeptonProcessWeighter(primary_physical_process,injector->GetPrimaryProcess(),earth_model)));
    std::map<LI::dataclasses::Particle::ParticleType,
             std::shared_ptr<LeptonProcessWeighter>
    > injector_sec_process_weighter_map;
    std::map<LI::dataclasses::Particle::ParticleType,
             std::shared_ptr<LI::dataclasses::InjectionProcess>
    > injector_sec_process_map = injector->GetSecondaryProcessMap();
    for(auto const & sec_phys_process : secondary_physical_processes) {
      try{
        std::shared_ptr<LI::dataclasses::InjectionProcess> sec_inj_process = injector_sec_process_map.at(sec_phys_process->primary_type);
        assert(sec_phys_process->MatchesHead(sec_inj_process)); // make sure cross section collection matches
        injector_sec_process_weighter_map[sec_phys_process->primary_type] = std::make_shared<LeptonProcessWeighter>(LeptonProcessWeighter(sec_phys_process,sec_inj_process,earth_model));
      } catch(const std::out_of_range& oor) {
        std::cout << "Out of Range error: " << oor.what() << '\n';
        std::cout << "Initialization Incomplete: Particle " <<  sec_phys_process->primary_type << " does not exist in injector\n";
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
    double generation_probability = injectors[idx]->GenerationProbability(tree);
    for(auto const datum : tree.tree) {
      std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds;
      if(datum->depth()==0) {
        physical_probability *= primary_process_weighters[idx]->PhysicalProbability(bounds, datum->record);
        bounds = injectors[idx]->InjectionBounds(datum->record);
      }
      else {
        try {
          double phys_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->PhysicalProbability(bounds, datum->record);
          physical_probability *= phys_prob;
          bounds = injectors[idx]->InjectionBounds(datum->record, datum->record.signature.primary_type);
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

LeptonTreeWeighter::LeptonTreeWeighter(std::vector<std::shared_ptr<InjectorBase>> injectors, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<LI::dataclasses::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<LI::dataclasses::PhysicalProcess>> secondary_physical_processes)
    : injectors(injectors)
    , earth_model(earth_model)
    , primary_physical_process(primary_physical_process)
    , secondary_physical_processes(secondary_physical_processes)
{
  Initialize();
}

//---------------
// class LeptonProcessWeighter
//---------------

// TODO: find distributions that cancel
void LeptonProcessWeighter::Initialize() {
  normalization = 1.0;
  for(auto physical_dist : phys_process->physical_distributions) {
    const LI::distributions::PhysicallyNormalizedDistribution* p = dynamic_cast<const LI::distributions::PhysicallyNormalizedDistribution*>(physical_dist.get());
    if(p) {
      if(p->IsNormalizationSet()) {
        normalization *= p->GetNormalization();
      }
    }
  }
}

double LeptonProcessWeighter::InteractionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> & bounds, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), earth_model->GetEarthCoordDirFromDetCoordDir(primary_direction));
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<LI::crosssections::CrossSection>>> const & cross_sections_by_target = phys_process->cross_sections->GetCrossSectionsByTarget();
    std::vector<LI::dataclasses::Particle::ParticleType> targets;
    targets.reserve(cross_sections_by_target.size());
    std::vector<double> total_cross_sections;
    double total_decay_width = phys_process->cross_sections->TotalDecayWidth(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & xs_list = target_xs.second;
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

    double total_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, bounds.second, targets, total_cross_sections, total_decay_width);
    double interaction_probability;
    if(total_interaction_depth < 1e-6) {
        interaction_probability = total_interaction_depth;
    } else {
        interaction_probability = one_minus_exp_of_negative(total_interaction_depth);
    }
    return interaction_probability;
}

double LeptonProcessWeighter::NormalizedPositionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds, LI::dataclasses::InteractionRecord const & record) const {
    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), primary_direction);
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<LI::crosssections::CrossSection>>> const & cross_sections_by_target = phys_process->cross_sections->GetCrossSectionsByTarget();

    unsigned int n_targets = cross_sections_by_target.size();

    std::vector<LI::dataclasses::Particle::ParticleType> targets; targets.reserve(n_targets);
    std::vector<double> total_cross_sections;
    double total_decay_width = phys_process->cross_sections->TotalDecayWidth(record);
    LI::dataclasses::InteractionRecord fake_record = record;
    for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
        fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & xs_list = target_xs.second;
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

    double total_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, bounds.second, targets, total_cross_sections, total_decay_width);
    double traversed_interaction_depth = earth_model->GetInteractionDepthInCGS(intersections, bounds.first, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections, total_decay_width);
    double interaction_density = earth_model->GetInteractionDensity(intersections, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), targets, total_cross_sections, total_decay_width);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }

    return prob_density;
}

double LeptonProcessWeighter::PhysicalProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> & bounds,
                                                  LI::dataclasses::InteractionRecord const & record ) const {
        double physical_probability = 1.0;
        double prob = InteractionProbability(bounds, record);
        physical_probability *= prob;
        prob = NormalizedPositionProbability(bounds, record);
        physical_probability *= prob;
        prob = LI::injection::CrossSectionProbability(earth_model, phys_process->cross_sections, record);
        physical_probability *= prob;
        for(auto physical_dist : phys_process->physical_distributions) {
          physical_probability *= physical_dist->GenerationProbability(earth_model, phys_process->cross_sections, record);
        }
        return normalization * physical_probability; 
}

// TODO: implement smart EventWeight function that cancels common distributions
double LeptonProcessWeighter::EventWeight(std::pair<LI::math::Vector3D, LI::math::Vector3D> bounds,
                                          LI::dataclasses::InteractionRecord const & record) const {
}

LeptonProcessWeighter::LeptonProcessWeighter(std::shared_ptr<LI::dataclasses::PhysicalProcess> phys_process,std::shared_ptr<LI::dataclasses::InjectionProcess> inj_process, std::shared_ptr<LI::detector::EarthModel> earth_model)
    : phys_process(phys_process)
    , inj_process(inj_process)
    , earth_model(earth_model)
{
  Initialize();
}

} // namespace injection
} // namespace LI
