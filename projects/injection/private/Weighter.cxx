#include "SIREN/injection/Weighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp
#include <initializer_list>                                       // for ini...
#include <iostream>                                               // for ope...
#include <set>                                                    // for set
#include <stdexcept>
#include <sstream>
#include <tuple>
#include <cassert>
#include <fstream>
#include <algorithm>
                                             // for out...
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
#include "SIREN/utilities/Errors.h"                       // for Con...


#include "SIREN/injection/Injector.h"

#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"

#include "SIREN/dataclasses/InteractionSignature.h"

#include <rk/rk.hh>

namespace siren {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class Weighter
//---------------

void Weighter::Initialize() {
    // Idempotent: clear any weighters from a prior Initialize() so this can be
    // called again after LoadWeighter() or after injectors are overwritten.
    primary_process_weighters.clear();
    secondary_process_weighter_maps.clear();
    int i = 0;
    primary_process_weighters.reserve(injectors.size());
    secondary_process_weighter_maps.reserve(injectors.size());
    for(auto const & injector : injectors) {
        if(!primary_physical_process->MatchesHead(injector->GetPrimaryProcess())) {
            std::ostringstream oss;
            oss << "Weighter::Initialize: primary physical process (primary type "
                << primary_physical_process->GetPrimaryType()
                << ") does not match the injector primary process (primary type "
                << injector->GetPrimaryProcess()->GetPrimaryType()
                << ") for injector " << i
                << " [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
        primary_process_weighters.push_back(std::make_shared<PrimaryProcessWeighter>(PrimaryProcessWeighter(primary_physical_process, injector->GetPrimaryProcess(), detector_model)));
        std::map<siren::dataclasses::ParticleType, std::shared_ptr<SecondaryProcessWeighter>>
            injector_sec_process_weighter_map;
        std::map<siren::dataclasses::ParticleType, std::shared_ptr<siren::injection::SecondaryInjectionProcess>>
            injector_sec_process_map = injector->GetSecondaryProcessMap();
        for(auto const & sec_phys_process : secondary_physical_processes) {
            try{
                std::shared_ptr<siren::injection::SecondaryInjectionProcess> sec_inj_process = injector_sec_process_map.at(sec_phys_process->GetPrimaryType());
                if(!sec_phys_process->MatchesHead(sec_inj_process)) { // make sure cross section collection matches
                    std::ostringstream oss;
                    oss << "Weighter::Initialize: secondary physical process (primary type "
                        << sec_phys_process->GetPrimaryType()
                        << ") does not match the injector secondary process for injector "
                        << i << " [siren-docs: errors#configuration]";
                    throw siren::utilities::ConfigurationError(oss.str());
                }
                injector_sec_process_weighter_map[sec_phys_process->GetPrimaryType()] =
                    std::make_shared<SecondaryProcessWeighter>(
                            SecondaryProcessWeighter(
                                sec_phys_process,
                                sec_inj_process,detector_model
                            )
                    );
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
        ++i;
    }
}

double Weighter::EventWeight(siren::dataclasses::InteractionTree const & tree) const {
    // The weight is given by
    //
    // [sum_{injectors i}
    //  x prod_{tree datum d}
    //  x (prod_{generation dist j} p_gen^{idj})
    //  / (prod_{physical dist j} p_phys^{idj}) ] ^-1
    //
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



    double inv_weight = 0;
    for(unsigned int idx = 0; idx < injectors.size(); ++idx) {
        double physical_probability = 1.0;
        // Seed with the number of events actually injected so the
        // weight normalizes by the realized sample size.  Fall back to the
        // requested count when weighting before/without generation (InjectedEvents
        // still 0), which keeps the requested-count normalization for that case.
        double generation_probability = injectors[idx]->InjectedEvents();
        if(injectors[idx]->InjectedEvents() == 0) {
            generation_probability = injectors[idx]->EventsToInject();
        }
        for(auto const & datum : tree.tree) {
            std::tuple<siren::math::Vector3D, siren::math::Vector3D> bounds;
            if(datum->is_root()) {
                bounds = injectors[idx]->PrimaryInjectionBounds(datum->record);
                physical_probability *= primary_process_weighters[idx]->PhysicalProbability(bounds, datum->record);
                generation_probability *= primary_process_weighters[idx]->GenerationProbability(*datum);
            }
            else {
                try {
                    bounds = injectors[idx]->SecondaryInjectionBounds(datum->record);
                    double phys_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->PhysicalProbability(bounds, datum->record);
                    double gen_prob = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type)->GenerationProbability(*datum);
                    physical_probability *= phys_prob;
                    generation_probability *= gen_prob;
                } catch(const std::out_of_range& oor) {
                    std::ostringstream oss;
                    oss << "Weighter::EventWeight: no secondary process weighter for secondary type "
                        << datum->record.signature.primary_type
                        << " in injector " << idx
                        << " (" << oor.what() << ") [siren-docs: errors#configuration]";
                    throw siren::utilities::ConfigurationError(oss.str());
                }
            }
        }
        // Weight-calculation guard: a finite
        // generation_probability <= 0 would contribute 0 to inv_weight and, via
        // the reciprocal below, produce an infinite weight; a non-finite physical
        // probability is equally unrecoverable.  Both are configuration/physics
        // failures the caller must not silently absorb.  physical_probability == 0
        // is NOT caught here: it drives inv_weight to +inf and yields a legitimate
        // 0.0-weight event.
        if(generation_probability <= 0.0 || !std::isfinite(physical_probability)) {
            std::ostringstream oss;
            oss << "Weighter::EventWeight: unusable probabilities for injector " << idx;
            if(!tree.tree.empty()) {
                oss << " (primary type "
                    << tree.tree.front()->record.signature.primary_type << ")";
            }
            oss << ": generation_probability=" << generation_probability
                << ", physical_probability=" << physical_probability
                << " [siren-docs: errors#weight-calc]";
            throw siren::utilities::WeightCalculationError(oss.str());
        }
        inv_weight += generation_probability / physical_probability;
    }
    return 1./inv_weight;
}

std::vector<std::shared_ptr<Injector>> const & Weighter::GetInjectors() const {
    return injectors;
}

std::shared_ptr<siren::detector::DetectorModel> Weighter::GetDetectorModel() const {
    return detector_model;
}

std::shared_ptr<siren::injection::PhysicalProcess> Weighter::GetPrimaryPhysicalProcess() const {
    return primary_physical_process;
}

std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> const & Weighter::GetSecondaryPhysicalProcesses() const {
    return secondary_physical_processes;
}

std::vector<double> Weighter::GetInteractionProbabilities(siren::dataclasses::InteractionTree const & tree, int i_inj) const {
    if(i_inj < 0 || static_cast<size_t>(i_inj) >= injectors.size()) {
        throw std::out_of_range("i_inj index out of range in GetInteractionProbabilities");
    }

    std::vector<double> int_probs;
    for(auto const & datum : tree.tree) {
        std::tuple<siren::math::Vector3D, siren::math::Vector3D> bounds;
        if(datum->is_root()) {
            bounds = injectors[i_inj]->PrimaryInjectionBounds(datum->record);
            int_probs.push_back(primary_process_weighters[i_inj]->InteractionProbability(bounds, datum->record));
        }
        else {
            try {
                bounds = injectors[i_inj]->SecondaryInjectionBounds(datum->record);
                int_probs.push_back(secondary_process_weighter_maps[i_inj].at(datum->record.signature.primary_type)->InteractionProbability(bounds, datum->record));
            } catch(const std::out_of_range& oor) {
                std::ostringstream oss;
                oss << "Weighter::GetInteractionProbabilities: no secondary process weighter for secondary type "
                    << datum->record.signature.primary_type
                    << " in injector " << i_inj
                    << " (" << oor.what() << ") [siren-docs: errors#configuration]";
                throw siren::utilities::ConfigurationError(oss.str());
            }
        }
    }
    return int_probs;
}

std::vector<double> Weighter::GetSurvivalProbabilities(siren::dataclasses::InteractionTree const & tree, int i_inj) const {
    if(i_inj < 0 || static_cast<size_t>(i_inj) >= injectors.size()) {
        throw std::out_of_range("i_inj index out of range in GetSurvivalProbabilities");
    }

    std::vector<double> survival_probs;
    for(auto const & datum : tree.tree) {
        std::tuple<siren::math::Vector3D, siren::math::Vector3D> bounds;
        if(datum->is_root()) {
            std::get<0>(bounds) = datum->record.primary_initial_position;
            std::get<1>(bounds) = std::get<0>(injectors[i_inj]->PrimaryInjectionBounds(datum->record));
            survival_probs.push_back(primary_process_weighters[i_inj]->SurvivalProbability(bounds, datum->record));
        }
        else {
            try {
                std::get<0>(bounds) = datum->record.primary_initial_position;
                std::get<1>(bounds) = std::get<0>(injectors[i_inj]->SecondaryInjectionBounds(datum->record));
                survival_probs.push_back(secondary_process_weighter_maps[i_inj].at(datum->record.signature.primary_type)->SurvivalProbability(bounds, datum->record));
            } catch(const std::out_of_range& oor) {
                std::ostringstream oss;
                oss << "Weighter::GetSurvivalProbabilities: no secondary process weighter for secondary type "
                    << datum->record.signature.primary_type
                    << " in injector " << i_inj
                    << " (" << oor.what() << ") [siren-docs: errors#configuration]";
                throw siren::utilities::ConfigurationError(oss.str());
            }
        }
    }
    return survival_probs;
}

void Weighter::SaveWeighter(std::string const & filename) const {
    std::ofstream os(filename+".siren_weighter", std::ios::binary);
    ::cereal::BinaryOutputArchive archive(os);
    this->save(archive,0);
}

void Weighter::LoadWeighter(std::string const & filename) {
    std::ifstream is(filename+".siren_weighter", std::ios::binary);
    ::cereal::BinaryInputArchive archive(is);
    // Read members in the same order Weighter::save writes them; the polymorphic
    // Injector / PhysicalProcess (and the cross sections and decays they hold)
    // are reconstructed via their registered cereal load hooks. Then rebuild the
    // per-process weighters.
    archive(::cereal::make_nvp("Injectors", injectors));
    archive(::cereal::make_nvp("DetectorModel", detector_model));
    archive(::cereal::make_nvp("PrimaryPhysicalProcess", primary_physical_process));
    archive(::cereal::make_nvp("SecondaryPhysicalProcesses", secondary_physical_processes));
    Initialize();
}

Weighter::Weighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> secondary_physical_processes)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(secondary_physical_processes)
{
    Initialize();
}

Weighter::Weighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process)
    : injectors(injectors)
      , detector_model(detector_model)
      , primary_physical_process(primary_physical_process)
      , secondary_physical_processes(std::vector<std::shared_ptr<siren::injection::PhysicalProcess>>())
{
    Initialize();
}

Weighter::Weighter(std::vector<std::shared_ptr<Injector>> _injectors, std::string filename) {
    LoadWeighter(filename);
    if(_injectors.size() > 0) {
        // overwrite the serialized injectors if the user have provided any
        injectors = _injectors;
    }
    Initialize();
}

} // namespace injection
} // namespace siren
