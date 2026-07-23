#include "SIREN/injection/Weighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp, isfinite
#include <initializer_list>                                       // for ini...
#include <iostream>                                               // for ope...
#include <limits>                                                 // for quiet_NaN
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
                std::ostringstream oss;
                oss << "Weighter::Initialize: secondary physical process (primary type "
                    << sec_phys_process->GetPrimaryType()
                    << ") has no matching injector secondary process for injector "
                    << i << " (" << oor.what() << ") [siren-docs: errors#configuration]";
                throw siren::utilities::ConfigurationError(oss.str());
            }
        }
        if(injector_sec_process_weighter_map.size() != injector_sec_process_map.size()) {
            std::ostringstream oss;
            oss << "Weighter::Initialize: no one-to-one mapping between injection ("
                << injector_sec_process_map.size() << ") and physical ("
                << injector_sec_process_weighter_map.size()
                << ") secondary processes for injector " << i
                << " [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
        secondary_process_weighter_maps.push_back(injector_sec_process_weighter_map);
        ++i;
    }
}

// Computes one tree datum's physical and generation factors for injector idx,
// via the same ProcessWeighter calls EventWeight consumes. depth is left at
// its default for non-root data; the caller fills it in from the owning tree.
// with_diagnostics adds the observation-only interaction_prob/position_prob.
VertexWeightFactors Weighter::ComputeVertexFactors(unsigned int idx,
        std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum,
        bool with_diagnostics) const {
    VertexWeightFactors factors;
    factors.injector_index = static_cast<int>(idx);
    factors.vertex_pdg = static_cast<int>(datum->record.signature.primary_type);
    std::tuple<siren::math::Vector3D, siren::math::Vector3D> bounds;
    if(datum->is_root()) {
        factors.depth = 0;
        bounds = injectors[idx]->PrimaryInjectionBounds(datum->record);
        factors.physical = primary_process_weighters[idx]->PhysicalProbability(bounds, datum->record);
        factors.generation = primary_process_weighters[idx]->GenerationProbability(*datum);
        if(with_diagnostics) {
            factors.interaction_prob = primary_process_weighters[idx]->InteractionProbability(bounds, datum->record);
            factors.position_prob = primary_process_weighters[idx]->NormalizedPositionProbability(bounds, datum->record);
        }
    } else {
        try {
            bounds = injectors[idx]->SecondaryInjectionBounds(datum->record);
            auto const & w = secondary_process_weighter_maps[idx].at(datum->record.signature.primary_type);
            factors.physical = w->PhysicalProbability(bounds, datum->record);
            factors.generation = w->GenerationProbability(*datum);
            if(with_diagnostics) {
                factors.interaction_prob = w->InteractionProbability(bounds, datum->record);
                factors.position_prob = w->NormalizedPositionProbability(bounds, datum->record);
            }
        } catch(const std::out_of_range& oor) {
            std::ostringstream oss;
            oss << "Weighter::ComputeVertexFactors: no secondary process weighter for secondary type "
                << datum->record.signature.primary_type
                << " in injector " << idx
                << " (" << oor.what() << ") [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
    }
    return factors;
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
            VertexWeightFactors factors = ComputeVertexFactors(idx, datum);
            physical_probability *= factors.physical;
            generation_probability *= factors.generation;
        }
        // Weight-calculation guard: a generation_probability that is <= 0 or
        // non-finite would produce an infinite or NaN weight via the reciprocal
        // below; a non-finite physical probability is equally unrecoverable.
        // Both are configuration/physics failures the caller must not silently
        // absorb.  physical_probability == 0 is NOT caught here: it drives
        // inv_weight to +inf and yields a legitimate 0.0-weight event.
        if(generation_probability <= 0.0 || !std::isfinite(generation_probability)
           || !std::isfinite(physical_probability)) {
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

// Per-vertex decomposition of EventWeight, using the same ComputeVertexFactors
// call per datum so the two paths cannot diverge. A generation-probability
// failure that would make EventWeight throw is instead flagged and yields a
// NaN total; a zero physical probability is flagged and yields a 0.0 total,
// matching EventWeight's non-throwing zero-weight case.
EventWeightBreakdown Weighter::EventWeightWithBreakdown(
        siren::dataclasses::InteractionTree const & tree) const {
    EventWeightBreakdown breakdown;
    double inv_weight = 0;
    bool usable = true;
    for(unsigned int idx = 0; idx < injectors.size(); ++idx) {
        double physical_probability = 1.0;
        double generation_probability = injectors[idx]->InjectedEvents();
        if(injectors[idx]->InjectedEvents() == 0) {
            generation_probability = injectors[idx]->EventsToInject();
        }
        for(auto const & datum : tree.tree) {
            VertexWeightFactors factors = ComputeVertexFactors(idx, datum, true);
            if(!datum->is_root()) {
                factors.depth = static_cast<int>(datum->depth(tree));
            }
            if(factors.generation <= 0.0 || !std::isfinite(factors.generation)) {
                factors.flags.push_back("generation density zero");
            }
            if(factors.physical == 0.0) {
                factors.flags.push_back("outside physical support (weight 0)");
            } else if(!std::isfinite(factors.physical)) {
                factors.flags.push_back("physical density non-finite");
            }
            physical_probability *= factors.physical;
            generation_probability *= factors.generation;
            breakdown.vertices.push_back(std::move(factors));
        }
        if(generation_probability <= 0.0 || !std::isfinite(generation_probability)
           || !std::isfinite(physical_probability)) {
            usable = false;
            continue;
        }
        inv_weight += generation_probability / physical_probability;
    }
    breakdown.total = usable ? (1.0 / inv_weight) : std::numeric_limits<double>::quiet_NaN();
    return breakdown;
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
