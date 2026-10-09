#include "SIREN/injection/Weighter.h"

#include <iterator>                                              // for ite...
#include <array>                                                  // for array
#include <cassert>                                                // for assert
#include <cmath>                                                  // for exp, isfinite
#include <cstdint>                                                // for uint32_t
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

template<typename ProcessPtr>
static void RecordMixtureDiagnostics(VertexWeightFactors & factors,
        ProcessPtr const & inj_process,
        std::vector<std::string> const & cancelled_names,
        siren::dataclasses::InteractionRecord const & record,
        std::shared_ptr<siren::detector::DetectorModel> const & detector_model) {
    factors.cancelled = cancelled_names;
    if(inj_process && inj_process->HasPhaseSpace(record.signature)) {
        auto ps = inj_process->GetPhaseSpace(record.signature);
        factors.channel_density_topology = ps->CommonTopology();
        factors.channel_density_measure = ps->CommonMeasure();
        std::vector<double> contributions = ps->DensityBreakdown(detector_model, record);
        for(std::size_t c = 0; c < contributions.size(); ++c) {
            factors.channel_densities["channel[" + std::to_string(c) + "]"] = contributions[c];
        }
    }
}

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
            RecordMixtureDiagnostics(factors,
                injectors[idx]->GetPrimaryProcess(),
                primary_process_weighters[idx]->GetCancelledDistributionNames(),
                datum->record, detector_model);
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
                RecordMixtureDiagnostics(factors,
                    injectors[idx]->GetSecondaryProcessMap().at(datum->record.signature.primary_type),
                    w->GetCancelledDistributionNames(),
                    datum->record, detector_model);
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

EventWeightBreakdown Weighter::PoolEventWeight(
        siren::dataclasses::InteractionTree const & tree,
        bool with_diagnostics,
        std::string * problem) const {
    // The weight pools every injector i that could have produced the tree:
    //
    //   w = 1 / sum_i [ N_i * prod_d p_gen(i, d) / prod_d p_phys(i, d) ]
    //
    // N_i is injector i's attempt count and d runs over the tree's vertices.
    // The physical interaction and position probabilities depend on each
    // injector's bounds. The position density's normalization equals the
    // interaction probability, so their product is the unnormalized position
    // density. An injector whose generation density vanishes at some vertex
    // cannot produce the tree and contributes nothing; the weight is undefined
    // only if no injector can produce it.
    EventWeightBreakdown breakdown;
    std::string first_problem;
    auto reject = [&first_problem](std::string const & reason) {
        if(first_problem.empty())
            first_problem = reason;
    };

    if(tree.tree.empty())
        reject("cannot weight an empty or failed event");

    double inv_weight = 0.0;
    bool zero_weight = false;
    bool any_generating = false;
    for(unsigned int idx = 0; idx < injectors.size() && !tree.tree.empty(); ++idx) {
        std::string const injector_name = "injector " + std::to_string(idx);
        // Before generation starts, the configured event count stands in for
        // the attempt count.
        double generation = injectors[idx]->InjectionAttempts();
        if(generation == 0)
            generation = injectors[idx]->EventsToInject();
        if(generation == 0)
            reject(injector_name + " has no attempts and no events to inject");
        double physical = 1.0;
        bool zero_generation_density = false;
        for(auto const & datum : tree.tree) {
            VertexWeightFactors factors = ComputeVertexFactors(idx, datum, with_diagnostics);
            if(!datum->is_root())
                factors.depth = static_cast<int>(datum->depth(tree));
            std::ostringstream where;
            where << injector_name << " at depth " << factors.depth
                  << " (pdg " << factors.vertex_pdg << ")";
            if(!std::isfinite(factors.generation) || factors.generation < 0.0) {
                std::string flag = std::isfinite(factors.generation)
                    ? "generation density negative" : "generation density non-finite";
                factors.flags.push_back(flag);
                reject(flag + " for " + where.str() + ": " + std::to_string(factors.generation));
            } else if(factors.generation == 0.0) {
                factors.flags.push_back("generation density zero");
                zero_generation_density = true;
            }
            if(!std::isfinite(factors.physical) || factors.physical < 0.0) {
                std::string flag = std::isfinite(factors.physical)
                    ? "physical density negative" : "physical density non-finite";
                factors.flags.push_back(flag);
                reject(flag + " for " + where.str() + ": " + std::to_string(factors.physical));
            } else if(factors.physical == 0.0) {
                factors.flags.push_back("outside physical support (weight 0)");
            }
            generation *= factors.generation;
            physical *= factors.physical;
            breakdown.vertices.push_back(std::move(factors));
        }
        std::ostringstream values;
        values << ": generation_probability=" << generation
               << ", physical_probability=" << physical;
        if(!std::isfinite(generation) || !std::isfinite(physical)
           || (generation == 0.0 && !zero_generation_density)) {
            // A product of finite positive factors overflowed or underflowed.
            breakdown.vertices.back().flags.push_back("unusable event probabilities");
            reject("unusable probabilities for " + injector_name + values.str());
            continue;
        }
        if(generation == 0.0)
            continue;
        any_generating = true;
        if(physical == 0.0) {
            zero_weight = true;
            continue;
        }
        inv_weight += generation / physical;
        if(!std::isfinite(inv_weight)) {
            breakdown.vertices.back().flags.push_back("inverse weight overflow");
            reject("inverse weight overflow at " + injector_name + values.str());
        }
    }
    if(!tree.tree.empty() && !any_generating) {
        if(!breakdown.vertices.empty())
            breakdown.vertices.back().flags.push_back("no injector can produce this event");
        reject("no injector can produce this event");
    }

    double weight = zero_weight ? 0.0 : 1.0 / inv_weight;
    if(first_problem.empty() && (!std::isfinite(weight) || weight < 0.0)) {
        breakdown.vertices.back().flags.push_back("unusable event weight");
        reject("unusable event weight " + std::to_string(weight)
               + " from inverse weight " + std::to_string(inv_weight));
    }
    breakdown.total = first_problem.empty()
        ? weight : std::numeric_limits<double>::quiet_NaN();
    if(problem)
        *problem = first_problem;
    return breakdown;
}

double Weighter::EventWeight(siren::dataclasses::InteractionTree const & tree) const {
    std::string problem;
    EventWeightBreakdown breakdown = PoolEventWeight(tree, false, &problem);
    if(!problem.empty()) {
        throw siren::utilities::WeightCalculationError(
            "Weighter::EventWeight: " + problem + " [siren-docs: errors#weight-calc]");
    }
    return breakdown.total;
}

EventWeightBreakdown Weighter::EventWeightWithBreakdown(
        siren::dataclasses::InteractionTree const & tree) const {
    return PoolEventWeight(tree, true, nullptr);
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

namespace {
// Header word marking a version-stamped weighter archive; headerless archives
// begin directly with the Injectors payload.
constexpr std::uint32_t kWeighterArchiveMagic = 0x53575447; // "SWGT"
} // anonymous namespace

void Weighter::SaveWeighter(std::string const & filename) const {
    std::string const path = filename + ".siren_weighter";
    std::ofstream os(path, std::ios::binary);
    if(!os) {
        throw std::runtime_error(
            "Failed to open weighter archive '" + path + "' for writing");
    }
    ::cereal::BinaryOutputArchive archive(os);
    std::uint32_t magic = kWeighterArchiveMagic;
    // ::cereal::detail::Version<Weighter> is what CEREAL_CLASS_VERSION registered,
    // so the header version can never drift from the class version.
    std::uint32_t version = ::cereal::detail::Version<Weighter>::version;
    archive(magic, version);
    this->save(archive, version);
}

void Weighter::LoadWeighter(std::string const & filename) {
    std::string const path = filename + ".siren_weighter";
    // A missing file otherwise surfaces as a cryptic cereal stream error; name it.
    {
        std::ifstream is(path, std::ios::binary);
        if(!is) {
            throw std::runtime_error(
                "Failed to load weighter archive '" + path + "': cannot open file");
        }
    }
    {
        std::ifstream is(path, std::ios::binary);
        ::cereal::BinaryInputArchive archive(is);
        std::uint32_t magic = 0;
        try {
            archive(magic);
        } catch(...) {
            // Too short to hold a magic word; fall through to the headerless path.
            magic = 0;
        }
        if(magic == kWeighterArchiveMagic) {
            try {
                std::uint32_t version = 0;
                archive(version);
                Weighter temp;
                temp.load(archive, version);
                *this = std::move(temp);
            } catch(std::exception const & e) {
                throw std::runtime_error(
                    "Failed to load weighter archive '" + path
                    + "': the headered parse failed: " + e.what());
            }
            Initialize();
            return;
        }
    }
    // Headerless (legacy version-0) archive: the first word was the Injectors
    // payload, not a magic word. Reparse the whole stream from the start as the
    // version-0 schema, again into a temporary.
    try {
        std::ifstream is(path, std::ios::binary);
        ::cereal::BinaryInputArchive archive(is);
        Weighter temp;
        temp.load(archive, 0);
        *this = std::move(temp);
    } catch(std::exception const & e) {
        throw std::runtime_error(
            "Failed to load weighter archive '" + path
            + "': not a headered archive and the headerless version-0 parse "
              "failed: " + e.what());
    }
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
