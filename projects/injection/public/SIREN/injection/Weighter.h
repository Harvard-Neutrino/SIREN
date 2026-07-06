#pragma once
#ifndef SIREN_Weighter_H
#define SIREN_Weighter_H

#include <map>                                           // for map
#include <string>                                        // for string
#include <tuple>
#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/InteractionTree.h"  // for InteractionT...
#include "SIREN/dataclasses/Particle.h"         // for Particle

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class SecondaryInjectionDistribution; } }
namespace siren { namespace distributions { class WeightableDistribution; } }
namespace siren { namespace injection { class Injector; } }
namespace siren { namespace injection { class PrimaryInjectionProcess; } }
namespace siren { namespace injection { class SecondaryInjectionProcess; } }
namespace siren { namespace injection { class PhysicalProcess; } }
namespace siren { namespace math { class Vector3D; } }

namespace siren {
namespace injection {

// Class handling weight calculation for a single pair of injection and physical processes
template<typename ProcessType>
class ProcessWeighter {
private:
    std::shared_ptr<siren::injection::PhysicalProcess> phys_process;
    std::shared_ptr<ProcessType> inj_process;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> unique_gen_distributions;
    std::vector<std::shared_ptr<siren::distributions::WeightableDistribution>> unique_phys_distributions;
    std::shared_ptr<siren::detector::DetectorModel> detector_model;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> const & GetInjectionDistributions();
    void Initialize();
    double normalization;
public:
    double InteractionProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const;
    double SurvivalProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const;
    double PhysicalProbability(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionRecord const & record) const;
    double GenerationProbability(siren::dataclasses::InteractionTreeDatum const & datum) const;
    double EventWeight(std::tuple<siren::math::Vector3D, siren::math::Vector3D> const & bounds, siren::dataclasses::InteractionTreeDatum const & datum) const;
    ProcessWeighter(std::shared_ptr<siren::injection::PhysicalProcess> phys_process, std::shared_ptr<ProcessType> inj_process, std::shared_ptr<siren::detector::DetectorModel> detector_model);

}; // ProcessWeighter

typedef ProcessWeighter<siren::injection::PrimaryInjectionProcess> PrimaryProcessWeighter;
typedef ProcessWeighter<siren::injection::SecondaryInjectionProcess> SecondaryProcessWeighter;

// One tree vertex's weight factors, recorded for diagnostics. generation and
// physical are the exact per-datum factors the event weight consumes; the other
// fields are observation-only and never fed back into the weight.
struct VertexWeightFactors {
    int injector_index = 0;
    int depth = 0;
    int primary_pdg = 0;
    double generation = 1.0;
    double physical = 1.0;
    double interaction_prob = 1.0;
    double position_prob = 1.0;
    std::map<std::string, double> channel_densities;
    std::vector<std::string> cancelled;
    std::vector<std::string> flags;
};

// Per-vertex decomposition of an event weight. total equals
// Weighter::EventWeight(tree): a zero-generation vertex (the path EventWeight
// throws on) is flagged and yields a NaN total, while a zero-physical vertex is
// flagged and yields a 0.0 total, matching EventWeight's zero-weight case.
struct EventWeightBreakdown {
    double total = 0.0;
    std::vector<VertexWeightFactors> vertices;
};

// Parent class for calculating event weights
// Assumes there is a unique secondary physical process for each particle type
class Weighter {
private:
    // Supplied by constructor
    std::vector<std::shared_ptr<Injector>> injectors;
    std::shared_ptr<siren::detector::DetectorModel> detector_model;
    std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process;
    std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> secondary_physical_processes;

    // Calculated upon initialization
    std::vector<std::shared_ptr<PrimaryProcessWeighter>> primary_process_weighters;
    std::vector<
      std::map<
        siren::dataclasses::ParticleType,
        std::shared_ptr<SecondaryProcessWeighter>
      >
    > secondary_process_weighter_maps;

    void Initialize();
    // with_diagnostics fills the observation-only interaction_prob/position_prob
    // fields; EventWeight leaves it false so it computes only the physical and
    // generation factors it consumes.
    VertexWeightFactors ComputeVertexFactors(unsigned int idx,
        std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum,
        bool with_diagnostics = false) const;
public:
    double EventWeight(siren::dataclasses::InteractionTree const & tree) const;
    EventWeightBreakdown EventWeightWithBreakdown(siren::dataclasses::InteractionTree const & tree) const;
    std::vector<std::shared_ptr<Injector>> const & GetInjectors() const;
    std::shared_ptr<siren::detector::DetectorModel> GetDetectorModel() const;
    std::shared_ptr<siren::injection::PhysicalProcess> GetPrimaryPhysicalProcess() const;
    std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> const & GetSecondaryPhysicalProcesses() const;
    std::vector<double> GetInteractionProbabilities(siren::dataclasses::InteractionTree const & tree, int i_inj = 0) const;
    std::vector<double> GetSurvivalProbabilities(siren::dataclasses::InteractionTree const & tree, int i_inj = 0) const;
    Weighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> secondary_physical_processes);
    Weighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process);
    Weighter(std::vector<std::shared_ptr<Injector>> injectors, std::string filename);
    void SaveWeighter(std::string const & filename) const;
    void LoadWeighter(std::string const & filename);

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Injectors", injectors));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            archive(::cereal::make_nvp("PrimaryPhysicalProcess", primary_physical_process));
            archive(::cereal::make_nvp("SecondaryPhysicalProcesses", secondary_physical_processes));
        } else {
            throw std::runtime_error("Weighter only supports version <= 0!");
        }
    }

    template<typename Archive>
    static void load_and_construct(Archive & archive, cereal::construct<Weighter> & construct, std::uint32_t const version) {
        if(version == 0) {
            std::vector<std::shared_ptr<Injector>> injectors;
            std::shared_ptr<siren::detector::DetectorModel> detector_model;
            std::shared_ptr<siren::injection::PhysicalProcess> primary_physical_process;
            std::vector<std::shared_ptr<siren::injection::PhysicalProcess>> secondary_physical_processes;
            archive(::cereal::make_nvp("Injectors", injectors));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            archive(::cereal::make_nvp("PrimaryPhysicalProcess", primary_physical_process));
            archive(::cereal::make_nvp("SecondaryPhysicalProcesses", secondary_physical_processes));
            construct(injectors, detector_model, primary_physical_process, secondary_physical_processes);
        } else {
            throw std::runtime_error("Weighter only supports version <= 0!");
        }
    }

}; // Weighter


} //namespace injection
} //namespace siren

#include "Weighter.tcc"

CEREAL_CLASS_VERSION(siren::injection::Weighter, 0);


#endif // SIREN_Weighter_H
