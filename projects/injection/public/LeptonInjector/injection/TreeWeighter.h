#pragma once
#ifndef LI_TreeWeighter_H
#define LI_TreeWeighter_H

#include <map>                                           // for map
#include <tuple>
#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/InteractionTree.h"  // for InteractionT...
#include "LeptonInjector/dataclasses/Particle.h"         // for Particle

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace LI { namespace distributions { class SecondaryInjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace injection { class Injector; } }
namespace LI { namespace injection { class PrimaryInjectionProcess; } }
namespace LI { namespace injection { class SecondaryInjectionProcess; } }
namespace LI { namespace injection { class PhysicalProcess; } }
namespace LI { namespace math { class Vector3D; } }

namespace LI {
namespace injection {

// Class handling weight calculation for a single pair of injection and physical processes
template<typename ProcessType>
class ProcessWeighter {
private:
    std::shared_ptr<LI::injection::PhysicalProcess> phys_process;
    std::shared_ptr<ProcessType> inj_process;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> unique_gen_distributions;
    std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>> unique_phys_distributions;
    std::shared_ptr<LI::detector::DetectorModel> detector_model;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> const & GetInjectionDistributions();
    void Initialize();
    double normalization;
public:
    double InteractionProbability(std::tuple<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::tuple<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double PhysicalProbability(std::tuple<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double GenerationProbability(LI::dataclasses::InteractionTreeDatum const & datum) const;
    double EventWeight(std::tuple<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionTreeDatum const & datum) const;
    ProcessWeighter(std::shared_ptr<LI::injection::PhysicalProcess> phys_process, std::shared_ptr<ProcessType> inj_process, std::shared_ptr<LI::detector::DetectorModel> detector_model);

}; // ProcessWeighter

typedef ProcessWeighter<LI::injection::PrimaryInjectionProcess> PrimaryProcessWeighter;
typedef ProcessWeighter<LI::injection::SecondaryInjectionProcess> SecondaryProcessWeighter;

// Parent class for calculating event weights
// Assumes there is a unique secondary physical process for each particle type
class LeptonTreeWeighter {
private:
    // Supplied by constructor
    std::vector<std::shared_ptr<Injector>> injectors;
    std::shared_ptr<LI::detector::DetectorModel> detector_model;
    std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process;
    std::vector<std::shared_ptr<LI::injection::PhysicalProcess>> secondary_physical_processes;

    // Calculated upon initialization
    std::vector<std::shared_ptr<PrimaryProcessWeighter>> primary_process_weighters;
    std::vector<
      std::map<
        LI::dataclasses::Particle::ParticleType,
        std::shared_ptr<SecondaryProcessWeighter>
      >
    > secondary_process_weighter_maps;

    void Initialize();
public:
    double EventWeight(LI::dataclasses::InteractionTree const & tree) const;
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<LI::injection::PhysicalProcess>> secondary_physical_processes);
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process);

}; // LeptonTreeWeighter


} //namespace injection
} //namespace LI

#include "TreeWeighter.tcc"

CEREAL_CLASS_VERSION(LI::injection::LeptonTreeWeighter, 0);


#endif // LI_TreeWeighter_H
