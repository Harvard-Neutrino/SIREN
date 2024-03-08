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

#include "SIREN/dataclasses/InteractionTree.h"  // for InteractionT...
#include "SIREN/dataclasses/Particle.h"         // for Particle

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class SecondaryInjectionDistribution; } }
namespace SI { namespace distributions { class WeightableDistribution; } }
namespace SI { namespace injection { class Injector; } }
namespace SI { namespace injection { class PrimaryInjectionProcess; } }
namespace SI { namespace injection { class SecondaryInjectionProcess; } }
namespace SI { namespace injection { class PhysicalProcess; } }
namespace SI { namespace math { class Vector3D; } }

namespace SI {
namespace injection {

// Class handling weight calculation for a single pair of injection and physical processes
template<typename ProcessType>
class ProcessWeighter {
private:
    std::shared_ptr<SI::injection::PhysicalProcess> phys_process;
    std::shared_ptr<ProcessType> inj_process;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> unique_gen_distributions;
    std::vector<std::shared_ptr<SI::distributions::WeightableDistribution>> unique_phys_distributions;
    std::shared_ptr<SI::detector::DetectorModel> detector_model;
    std::vector<std::shared_ptr<typename ProcessType::InjectionType>> const & GetInjectionDistributions();
    void Initialize();
    double normalization;
public:
    double InteractionProbability(std::tuple<SI::math::Vector3D, SI::math::Vector3D> const & bounds, SI::dataclasses::InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::tuple<SI::math::Vector3D, SI::math::Vector3D> const & bounds, SI::dataclasses::InteractionRecord const & record) const;
    double PhysicalProbability(std::tuple<SI::math::Vector3D, SI::math::Vector3D> const & bounds, SI::dataclasses::InteractionRecord const & record) const;
    double GenerationProbability(SI::dataclasses::InteractionTreeDatum const & datum) const;
    double EventWeight(std::tuple<SI::math::Vector3D, SI::math::Vector3D> const & bounds, SI::dataclasses::InteractionTreeDatum const & datum) const;
    ProcessWeighter(std::shared_ptr<SI::injection::PhysicalProcess> phys_process, std::shared_ptr<ProcessType> inj_process, std::shared_ptr<SI::detector::DetectorModel> detector_model);

}; // ProcessWeighter

typedef ProcessWeighter<SI::injection::PrimaryInjectionProcess> PrimaryProcessWeighter;
typedef ProcessWeighter<SI::injection::SecondaryInjectionProcess> SecondaryProcessWeighter;

// Parent class for calculating event weights
// Assumes there is a unique secondary physical process for each particle type
class LeptonTreeWeighter {
private:
    // Supplied by constructor
    std::vector<std::shared_ptr<Injector>> injectors;
    std::shared_ptr<SI::detector::DetectorModel> detector_model;
    std::shared_ptr<SI::injection::PhysicalProcess> primary_physical_process;
    std::vector<std::shared_ptr<SI::injection::PhysicalProcess>> secondary_physical_processes;

    // Calculated upon initialization
    std::vector<std::shared_ptr<PrimaryProcessWeighter>> primary_process_weighters;
    std::vector<
      std::map<
        SI::dataclasses::ParticleType,
        std::shared_ptr<SecondaryProcessWeighter>
      >
    > secondary_process_weighter_maps;

    void Initialize();
public:
    double EventWeight(SI::dataclasses::InteractionTree const & tree) const;
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<SI::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<SI::injection::PhysicalProcess>> secondary_physical_processes);
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<SI::injection::PhysicalProcess> primary_physical_process);

}; // LeptonTreeWeighter


} //namespace injection
} //namespace SI

#include "TreeWeighter.tcc"

CEREAL_CLASS_VERSION(SI::injection::LeptonTreeWeighter, 0);


#endif // LI_TreeWeighter_H
