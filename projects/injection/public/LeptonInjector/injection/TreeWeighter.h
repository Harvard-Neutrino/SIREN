#pragma once
#ifndef LI_TreeWeighter_H
#define LI_TreeWeighter_H

#include <map>                                           // for map
#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector
#include <utility>                                       // for pair

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/InteractionTree.h"  // for InteractionT...
#include "LeptonInjector/dataclasses/Particle.h"         // for Particle

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace detector { class DetectorModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class WeightableDistribution; } }
namespace LI { namespace injection { class Injector; } }
namespace LI { namespace injection { struct InjectionProcess; } }
namespace LI { namespace injection { struct PhysicalProcess; } }
namespace LI { namespace math { class Vector3D; } }

namespace LI {
namespace injection {

// Class handling weight calculation for a single pair of injection and physical processes
class LeptonProcessWeighter {
private:
    std::shared_ptr<LI::injection::PhysicalProcess> phys_process;
    std::shared_ptr<LI::injection::InjectionProcess> inj_process;
    std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> unique_gen_distributions;
    std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>> unique_phys_distributions;
    std::shared_ptr<LI::detector::DetectorModel> earth_model;
    void Initialize();
    double normalization;
public:
    double InteractionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double NormalizedPositionProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double PhysicalProbability(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionRecord const & record) const;
    double GenerationProbability(LI::dataclasses::InteractionTreeDatum const & datum) const;
    double EventWeight(std::pair<LI::math::Vector3D, LI::math::Vector3D> const & bounds, LI::dataclasses::InteractionTreeDatum const & datum) const;
    LeptonProcessWeighter(std::shared_ptr<LI::injection::PhysicalProcess> phys_process, std::shared_ptr<LI::injection::InjectionProcess> inj_process, std::shared_ptr<LI::detector::DetectorModel> earth_model);

}; // LeptonProcessWeighter

// Parent class for calculating event weights
// Assumes there is a unique secondary physical process for each particle type
class LeptonTreeWeighter {
private:
    // Supplied by constructor
    std::vector<std::shared_ptr<Injector>> injectors;
    std::shared_ptr<LI::detector::DetectorModel> earth_model;
    std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process;
    std::vector<std::shared_ptr<LI::injection::PhysicalProcess>> secondary_physical_processes;

    // Calculated upon initialization
    std::vector<std::shared_ptr<LeptonProcessWeighter>> primary_process_weighters;
    std::vector<
      std::map<
        LI::dataclasses::Particle::ParticleType,
        std::shared_ptr<LeptonProcessWeighter>
      >
    > secondary_process_weighter_maps;

    void Initialize();
public:
    double EventWeight(LI::dataclasses::InteractionTree const & tree) const;
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> earth_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process, std::vector<std::shared_ptr<LI::injection::PhysicalProcess>> secondary_physical_processes);
    LeptonTreeWeighter(std::vector<std::shared_ptr<Injector>> injectors, std::shared_ptr<LI::detector::DetectorModel> earth_model, std::shared_ptr<LI::injection::PhysicalProcess> primary_physical_process);

}; // LeptonTreeWeighter


} //namespace injection
} //namespace LI

CEREAL_CLASS_VERSION(LI::injection::LeptonTreeWeighter, 0);


#endif // LI_TreeWeighter_H
