#pragma once
#ifndef LI_SIREN_H
#define LI_SIREN_H

#include <map>                                             // for map
#include <set>                                             // for set
#include <tuple>
#include <string>
#include <memory>
#include <vector>                                          // for vector
#include <cstdint>                                         // for uint32_t
#include <utility>
#include <stddef.h>                                        // for NULL
#include <stdexcept>                                       // for runtime_error
#include <functional>
#include <cereal/cereal.hpp>                               // for make_nvp

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

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/dataclasses/InteractionTree.h"    // for Interactio...
#include "SIREN/dataclasses/Particle.h"           // for Particle

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace detector { class DetectorModel; } }
namespace SI { namespace distributions { class PrimaryInjectionDistribution; } }
namespace SI { namespace distributions { class VertexPositionDistribution; } }
namespace SI { namespace distributions { class SecondaryVertexPositionDistribution; } }
namespace SI { namespace geometry { class Geometry; } }
namespace SI { namespace injection { class PrimaryInjectionProcess; } }
namespace SI { namespace injection { class SecondaryInjectionProcess; } }
namespace SI { namespace math { class Vector3D; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace injection {

class Injector {
friend cereal::access;
public:
    virtual ~Injector() {};
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<SI::utilities::LI_random> random;
    std::shared_ptr<SI::detector::DetectorModel> detector_model;
    // This funciton returns true if the given datum is the last entry to be saved in a tree
    std::function<bool(std::shared_ptr<SI::dataclasses::InteractionTreeDatum>, size_t)> stopping_condition;
    Injector();
private:
    std::shared_ptr<injection::PrimaryInjectionProcess> primary_process;
    std::shared_ptr<distributions::VertexPositionDistribution> primary_position_distribution;
    std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes;
    std::vector<std::shared_ptr<distributions::SecondaryVertexPositionDistribution>> secondary_position_distributions;
    std::map<SI::dataclasses::ParticleType,std::shared_ptr<SI::injection::SecondaryInjectionProcess>> secondary_process_map;
    std::map<SI::dataclasses::ParticleType,std::shared_ptr<distributions::SecondaryVertexPositionDistribution>> secondary_position_distribution_map;
public:
    // Constructors
    Injector(unsigned int events_to_inject, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<SI::utilities::LI_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::shared_ptr<SI::utilities::LI_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<SI::utilities::LI_random> random);

    void SetStoppingCondition(std::function<bool(std::shared_ptr<SI::dataclasses::InteractionTreeDatum>, size_t)> f_in) {stopping_condition = f_in;}
    std::shared_ptr<distributions::VertexPositionDistribution> FindPrimaryVertexDistribution(std::shared_ptr<SI::injection::PrimaryInjectionProcess> process);
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> FindSecondaryVertexDistribution(std::shared_ptr<SI::injection::SecondaryInjectionProcess> process);
    void SetPrimaryProcess(std::shared_ptr<SI::injection::PrimaryInjectionProcess> primary);
    std::shared_ptr<SI::injection::PrimaryInjectionProcess> GetPrimaryProcess() {return primary_process;}
    std::vector<std::shared_ptr<SI::injection::SecondaryInjectionProcess>> GetSecondaryProcesses() {return secondary_processes;}
    std::map<SI::dataclasses::ParticleType,std::shared_ptr<SI::injection::SecondaryInjectionProcess>> GetSecondaryProcessMap() {return secondary_process_map;}
    void AddSecondaryProcess(std::shared_ptr<SI::injection::SecondaryInjectionProcess> secondary);
    virtual SI::dataclasses::InteractionRecord NewRecord() const; // set primary type from primary process;
    void SetRandom(std::shared_ptr<SI::utilities::LI_random> random);
    virtual void SampleCrossSection(SI::dataclasses::InteractionRecord & record) const;
    virtual void SampleCrossSection(SI::dataclasses::InteractionRecord & record,
                                    std::shared_ptr<SI::interactions::InteractionCollection> interactions) const;
    SI::dataclasses::InteractionRecord SampleSecondaryProcess(SI::dataclasses::SecondaryDistributionRecord & secondary_record) const;
    SI::dataclasses::InteractionTree GenerateEvent();
    virtual std::string Name() const;
    virtual double SecondaryGenerationProbability(std::shared_ptr<SI::dataclasses::InteractionTreeDatum> const & datum) const;
    virtual double SecondaryGenerationProbability(std::shared_ptr<SI::dataclasses::InteractionTreeDatum> const & datum, std::shared_ptr<SI::injection::SecondaryInjectionProcess> process) const;
    virtual double GenerationProbability(SI::dataclasses::InteractionTree const & tree) const;
    virtual double GenerationProbability(std::shared_ptr<SI::dataclasses::InteractionTreeDatum> const & datum, std::shared_ptr<SI::injection::PrimaryInjectionProcess> process = NULL) const;
    virtual double GenerationProbability(SI::dataclasses::InteractionRecord const & record, std::shared_ptr<SI::injection::PrimaryInjectionProcess> process = NULL) const;
    virtual std::set<std::vector<std::string>> DensityVariables() const;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> SecondaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const;
    virtual std::vector<std::shared_ptr<SI::distributions::PrimaryInjectionDistribution>> GetPrimaryInjectionDistributions() const;
    virtual std::shared_ptr<SI::detector::DetectorModel> GetDetectorModel() const;
    virtual std::shared_ptr<SI::interactions::InteractionCollection> GetInteractions() const;
    unsigned int InjectedEvents() const;
    unsigned int EventsToInject() const;
    operator bool() const;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            //archive(::cereal::make_nvp("StoppingCondition", stopping_condition));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("Injector only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            //archive(::cereal::make_nvp("StoppingCondition", stopping_condition));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("Injector only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace SI

CEREAL_CLASS_VERSION(SI::injection::Injector, 0);

#endif // LI_SIREN_H

