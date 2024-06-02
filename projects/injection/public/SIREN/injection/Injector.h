#pragma once
#ifndef SIREN_Injector_H
#define SIREN_Injector_H

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
#include <cereal/types/memory.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/InteractionRecord.h"  // for Interactio...
#include "SIREN/dataclasses/InteractionTree.h"    // for Interactio...
#include "SIREN/dataclasses/Particle.h"           // for Particle
#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h" // for Secondary...
#include "SIREN/interactions/pyDarkNewsCrossSection.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class VertexPositionDistribution; } }
namespace siren { namespace distributions { class SecondaryVertexPositionDistribution; } }
namespace siren { namespace geometry { class Geometry; } }
namespace siren { namespace injection { class PrimaryInjectionProcess; } }
namespace siren { namespace injection { class SecondaryInjectionProcess; } }
namespace siren { namespace math { class Vector3D; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace injection {

class Injector {
friend cereal::access;
public:
    virtual ~Injector() {};
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<siren::utilities::SIREN_random> random;
    std::shared_ptr<siren::detector::DetectorModel> detector_model;
    // This function returns true if the given secondary index i of the datum should not be simulated
    // Defaults to no secondary interactions being saved
    std::function<bool(std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> stopping_condition= [&](std::shared_ptr<siren::dataclasses::InteractionTreeDatum> datum, size_t i) {
        return true;
    };
    Injector();
private:
    std::shared_ptr<injection::PrimaryInjectionProcess> primary_process;
    std::shared_ptr<distributions::VertexPositionDistribution> primary_position_distribution;
    std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes;
    std::vector<std::shared_ptr<distributions::SecondaryVertexPositionDistribution>> secondary_position_distributions;
    std::map<siren::dataclasses::ParticleType,std::shared_ptr<siren::injection::SecondaryInjectionProcess>> secondary_process_map;
    std::map<siren::dataclasses::ParticleType,std::shared_ptr<distributions::SecondaryVertexPositionDistribution>> secondary_position_distribution_map;
public:
    // Constructors
    Injector(unsigned int events_to_inject, std::string filename, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<siren::utilities::SIREN_random> random);

    void SetStoppingCondition(std::function<bool(std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> f_in) {stopping_condition = f_in;}
    std::shared_ptr<distributions::VertexPositionDistribution> FindPrimaryVertexDistribution(std::shared_ptr<siren::injection::PrimaryInjectionProcess> process);
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> FindSecondaryVertexDistribution(std::shared_ptr<siren::injection::SecondaryInjectionProcess> process);
    void SetPrimaryProcess(std::shared_ptr<siren::injection::PrimaryInjectionProcess> primary);
    std::shared_ptr<siren::injection::PrimaryInjectionProcess> GetPrimaryProcess() {return primary_process;}
    std::vector<std::shared_ptr<siren::injection::SecondaryInjectionProcess>> GetSecondaryProcesses() {return secondary_processes;}
    std::map<siren::dataclasses::ParticleType,std::shared_ptr<siren::injection::SecondaryInjectionProcess>> GetSecondaryProcessMap() {return secondary_process_map;}
    void AddSecondaryProcess(std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary);
    virtual siren::dataclasses::InteractionRecord NewRecord() const; // set primary type from primary process;
    void SetRandom(std::shared_ptr<siren::utilities::SIREN_random> random);
    virtual void SampleCrossSection(siren::dataclasses::InteractionRecord & record) const;
    virtual void SampleCrossSection(siren::dataclasses::InteractionRecord & record,
                                    std::shared_ptr<siren::interactions::InteractionCollection> interactions) const;
    siren::dataclasses::InteractionRecord SampleSecondaryProcess(siren::dataclasses::SecondaryDistributionRecord & secondary_record) const;
    siren::dataclasses::InteractionTree GenerateEvent();
    virtual std::string Name() const;
    virtual double SecondaryGenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum) const;
    virtual double SecondaryGenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum, std::shared_ptr<siren::injection::SecondaryInjectionProcess> process) const;
    virtual double GenerationProbability(siren::dataclasses::InteractionTree const & tree) const;
    virtual double GenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum, std::shared_ptr<siren::injection::PrimaryInjectionProcess> process = NULL) const;
    virtual double GenerationProbability(siren::dataclasses::InteractionRecord const & record, std::shared_ptr<siren::injection::PrimaryInjectionProcess> process = NULL) const;
    virtual std::set<std::vector<std::string>> DensityVariables() const;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> PrimaryInjectionBounds(siren::dataclasses::InteractionRecord const & interaction) const;
    virtual std::tuple<siren::math::Vector3D, siren::math::Vector3D> SecondaryInjectionBounds(siren::dataclasses::InteractionRecord const & interaction) const;
    virtual std::vector<std::shared_ptr<siren::distributions::PrimaryInjectionDistribution>> GetPrimaryInjectionDistributions() const;
    virtual std::shared_ptr<siren::detector::DetectorModel> GetDetectorModel() const;
    virtual std::shared_ptr<siren::interactions::InteractionCollection> GetInteractions() const;
    unsigned int InjectedEvents() const;
    unsigned int EventsToInject() const;
    void ResetInjectedEvents();
    operator bool() const;
    void SaveInjector(std::string const & filename) const;
    void LoadInjector(std::string const & filename);

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            // archive(::cereal::make_nvp("SIRENRandom", random));
            // std::cout << "saved SIRENRandom\n";
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("Injector only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            std::shared_ptr<injection::PrimaryInjectionProcess> _primary_process;
            std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> _secondary_processes;

            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            // archive(::cereal::make_nvp("SIRENRandom", random));
            // std::cout << "loaded SIRENRandom\n";
            archive(::cereal::make_nvp("PrimaryProcess", _primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", _secondary_processes));
            SetPrimaryProcess(_primary_process);
            for(auto secondary_process : _secondary_processes) {
                AddSecondaryProcess(secondary_process);
            }
        } else {
            throw std::runtime_error("Injector only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::Injector, 0);

#endif // SIREN_Injector_H

