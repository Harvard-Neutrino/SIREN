#pragma once
#ifndef LI_LeptonInjector_H
#define LI_LeptonInjector_H

#include <map>                                             // for map
#include <set>                                             // for set
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

#include "LeptonInjector/dataclasses/InteractionRecord.h"  // for Interactio...
#include "LeptonInjector/dataclasses/InteractionTree.h"    // for Interactio...
#include "LeptonInjector/dataclasses/Particle.h"           // for Particle

namespace LI { namespace crosssections { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct DecayRecord; } }
namespace LI { namespace detector { class EarthModel; } }
namespace LI { namespace distributions { class InjectionDistribution; } }
namespace LI { namespace distributions { class VertexPositionDistribution; } }
namespace LI { namespace geometry { class Geometry; } }
namespace LI { namespace injection { struct InjectionProcess; } }
namespace LI { namespace math { class Vector3D; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace injection {

class InjectorBase {
friend cereal::access;
public:
    virtual ~InjectorBase() {};
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<LI::utilities::LI_random> random;
    std::shared_ptr<LI::detector::EarthModel> earth_model;
    // This funciton returns true if the given datum is the last entry to be saved in a tree
    std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> stopping_condition;
    InjectorBase();
private:
    std::shared_ptr<injection::InjectionProcess> primary_process;
    std::shared_ptr<distributions::VertexPositionDistribution> primary_position_distribution;
    std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes;
    std::vector<std::shared_ptr<distributions::VertexPositionDistribution>> secondary_position_distributions;
    std::map<LI::dataclasses::Particle::ParticleType,std::shared_ptr<LI::injection::InjectionProcess>> secondary_process_map;
    std::map<LI::dataclasses::Particle::ParticleType,std::shared_ptr<distributions::VertexPositionDistribution>> secondary_position_distribution_map;
public:
    // Constructors
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<LI::utilities::LI_random> random);
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::shared_ptr<LI::utilities::LI_random> random);
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes, std::shared_ptr<LI::utilities::LI_random> random);

    void SetStoppingCondition(std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> f_in) {stopping_condition = f_in;}
    std::shared_ptr<distributions::VertexPositionDistribution> FindPositionDistribution(std::shared_ptr<LI::injection::InjectionProcess> process);
    void SetPrimaryProcess(std::shared_ptr<LI::injection::InjectionProcess> primary);
    std::shared_ptr<LI::injection::InjectionProcess> GetPrimaryProcess() {return primary_process;}
    std::vector<std::shared_ptr<LI::injection::InjectionProcess>> GetSecondaryProcesses() {return secondary_processes;}
    std::map<LI::dataclasses::Particle::ParticleType,std::shared_ptr<LI::injection::InjectionProcess>> GetSecondaryProcessMap() {return secondary_process_map;}
    void AddSecondaryProcess(std::shared_ptr<LI::injection::InjectionProcess> secondary);
    virtual LI::dataclasses::InteractionRecord NewRecord() const; // set primary type from primary process;
    void SetRandom(std::shared_ptr<LI::utilities::LI_random> random);
    virtual void SampleCrossSection(LI::dataclasses::InteractionRecord & record) const;
    virtual void SampleCrossSection(LI::dataclasses::InteractionRecord & record,
                                    std::shared_ptr<LI::crosssections::InteractionCollection> cross_sections) const;
    virtual void SampleNeutrissimoDecay(LI::dataclasses::InteractionRecord const & interaction, LI::dataclasses::DecayRecord & decay, double width, double alpha_gen, double alpha_phys, LI::geometry::Geometry *fiducial, double buffer) const;
    virtual void SamplePairProduction(LI::dataclasses::DecayRecord const & decay, LI::dataclasses::InteractionRecord & pairprod) const;
    bool SampleSecondaryProcess(unsigned int idx,
                                std::shared_ptr<LI::dataclasses::InteractionTreeDatum> parent,
                                LI::dataclasses::InteractionTreeDatum & datum);
    LI::dataclasses::InteractionTree GenerateEvent();
    virtual std::string Name() const;
    virtual double SecondaryGenerationProbability(std::shared_ptr<LI::dataclasses::InteractionTreeDatum> const & datum) const;
    virtual double GenerationProbability(LI::dataclasses::InteractionTree const & tree) const;
    virtual double GenerationProbability(std::shared_ptr<LI::dataclasses::InteractionTreeDatum> const & datum, std::shared_ptr<LI::injection::InjectionProcess> process = NULL) const;
    virtual double GenerationProbability(LI::dataclasses::InteractionRecord const & record, std::shared_ptr<LI::injection::InjectionProcess> process = NULL) const;
    virtual std::set<std::vector<std::string>> DensityVariables() const;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionTreeDatum const & datum, LI::dataclasses::Particle::ParticleType const & primary_type) const;
    virtual std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> GetInjectionDistributions() const;
    virtual std::shared_ptr<LI::detector::EarthModel> GetEarthModel() const;
    virtual std::shared_ptr<LI::crosssections::InteractionCollection> GetCrossSections() const;
    unsigned int InjectedEvents() const;
    unsigned int EventsToInject() const;
    operator bool() const;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            //archive(::cereal::make_nvp("StoppingCondition", stopping_condition));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            //archive(::cereal::make_nvp("StoppingCondition", stopping_condition));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::InjectorBase, 0);

#endif // LI_LeptonInjector_H

