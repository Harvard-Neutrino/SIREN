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
#include "SIREN/injection/FailureLedger.h"
#include "SIREN/interactions/pyDarkNewsCrossSection.h"

namespace siren { namespace interactions { class Interaction; } }
namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace distributions { class PrimaryInjectionDistribution; } }
namespace siren { namespace distributions { class VertexPositionDistribution; } }
namespace siren { namespace distributions { class SecondaryVertexPositionDistribution; } }
namespace siren { namespace geometry { class Geometry; } }
namespace siren { namespace injection { class PrimaryInjectionProcess; } }
namespace siren { namespace injection { class SecondaryInjectionProcess; } }
namespace siren { namespace injection { struct MultiChannelPhaseSpace; } }
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
    unsigned int injection_attempts = 0;
    unsigned int injected_events = 0;
    unsigned int failed_events = 0;
    unsigned int unregistered_secondary_count_ = 0;
    FailureLedger failure_ledger_;
    std::string last_failure_message_;
    siren::dataclasses::InteractionTree last_failed_tree_;
    std::shared_ptr<siren::utilities::SIREN_random> random;
    std::shared_ptr<siren::detector::DetectorModel> detector_model;
    // An unset callback preserves the default of stopping all secondaries.
    std::function<bool(siren::dataclasses::InteractionTree const &,
        std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> stopping_condition;
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
    // Loads an injector archive. An archive that stores RNG state continues
    // that stream; `random` is used only for archives written before then.
    Injector(unsigned int events_to_inject, std::string filename, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::shared_ptr<siren::utilities::SIREN_random> random);
    Injector(unsigned int events_to_inject, std::shared_ptr<siren::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<siren::utilities::SIREN_random> random);

    void SetStoppingCondition(std::function<bool(siren::dataclasses::InteractionTree const &, std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> f_in) {stopping_condition = f_in;}
    std::function<bool(siren::dataclasses::InteractionTree const &, std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> GetStoppingCondition() {return stopping_condition;}
    std::shared_ptr<distributions::VertexPositionDistribution> FindPrimaryVertexDistribution(std::shared_ptr<siren::injection::PrimaryInjectionProcess> process);
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> FindSecondaryVertexDistribution(std::shared_ptr<siren::injection::SecondaryInjectionProcess> process);
    void SetPrimaryProcess(std::shared_ptr<siren::injection::PrimaryInjectionProcess> primary);
    std::shared_ptr<siren::injection::PrimaryInjectionProcess> GetPrimaryProcess() {return primary_process;}
    std::vector<std::shared_ptr<siren::injection::SecondaryInjectionProcess>> GetSecondaryProcesses() {return secondary_processes;}
    std::map<siren::dataclasses::ParticleType,std::shared_ptr<siren::injection::SecondaryInjectionProcess>> GetSecondaryProcessMap() {return secondary_process_map;}
    void AddSecondaryProcess(std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary);
    void SetSecondaryProcesses(std::vector<std::shared_ptr<siren::injection::SecondaryInjectionProcess>> secondary_processes);
    virtual siren::dataclasses::InteractionRecord NewRecord() const; // set primary type from primary process;
    void SetRandom(std::shared_ptr<siren::utilities::SIREN_random> random);
    std::shared_ptr<siren::utilities::SIREN_random> GetRandom() const;
    virtual void SampleCrossSection(siren::dataclasses::InteractionRecord & record) const;
    virtual void SampleCrossSection(siren::dataclasses::InteractionRecord & record,
                                    std::shared_ptr<siren::interactions::InteractionCollection> interactions) const;

    // Select one concrete (interaction model, signature) entry proportional
    // to its rate. Sets record.signature and record.target_mass and returns
    // the selected model for exact downstream dispatch.
    virtual std::shared_ptr<siren::interactions::Interaction> SelectChannel(
        siren::dataclasses::InteractionRecord & record,
        std::shared_ptr<siren::interactions::InteractionCollection> interactions) const;

    // Sample from the concrete model returned by SelectChannel. Used as the
    // fallback when no PhaseSpaceChannel is registered for the signature.
    virtual void SampleMatchingFinalState(
        siren::dataclasses::InteractionRecord & record,
        std::shared_ptr<siren::interactions::Interaction> selected_interaction) const;

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
    virtual void SetDetectorModel(std::shared_ptr<siren::detector::DetectorModel> detector_model);
    virtual std::shared_ptr<siren::interactions::InteractionCollection> GetInteractions() const;
    unsigned int InjectedEvents() const;
    unsigned int InjectionAttempts() const;
    unsigned int EventsToInject() const;
    unsigned int FailedEvents() const;
    unsigned int UnregisteredSecondaryCount() const;
    std::string GetLastFailureMessage() const;
    siren::dataclasses::InteractionTree const & GetLastFailedTree() const;
    FailureLedger const & GetFailureLedger() const;
    void ResetInjectedEvents(unsigned int events_to_inject);
    void ResetInjectedEvents();

    operator bool() const;
    void SaveInjector(std::string const & filename) const;
    void LoadInjector(std::string const & filename);

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (stopping_condition) {
            throw std::runtime_error(
                "Injector archives cannot preserve a stopping-condition callback");
        }
        if(version <= 1) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectionAttempts", injection_attempts));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            if(version >= 1) {
                archive(::cereal::make_nvp("FailedEvents", failed_events));
            }
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            if(version >= 1) {
                archive(::cereal::make_nvp("SIRENRandom", random));
            }
            archive(::cereal::make_nvp("PrimaryProcess", primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", secondary_processes));
        } else {
            throw std::runtime_error("Injector only supports version <= 1!");
        }
    }

    // Rebuilds processes via SetPrimaryProcess/AddSecondaryProcess, so a corrupt or
    // incompatible archive throws siren::utilities::AddProcessFailure instead of exiting.
    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version <= 1) {
            std::shared_ptr<injection::PrimaryInjectionProcess> _primary_process;
            std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> _secondary_processes;

            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectionAttempts", injection_attempts));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            // Version 0 archives have no failed-event count or RNG state;
            // LoadInjector then keeps the caller's engine.
            if(version >= 1) {
                archive(::cereal::make_nvp("FailedEvents", failed_events));
            }
            archive(::cereal::make_nvp("DetectorModel", detector_model));
            if(version >= 1) {
                archive(::cereal::make_nvp("SIRENRandom", random));
            }
            archive(::cereal::make_nvp("PrimaryProcess", _primary_process));
            archive(::cereal::make_nvp("SecondaryProcesses", _secondary_processes));
            SetPrimaryProcess(_primary_process);
            for(auto secondary_process : _secondary_processes) {
                AddSecondaryProcess(secondary_process);
            }
        } else {
            throw std::runtime_error("Injector only supports version <= 1!");
        }
    }
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::Injector, 1);

#endif // SIREN_Injector_H
