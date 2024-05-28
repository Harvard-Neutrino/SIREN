#pragma once
#ifndef SIREN_Process_H
#define SIREN_Process_H

#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"         // for Particle
#include "SIREN/distributions/Distributions.h"  // for InjectionDis...
#include "SIREN/interactions/InteractionCollection.h"

namespace siren {
namespace injection {

class Process {
private:
    siren::dataclasses::ParticleType primary_type;
    std::shared_ptr<interactions::InteractionCollection> interactions;
public:
    Process() = default;
    Process(siren::dataclasses::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    Process(Process const & other);
    Process(Process && other);
    Process & operator=(Process const & other);
    Process & operator=(Process && other);
    virtual ~Process() = default;

    void SetInteractions(std::shared_ptr<interactions::InteractionCollection> _interactions);
    std::shared_ptr<interactions::InteractionCollection> GetInteractions() const;
    void SetPrimaryType(siren::dataclasses::ParticleType _primary_type);
    siren::dataclasses::ParticleType GetPrimaryType() const;

    bool operator==(Process const & other) const;
    bool MatchesHead(std::shared_ptr<Process> const & other) const; // required to compared instances of derived classs
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryType", primary_type));
            archive(::cereal::make_nvp("Interactions", interactions));
        } else {
            throw std::runtime_error("Process only supports version <= 0!");
        }
    };
};

class PhysicalProcess : public Process {
protected:
    std::vector<std::shared_ptr<distributions::WeightableDistribution>> physical_distributions;
public:
    PhysicalProcess() = default;
    PhysicalProcess(siren::dataclasses::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    PhysicalProcess(PhysicalProcess const & other);
    PhysicalProcess(PhysicalProcess && other);
    PhysicalProcess & operator=(PhysicalProcess const & other);
    PhysicalProcess & operator=(PhysicalProcess && other);
    virtual ~PhysicalProcess() = default;
    virtual void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist);
    std::vector<std::shared_ptr<distributions::WeightableDistribution>> const & GetPhysicalDistributions() const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PhysicalDistributions", physical_distributions));
            archive(cereal::virtual_base_class<Process>(this));
        } else {
            throw std::runtime_error("PhysicalProcess only supports version <= 0!");
        }
    };
};

class PrimaryInjectionProcess : public PhysicalProcess {
protected:
    std::vector<std::shared_ptr<distributions::PrimaryInjectionDistribution>> primary_injection_distributions;
public:
    typedef distributions::PrimaryInjectionDistribution InjectionType;
    PrimaryInjectionProcess() = default;
    PrimaryInjectionProcess(siren::dataclasses::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    PrimaryInjectionProcess(PrimaryInjectionProcess const & other);
    PrimaryInjectionProcess(PrimaryInjectionProcess && other);
    PrimaryInjectionProcess & operator=(PrimaryInjectionProcess const & other);
    PrimaryInjectionProcess & operator=(PrimaryInjectionProcess && other);
    virtual ~PrimaryInjectionProcess() = default;
    virtual void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) override;
    virtual void AddPrimaryInjectionDistribution(std::shared_ptr<distributions::PrimaryInjectionDistribution> dist);
    std::vector<std::shared_ptr<distributions::PrimaryInjectionDistribution>> const & GetPrimaryInjectionDistributions() const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryInjectionDistributions", primary_injection_distributions));
            archive(cereal::virtual_base_class<PhysicalProcess>(this));
        } else {
            throw std::runtime_error("PrimaryInjectionProcess only supports version <= 0!");
        }
    };
};

class SecondaryInjectionProcess : public PhysicalProcess {
protected:
    std::vector<std::shared_ptr<distributions::SecondaryInjectionDistribution>> secondary_injection_distributions;
public:
    typedef distributions::SecondaryInjectionDistribution InjectionType;
    SecondaryInjectionProcess() = default;
    SecondaryInjectionProcess(siren::dataclasses::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    SecondaryInjectionProcess(SecondaryInjectionProcess const & other);
    SecondaryInjectionProcess(SecondaryInjectionProcess && other);
    SecondaryInjectionProcess & operator=(SecondaryInjectionProcess const & other);
    SecondaryInjectionProcess & operator=(SecondaryInjectionProcess && other);
    virtual ~SecondaryInjectionProcess() = default;
    virtual void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) override;
    virtual void AddSecondaryInjectionDistribution(std::shared_ptr<distributions::SecondaryInjectionDistribution> dist);
    std::vector<std::shared_ptr<distributions::SecondaryInjectionDistribution>> const & GetSecondaryInjectionDistributions() const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("SecondaryInjectionDistributions", secondary_injection_distributions));
            archive(cereal::virtual_base_class<PhysicalProcess>(this));
        } else {
            throw std::runtime_error("SecondaryInjectionProcess only supports version <= 0!");
        }
    };
};

} // namespace injection
} // namespace siren

CEREAL_CLASS_VERSION(siren::injection::Process, 0);

CEREAL_CLASS_VERSION(siren::injection::PhysicalProcess, 0);
CEREAL_REGISTER_TYPE(siren::injection::PhysicalProcess);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::injection::Process, siren::injection::PhysicalProcess);

CEREAL_CLASS_VERSION(siren::injection::SecondaryInjectionProcess, 0);
CEREAL_REGISTER_TYPE(siren::injection::SecondaryInjectionProcess);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::injection::PhysicalProcess, siren::injection::SecondaryInjectionProcess);

CEREAL_CLASS_VERSION(siren::injection::PrimaryInjectionProcess, 0);
CEREAL_REGISTER_TYPE(siren::injection::PrimaryInjectionProcess);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::injection::PhysicalProcess, siren::injection::PrimaryInjectionProcess);

#endif // SIREN_Process_H
