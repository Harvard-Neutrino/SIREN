#pragma once
#ifndef LI_Process_H
#define LI_Process_H

#include <memory>                                        // for shared_ptr
#include <vector>                                        // for vector
#include <cstdint>                                       // for uint32_t
#include <stdexcept>                                     // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/Particle.h"         // for Particle
#include "LeptonInjector/distributions/Distributions.h"  // for InjectionDis...
#include "LeptonInjector/distributions/secondary/SecondaryInjectionDistribution.h"  // for InjectionDis...

namespace LI { namespace interactions { class InteractionCollection; } }

namespace LI {
namespace injection {

class Process {
private:
    LI::dataclasses::Particle::ParticleType primary_type;
    std::shared_ptr<interactions::InteractionCollection> interactions;
public:
    Process() = default;
    Process(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    Process(Process const & other);
    Process(Process && other);
    Process & operator=(Process const & other);
    Process & operator=(Process && other);
    virtual ~Process() = default;

    void SetInteractions(std::shared_ptr<interactions::InteractionCollection> _interactions);
    std::shared_ptr<interactions::InteractionCollection> GetInteractions() const;
    void SetPrimaryType(LI::dataclasses::Particle::ParticleType _primary_type);
    LI::dataclasses::Particle::ParticleType GetPrimaryType() const;

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
    PhysicalProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
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

class InjectionProcess : public PhysicalProcess {
protected:
    std::vector<std::shared_ptr<distributions::InjectionDistribution>> injection_distributions;
public:
    InjectionProcess() = default;
    InjectionProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    InjectionProcess(InjectionProcess const & other);
    InjectionProcess(InjectionProcess && other);
    InjectionProcess & operator=(InjectionProcess const & other);
    InjectionProcess & operator=(InjectionProcess && other);
    virtual ~InjectionProcess() = default;
    virtual void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) override;
    virtual void AddInjectionDistribution(std::shared_ptr<distributions::InjectionDistribution> dist);
    std::vector<std::shared_ptr<distributions::InjectionDistribution>> const & GetInjectionDistributions() const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InjectionDistributions", injection_distributions));
            archive(cereal::virtual_base_class<PhysicalProcess>(this));
        } else {
            throw std::runtime_error("InjectionProcess only supports version <= 0!");
        }
    };
};

class SecondaryInjectionProcess : public InjectionProcess {
protected:
    std::vector<std::shared_ptr<distributions::SecondaryInjectionDistribution>> secondary_injection_distributions;
public:
    SecondaryInjectionProcess() = default;
    SecondaryInjectionProcess(LI::dataclasses::Particle::ParticleType _primary_type, std::shared_ptr<interactions::InteractionCollection> _interactions);
    SecondaryInjectionProcess(SecondaryInjectionProcess const & other);
    SecondaryInjectionProcess(SecondaryInjectionProcess && other);
    SecondaryInjectionProcess & operator=(SecondaryInjectionProcess const & other);
    SecondaryInjectionProcess & operator=(SecondaryInjectionProcess && other);
    virtual ~SecondaryInjectionProcess() = default;
    virtual void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) override;
    virtual void AddInjectionDistribution(std::shared_ptr<distributions::InjectionDistribution> dist) override;
    virtual void AddSecondaryInjectionDistribution(std::shared_ptr<distributions::SecondaryInjectionDistribution> dist);
    std::vector<std::shared_ptr<distributions::SecondaryInjectionDistribution>> const & GetSecondaryInjectionDistributions() const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("SecondaryInjectionDistributions", injection_distributions));
            archive(cereal::virtual_base_class<InjectionProcess>(this));
        } else {
            throw std::runtime_error("SecondaryInjectionProcess only supports version <= 0!");
        }
    };
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::Process, 0);

#endif // LI_Process_H
