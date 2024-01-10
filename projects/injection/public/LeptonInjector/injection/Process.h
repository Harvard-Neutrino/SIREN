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

namespace LI { namespace crosssections { class InteractionCollection; } }

namespace LI {
namespace injection {

struct Process {
    LI::dataclasses::Particle::ParticleType primary_type;
    std::shared_ptr<crosssections::InteractionCollection> cross_sections;
    void SetCrossSections(std::shared_ptr<crosssections::InteractionCollection> _cross_sections) {cross_sections = _cross_sections;}
    bool operator==(Process const & other) const;
    bool MatchesHead(std::shared_ptr<Process> const & other) const; // required to compared instances of derived structs
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryType", primary_type));
            archive(::cereal::make_nvp("CrossSections", cross_sections));
        } else {
            throw std::runtime_error("Process only supports version <= 0!");
        }
    };
};

struct InjectionProcess : Process {
    std::vector<std::shared_ptr<distributions::InjectionDistribution>> injection_distributions;
    void AddInjectionDistribution(std::shared_ptr<distributions::InjectionDistribution> dist) {
      for(auto _dist: injection_distributions) {
        if((*_dist) == (*dist)) return;
      }
      injection_distributions.push_back(dist);
    }
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InjectionDistributions", injection_distributions));
            archive(cereal::virtual_base_class<Process>(this));
        } else {
            throw std::runtime_error("InjectionProcess only supports version <= 0!");
        }
    };
};

struct PhysicalProcess : Process{
    std::vector<std::shared_ptr<distributions::WeightableDistribution>> physical_distributions;
    void AddPhysicalDistribution(std::shared_ptr<distributions::WeightableDistribution> dist) {
      for(auto _dist: physical_distributions) {
        if((*_dist) == (*dist)) return;
      }
      physical_distributions.push_back(dist);
    }
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

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::Process, 0);

#endif // LI_Process_H
