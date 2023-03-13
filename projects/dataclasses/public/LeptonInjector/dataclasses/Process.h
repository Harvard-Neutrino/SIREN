#pragma once
#ifndef LI_Process_H
#define LI_Process_H

#include <array>
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionTree.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/crosssections/CrossSection.h"

namespace LI {
namespace dataclasses {

struct Process {
    LI::dataclasses::Particle::ParticleType primary_type;
    std::shared_ptr<crosssections::CrossSectionCollection> cross_sections;
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

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::Process, 0);

#endif // LI_Process_H
