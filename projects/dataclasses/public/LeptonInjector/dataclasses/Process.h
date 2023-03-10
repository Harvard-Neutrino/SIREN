#pragma once
#ifndef LI_Processes_H
#define LI_Processes_H

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

struct InjectionProcess
    LI::dataclasses::Particle::ParticleType primary_types;
    std::shared_ptr<crosssections::CrossSectionCollection> cross_sections;
    std::shared_ptr<distributions::PrimaryInjector> primary_injector;
    std::vector<std::shared_ptr<distributions::InjectionDistribution>> injeciton_distributions;
    // This funciton returns true if the given datum is the last entry to be saved in a tree
    std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> stopping_condition;
};

struct PhysicalProcess
    LI::dataclasses::Particle::ParticleType primary_types;
    std::shared_ptr<crosssections::CrossSectionCollection> cross_sections;
    std::vector<std::shared_ptr<distributions::InjectionDistribution>> physical_distributions;
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::Processes, 0);

#endif // LI_Processes_H
