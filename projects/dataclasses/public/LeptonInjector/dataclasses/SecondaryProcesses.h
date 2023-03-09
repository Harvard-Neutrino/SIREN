#pragma once
#ifndef LI_SecondaryProcesses_H
#define LI_SecondaryProcesses_H

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
#include "LeptonInjector/crosssections/CrossSection.h"

namespace LI {
namespace dataclasses {

struct SecondaryProcesses {
    std::vector<LI::dataclasses::Particle::ParticleType> primary_types;
    std::vector<std::shared_ptr<crosssections::CrossSectionCollection>> processes;
    // This funciton returns true if the given datum is the last entry to be saved in a tree
    std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> stopping_condition;
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::SecondaryProcesses, 0);

#endif // LI_SecondaryProcesses_H
