#pragma once
#ifndef SIREN_Interaction_H
#define SIREN_Interaction_H

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

namespace siren {
namespace interactions {

class Interaction {
public:
    virtual ~Interaction() = default;
};

};
};

CEREAL_CLASS_VERSION(siren::interactions::Interaction, 0);

#endif // SIREN_Interaction_H
