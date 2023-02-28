#pragma once
#ifndef LI_InteractionRecord_H
#define LI_InteractionRecord_H

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

namespace LI {
namespace dataclasses {

struct InteractionRecord {
    dataclasses::InteractionSignature signature;
    double primary_mass = 0;
    std::array<double, 4> primary_momentum = {0, 0, 0, 0};
    double primary_helicity = 0;
    double target_mass = 0;
    std::array<double, 4> target_momentum = {0, 0, 0, 0};
    double target_helicity = 0;
    std::array<double, 3> interaction_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> interaction_parameters;
    bool operator==(InteractionRecord const & other) const;
    friend std::ostream& operator<<(std::ostream& os, InteractionRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InteractionSignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("TargetMass", target_mass));
            archive(::cereal::make_nvp("TargetMomentum", target_momentum));
            archive(::cereal::make_nvp("TargetHelicity", target_helicity));
            archive(::cereal::make_nvp("InteractionVertex", interaction_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("InteractionParameters", interaction_parameters));
        } else {
            throw std::runtime_error("InteractionRecord only supports version <= 0!");
        }
    };
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::InteractionRecord, 0);

#endif // LI_InteractionRecord_H
