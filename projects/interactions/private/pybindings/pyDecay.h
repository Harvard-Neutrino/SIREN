#pragma once
#ifndef SIREN_pyDecay_H
#define SIREN_pyDecay_H

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline
#include "SIREN/utilities/Random.h" // for SIREN_random

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for Decay
class pyDecay : public Decay,Pybind11Trampoline<Decay, pyDecay> {
public:
    using Decay::Decay;
    pyDecay(Decay && parent) : Decay(std::move(parent)) {}
    pybind11::object self;

    double TotalDecayLength(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayLength,
            "TotalDecayLength",
            interaction
        )
    }

    double TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayLengthForFinalState,
            "TotalDecayLengthForFinalState",
            interaction
        )
    }
    
    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            interaction
        )
    }

    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            interaction
        )
    }

    double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            interaction
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            void,
            SampleFinalState,
            "SampleFinalState",
            record,
            random
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type
        )
    }

    std::vector<std::string> DensityVariables() const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<std::string>,
            DensityVariables,
            "DensityVariables"
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE_PURE(
            self,
            Decay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            record
        )
    }

}; // class pyDecay

} // namespace interactions
} // namespace siren

#endif // SIREN_pyDecay_H

