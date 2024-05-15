#pragma once
#ifndef SIREN_Decay_H
#define SIREN_Decay_H

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
//namespace siren { namespace dataclasses { struct InteractionSignature; } }
//namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class Decay{
friend cereal::access;
private:
public:
    Decay();
    virtual ~Decay() {};
    bool operator==(Decay const & other) const;
    virtual bool equal(Decay const & other) const = 0;
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const = 0;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayLength(siren::dataclasses::InteractionRecord const & record) const;
    virtual double TotalDecayLengthForFinalState(siren::dataclasses::InteractionRecord const & record) const;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const = 0;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};

}; // class Decay

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

CEREAL_CLASS_VERSION(siren::interactions::Decay, 0);

#endif // SIREN_Decay_H
