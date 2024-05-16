#pragma once
#ifndef SIREN_CrossSection_H
#define SIREN_CrossSection_H

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

class CrossSection {
friend cereal::access;
private:
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const;
public:
    CrossSection();
    virtual ~CrossSection() {};
    bool operator==(CrossSection const & other) const;
    virtual bool equal(CrossSection const & other) const = 0;
    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalCrossSectionAllFinalStates(dataclasses::InteractionRecord const &) const;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const = 0;

    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const = 0;
    virtual std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;

    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};
}; // class CrossSection

// Trampoline class for CrossSection
class pyCrossSection : public CrossSection, public Pybind11Trampoline<CrossSection, pyCrossSection>  {
public:
    using CrossSection::CrossSection;
    pyCrossSection(CrossSection && parent) : CrossSection(std::move(parent)) {}
    pybind11::object self;

    bool equal(CrossSection const & other) const override { 
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            bool,
            equal,
            "equal",
            other
        )
    }

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            interaction
        )
    }

    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            CrossSection,
            double,
            TotalCrossSectionAllFinalStates,
            "TotalCrossSectionAllFinalStates",
            record
        )
    }

    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            interaction
        )
    }

    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            interaction
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            record,
            random
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            record
        )
    }

    std::vector<std::string> DensityVariables() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<std::string>,
            CrossSection,
            DensityVariables
        );
    }
}; // class pyCrossSection

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::CrossSection, 0);
using CrossSectionTrampolineType = Pybind11Trampoline<siren::interactions::CrossSection,siren::interactions::pyCrossSection>;
RegisterTrampolineCerealMethods(siren::interactions::CrossSection, siren::interactions::pyCrossSection, CrossSectionTrampolineType)

#endif // SIREN_CrossSection_H

