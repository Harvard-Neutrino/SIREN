#pragma once
#ifndef SIREN_pyDarkNewsCrossSection_H
#define SIREN_pyDarkNewsCrossSection_H

#include <memory>
#include <string>
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/interactions/DarkNewsCrossSection.h"  // for DarkNewsCrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// Trampoline class for DarkNewsCrossSection
class pyDarkNewsCrossSection : public DarkNewsCrossSection, public Pybind11Trampoline<DarkNewsCrossSection, pyDarkNewsCrossSection> {
public:
    using DarkNewsCrossSection::DarkNewsCrossSection;
    pyDarkNewsCrossSection(DarkNewsCrossSection && parent) : DarkNewsCrossSection(std::move(parent)) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
    }
    pyDarkNewsCrossSection(DarkNewsCrossSection const & parent) : DarkNewsCrossSection(parent) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
    }
    pybind11::object self;

    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            CrossSection,
            double,
            TotalCrossSectionAllFinalStates,
            "TotalCrossSectionAllFinalStates",
            std::cref(record)
        )
    }

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            std::cref(interaction)
        )
    }

    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            primary,
            energy,
            target
        )
    }

    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            std::cref(interaction)
        )
    }

    double DifferentialCrossSection(siren::dataclasses::ParticleType primary, siren::dataclasses::ParticleType target, double energy, double Q2) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            primary,
            target,
            energy,
            Q2
        )
    }

    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            std::cref(interaction)
        )
    }

    double Q2Min(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Min,
            "Q2Min",
            std::cref(interaction)
        )
    }

    double Q2Max(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Max,
            "Q2Max",
            std::cref(interaction)
        )
    }

    double TargetMass(dataclasses::ParticleType const & target_type) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TargetMass,
            "TargetMass",
            std::cref(target_type)
        )
    }

    std::vector<double> SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryMasses,
            "SecondaryMasses",
            std::cref(secondary_types)
        )
    }

    std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const & record) const override{
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryHelicities,
            "SecondaryHelicities",
            std::cref(record)
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
        )
    }

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        archive(cereal::base_class<Pybind11Trampoline<DarkNewsCrossSection, pyDarkNewsCrossSection>>(this));
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        archive(cereal::base_class<Pybind11Trampoline<DarkNewsCrossSection, pyDarkNewsCrossSection>>(this));
    }

};

} // namespace interactions
} // namespace siren

using DarkNewsCrossSectionTrampolineType = Pybind11Trampoline<siren::interactions::DarkNewsCrossSection,siren::interactions::pyDarkNewsCrossSection>;
RegisterTrampolineCerealMethods(siren::interactions::DarkNewsCrossSection, siren::interactions::pyDarkNewsCrossSection, DarkNewsCrossSectionTrampolineType)

#endif // SIREN_pyDarkNewsCrossSection_H

