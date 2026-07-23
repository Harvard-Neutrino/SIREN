#pragma once
#ifndef SIREN_TrivialCrossSection_H
#define SIREN_TrivialCrossSection_H

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/interactions/CrossSection.h"  // for CrossSection
#include "SIREN/dataclasses/Particle.h"        // for Particle

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// A configurable pass-through cross section. The final state is the primary
// and the target unchanged, with unit final-state probability, so the physics
// content is the total cross section alone. Intended for two roles: a
// generation stand-in whose scale cancels in the event weight, and a physical
// cross section whose total is tabulated from an external generator and
// applied at weighting time.
class TrivialCrossSection : public CrossSection {
friend cereal::access;
private:
    // sigma(E) per primary type: ascending energies in GeV, values in cm^2.
    // Linear interpolation between knots; zero below the first knot; clamped
    // to the last value above the last knot.
    std::map<siren::dataclasses::ParticleType, std::vector<double>> energies_;
    std::map<siren::dataclasses::ParticleType, std::vector<double>> cross_sections_;
    std::vector<siren::dataclasses::ParticleType> targets_;

    TrivialCrossSection() = default;
public:
    // Constant cross section for every listed primary.
    TrivialCrossSection(double cross_section_cm2,
                        std::vector<siren::dataclasses::ParticleType> primary_types,
                        std::vector<siren::dataclasses::ParticleType> target_types);
    // Tabulated cross section per primary: {primary: (energies_GeV, sigma_cm2)}.
    TrivialCrossSection(std::map<siren::dataclasses::ParticleType,
                                 std::pair<std::vector<double>, std::vector<double>>> tables,
                        std::vector<siren::dataclasses::ParticleType> target_types);

    virtual bool equal(CrossSection const & other) const override;

    double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy) const;
    double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random> random) const override;

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override;

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

    virtual std::vector<std::string> DensityVariables() const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Energies", energies_));
            archive(::cereal::make_nvp("CrossSections", cross_sections_));
            archive(::cereal::make_nvp("Targets", targets_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("TrivialCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Energies", energies_));
            archive(::cereal::make_nvp("CrossSections", cross_sections_));
            archive(::cereal::make_nvp("Targets", targets_));
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("TrivialCrossSection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::TrivialCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::TrivialCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::TrivialCrossSection);

#endif // SIREN_TrivialCrossSection_H
