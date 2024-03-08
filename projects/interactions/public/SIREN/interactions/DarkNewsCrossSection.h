#pragma once
#ifndef LI_DarkNewsCrossSection_H
#define LI_DarkNewsCrossSection_H

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
#include "SIREN/dataclasses/Particle.h"        // for Particlev

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace dataclasses { struct InteractionSignature; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace interactions {

class DarkNewsCrossSection : public CrossSection {
friend cereal::access;
private:
public:

    double m_ups;
    double m_target;
    double h_ups;
    double h_target;

    DarkNewsCrossSection();

    virtual pybind11::object get_representation();

    virtual bool equal(CrossSection const & other) const override;

    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double TotalCrossSection(SI::dataclasses::ParticleType primary, double energy, SI::dataclasses::ParticleType target) const;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialCrossSection(SI::dataclasses::ParticleType primary, SI::dataclasses::ParticleType target, double energy, double Q2) const; 
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const override;
    virtual double Q2Min(dataclasses::InteractionRecord const &) const;
    virtual double Q2Max(dataclasses::InteractionRecord const &) const;
    virtual double TargetMass(dataclasses::ParticleType const &) const;
    virtual std::vector<double> SecondaryMasses(std::vector<dataclasses::ParticleType> const &) const;
    virtual std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const &) const;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<SI::utilities::LI_random> random) const override;

    virtual std::vector<SI::dataclasses::ParticleType> GetPossibleTargets() const override = 0; // Requires Python-side implementation
    virtual std::vector<SI::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(SI::dataclasses::ParticleType primary_type) const override = 0; // Requires Python-side implementation
    virtual std::vector<SI::dataclasses::ParticleType> GetPossiblePrimaries() const override = 0; // Requires Python-side implementation
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override = 0; // Requires Python-side implementation
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(SI::dataclasses::ParticleType primary_type, SI::dataclasses::ParticleType target_type) const override = 0; // Requires Python-side implementation

    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;

public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        if(version == 0) {
            archive(cereal::virtual_base_class<CrossSection>(this));
        } else {
            throw std::runtime_error("DarkNewsCrossSection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace SI

CEREAL_CLASS_VERSION(SI::interactions::DarkNewsCrossSection, 0);
CEREAL_REGISTER_TYPE(SI::interactions::DarkNewsCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::interactions::CrossSection, SI::interactions::DarkNewsCrossSection);

#endif // LI_DarkNewsCrossSection_H
