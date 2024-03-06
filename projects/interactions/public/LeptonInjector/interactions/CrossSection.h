#pragma once
#ifndef LI_CrossSection_H
#define LI_CrossSection_H

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

#include "LeptonInjector/dataclasses/Particle.h"  // for Particle

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace LI { namespace dataclasses { struct InteractionSignature; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace interactions {

class CrossSection {
friend cereal::access;
private:
    void SampleFinalState(dataclasses::InteractionRecord &, std::shared_ptr<LI::utilities::LI_random>) const;
public:
    CrossSection();
    virtual ~CrossSection() {};
    bool operator==(CrossSection const & other) const;
    virtual bool equal(CrossSection const & other) const = 0;
    virtual double TotalCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalCrossSectionAllFinalStates(dataclasses::InteractionRecord const &) const;
    virtual double DifferentialCrossSection(dataclasses::InteractionRecord const &) const = 0;
    virtual double InteractionThreshold(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<LI::utilities::LI_random>) const = 0;

    virtual std::vector<LI::dataclasses::ParticleType> GetPossibleTargets() const = 0;
    virtual std::vector<LI::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::ParticleType primary_type) const = 0;
    virtual std::vector<LI::dataclasses::ParticleType> GetPossiblePrimaries() const = 0;
    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;

    virtual std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::ParticleType primary_type, LI::dataclasses::ParticleType target_type) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};
};

} // namespace interactions
} // namespace LI

CEREAL_CLASS_VERSION(LI::interactions::CrossSection, 0);

#endif // LI_CrossSection_H

