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
#include "SIREN/dataclasses/PhaseSpaceConvention.h"
#include "SIREN/utilities/Random.h" // for SIREN_random
#include "SIREN/interactions/Interaction.h" // for Interaction

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class Decay : public Interaction {
friend cereal::access;
private:
public:
    Decay();
    virtual ~Decay() {};
    bool operator==(Decay const & other) const;
    virtual bool equal(Decay const & other) const = 0;
    // Optional whole-parent lifetime width. Zero keeps the legacy additive
    // contract; positive declarations in one collection must agree. Owned
    // partial widths remain additive through TotalDecayWidthAllFinalStates.
    virtual double ParentDecayWidth(dataclasses::InteractionRecord const &) const { return 0; }
    virtual double TotalDecayWidthAllFinalStates(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const = 0;
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual double TotalDecayLengthAllFinalStates(siren::dataclasses::InteractionRecord const & record) const;
    virtual double TotalDecayLength(siren::dataclasses::InteractionRecord const & record) const;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const = 0;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const = 0;
    virtual double SampleDecayTime(siren::dataclasses::CrossSectionDistributionRecord const & record, std::shared_ptr<siren::utilities::SIREN_random> random) const;
    virtual std::vector<double> SecondaryMasses(std::vector<siren::dataclasses::ParticleType> const & secondary_types) const;
    virtual std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const & record) const;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const = 0;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const = 0;
    virtual std::vector<std::string> DensityVariables() const = 0;
    // The topology and measure FinalStateProbability is differential in.
    // Undeclared measures are Unspecified, and an Unspecified density is never
    // converted or mixed with declared channels. MeasureForSignature defaults
    // to Measure(); override it only when channels use different measures.
    virtual siren::dataclasses::PhaseSpaceTopology Topology() const;
    virtual siren::dataclasses::PhaseSpaceMeasure Measure() const;
    virtual siren::dataclasses::PhaseSpaceTopology TopologyForSignature(
        siren::dataclasses::InteractionSignature const & signature) const;
    virtual siren::dataclasses::PhaseSpaceMeasure MeasureForSignature(
        siren::dataclasses::InteractionSignature const & signature) const;

    // An opaque final state: the model samples it but cannot evaluate its
    // density, so its FinalStateProbability for the signature is a placeholder
    // that is only correct where it cancels. The weighter then requires the
    // generation and physical processes to produce the signature with the same
    // final-state sampler and nothing else, and refuses a phase-space proposal
    // for it. FinalStateSamplerKey names the sampler: different model objects
    // that return the same non-empty key declare that they sample the same
    // final-state distribution (for example one generator wrapped with two
    // total cross sections). An empty key matches only the same object.
    virtual bool FinalStateDensityEvaluable(
        siren::dataclasses::InteractionSignature const & signature) const;
    virtual std::string FinalStateSamplerKey(
        siren::dataclasses::InteractionSignature const & signature) const;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};

}; // class Decay

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::Decay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::Decay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Interaction, siren::interactions::Decay);

#endif // SIREN_Decay_H
