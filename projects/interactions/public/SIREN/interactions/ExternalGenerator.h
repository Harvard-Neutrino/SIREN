#pragma once
#ifndef SIREN_ExternalGenerator_H
#define SIREN_ExternalGenerator_H

#include <memory>                                         // for shared_ptr
#include <string>                                         // for string
#include <vector>                                         // for vector
#include <cstdint>                                        // for uint32_t
#include <stdexcept>                                      // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/interactions/CrossSection.h"               // for CrossSection
#include "SIREN/dataclasses/InteractionSignature.h"        // for InteractionSignature
#include "SIREN/dataclasses/ParticleType.h"                // for ParticleType

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

// An external event generator, seen from SIREN. It produces a set of
// channels; each channel is the signature SIREN selects and weights (for
// example nu_mu + Ar40 -> mu- + Hadrons). For a selected channel the
// generator fills the vertex as it likes: the realised particles, their
// statuses and parents, integer and string attributes, and a sampler factor
// when it does not draw exactly from the physical density. Its final-state
// density is not evaluable, so the cross section that wraps it declares it
// opaque and names it with SamplerKey().
//
// Energies are lab primary energies in GeV and cross sections are in cm^2.
class ExternalGenerator {
public:
    virtual ~ExternalGenerator() = default;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetChannels() const = 0;
    virtual double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                                     double primary_energy) const = 0;
    virtual double Threshold(siren::dataclasses::InteractionSignature const & channel) const;
    virtual void GenerateFinalState(siren::dataclasses::CrossSectionDistributionRecord & record,
                                    std::shared_ptr<siren::utilities::SIREN_random> random) const = 0;
    virtual std::string SamplerKey() const = 0;
};

// A generator that works by single trials against an envelope, as an
// unweighting generator does: a trial either emits an event with weight
// W > 0 or is rejected (W = 0), and
//
//     E[W 1{emit}] * Envelope(channel, E) = TotalCrossSection(channel, E)
//
// at every energy. Keep() retains the last emitted trial and Fill() writes
// the retained one into a record (everything but the sampler factor).
class TrialGenerator {
public:
    virtual ~TrialGenerator() = default;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetChannels() const = 0;
    virtual double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                                     double primary_energy) const = 0;
    virtual double Threshold(siren::dataclasses::InteractionSignature const & channel) const;
    virtual double Envelope(siren::dataclasses::InteractionSignature const & channel,
                            double primary_energy) const = 0;
    virtual double Trial(siren::dataclasses::CrossSectionDistributionRecord const & record,
                         std::shared_ptr<siren::utilities::SIREN_random> random) const = 0;
    virtual void Keep() const = 0;
    virtual void Fill(siren::dataclasses::CrossSectionDistributionRecord & record) const = 0;
    virtual std::string SamplerKey() const = 0;
};

// Turns a TrialGenerator into an ExternalGenerator by fixed-batch
// resampling. Each vertex runs M trials, keeps one emitted trial chosen with
// probability W_j / sum W, and gets the sampler factor
//
//     (sum_j W_j / M) * Envelope / TotalCrossSection.
//
// The factor has expectation one at each energy and the weighted final state
// is the generator's physical distribution, for any M. A batch with no
// emitted trial is a valid zero-weight event: its particles are set to zero
// momentum and its factor to zero. The vertex records the batch size and the
// number of emitted trials as the integer attributes batch_trials and
// batch_emitted.
class BatchResampler : public ExternalGenerator {
public:
    BatchResampler(std::shared_ptr<TrialGenerator> trials, unsigned int batch_size);
    std::vector<siren::dataclasses::InteractionSignature> GetChannels() const override;
    double TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                             double primary_energy) const override;
    double Threshold(siren::dataclasses::InteractionSignature const & channel) const override;
    void GenerateFinalState(siren::dataclasses::CrossSectionDistributionRecord & record,
                            std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::string SamplerKey() const override;
    std::shared_ptr<TrialGenerator> GetTrialGenerator() const { return trials_; }
    unsigned int GetBatchSize() const { return batch_size_; }
private:
    std::shared_ptr<TrialGenerator> trials_;
    unsigned int batch_size_;
};

// A CrossSection that hands its channels' final states to an external
// generator. Total cross sections come from the generator; the final-state
// density is opaque (FinalStateProbability is a constant that the weighter
// requires to cancel) and named by the generator's SamplerKey().
//
// The generator is external state, so this cross section cannot be archived.
class ExternalGeneratorCrossSection : public CrossSection {
friend cereal::access;
public:
    explicit ExternalGeneratorCrossSection(std::shared_ptr<ExternalGenerator> generator);
    std::shared_ptr<ExternalGenerator> GetGenerator() const { return generator_; }

    bool equal(CrossSection const & other) const override;
    double TotalCrossSection(dataclasses::InteractionRecord const & record) const override;
    double DifferentialCrossSection(dataclasses::InteractionRecord const & record) const override;
    double InteractionThreshold(dataclasses::InteractionRecord const & record) const override;
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record,
                          std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override;
    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(
        siren::dataclasses::ParticleType primary_type) const override;
    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    std::vector<dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType primary_type,
        siren::dataclasses::ParticleType target_type) const override;
    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
    std::vector<std::string> DensityVariables() const override;
    bool FinalStateDensityEvaluable(
        siren::dataclasses::InteractionSignature const & signature) const override;
    std::string FinalStateSamplerKey(
        siren::dataclasses::InteractionSignature const & signature) const override;

    template<typename Archive>
    void save(Archive &, std::uint32_t const) const {
        throw std::runtime_error("ExternalGeneratorCrossSection cannot be archived: "
            "its generator is external state; rebuild it from the generator's own configuration");
    }
    template<typename Archive>
    void load(Archive &, std::uint32_t const) {
        throw std::runtime_error("ExternalGeneratorCrossSection cannot be loaded from an archive");
    }
private:
    ExternalGeneratorCrossSection() = default;  // for cereal registration only
    bool HasChannel(siren::dataclasses::InteractionSignature const & signature) const;
    std::shared_ptr<ExternalGenerator> generator_;
    std::vector<siren::dataclasses::InteractionSignature> channels_;
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::ExternalGeneratorCrossSection, 0);
CEREAL_REGISTER_TYPE(siren::interactions::ExternalGeneratorCrossSection);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::CrossSection, siren::interactions::ExternalGeneratorCrossSection);

#endif // SIREN_ExternalGenerator_H
