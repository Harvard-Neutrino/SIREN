#include "SIREN/interactions/ExternalGenerator.h"

#include <cmath>                                          // for isfinite
#include <set>                                            // for set
#include <sstream>                                        // for ostringstream
#include <algorithm>                                      // for find

#include "SIREN/dataclasses/InteractionRecord.h"          // for CrossSectionDistributionRecord
#include "SIREN/utilities/Random.h"                       // for SIREN_random

namespace siren {
namespace interactions {

namespace {

std::string Describe(siren::dataclasses::InteractionSignature const & s) {
    std::ostringstream oss;
    oss << s.primary_type << " + " << s.target_type << " ->";
    for(auto const & t : s.secondary_types) oss << ' ' << t;
    return oss.str();
}

// A zero-weight event still needs well-defined secondaries.
void ZeroSecondaries(siren::dataclasses::CrossSectionDistributionRecord & record) {
    for(auto & secondary : record.GetSecondaryParticleRecords()) {
        secondary.SetMass(0.0);
        secondary.SetFourMomentum({0.0, 0.0, 0.0, 0.0});
        secondary.SetHelicity(0.0);
    }
}

} // namespace

double ExternalGenerator::Threshold(siren::dataclasses::InteractionSignature const &) const {
    return 0.0;
}

double TrialGenerator::Threshold(siren::dataclasses::InteractionSignature const &) const {
    return 0.0;
}

/////////////////////////////////////////

BatchResampler::BatchResampler(std::shared_ptr<TrialGenerator> trials, unsigned int batch_size)
    : trials_(std::move(trials)), batch_size_(batch_size)
{
    if(!trials_)
        throw std::invalid_argument("BatchResampler: no trial generator");
    if(batch_size_ == 0)
        throw std::invalid_argument("BatchResampler: the batch size must be at least one");
}

std::vector<siren::dataclasses::InteractionSignature> BatchResampler::GetChannels() const {
    return trials_->GetChannels();
}

double BatchResampler::TotalCrossSection(siren::dataclasses::InteractionSignature const & channel,
                                         double primary_energy) const {
    return trials_->TotalCrossSection(channel, primary_energy);
}

double BatchResampler::Threshold(siren::dataclasses::InteractionSignature const & channel) const {
    return trials_->Threshold(channel);
}

void BatchResampler::GenerateFinalState(siren::dataclasses::CrossSectionDistributionRecord & record,
                                        std::shared_ptr<siren::utilities::SIREN_random> random) const {
    double const energy = record.GetPrimaryMomentum()[0];
    double const sigma = trials_->TotalCrossSection(record.GetSignature(), energy);
    double const envelope = trials_->Envelope(record.GetSignature(), energy);
    if(!(sigma > 0.0) || !std::isfinite(sigma) || !(envelope > 0.0) || !std::isfinite(envelope)) {
        throw std::runtime_error("BatchResampler: the channel " + Describe(record.GetSignature())
            + " was selected at E = " + std::to_string(energy) + " GeV but has total cross section "
            + std::to_string(sigma) + " and envelope " + std::to_string(envelope)
            + "; both must be positive and finite where the channel can be selected");
    }
    double sum = 0.0;
    std::int64_t emitted = 0;
    for(unsigned int j = 0; j < batch_size_; ++j) {
        double const w = trials_->Trial(record, random);
        if(!std::isfinite(w) || w < 0.0) {
            throw std::runtime_error("BatchResampler: trial weight " + std::to_string(w)
                + " for " + Describe(record.GetSignature()) + " is negative or non-finite");
        }
        if(w == 0.0) continue;
        ++emitted;
        sum += w;
        // Weighted reservoir of size one: trial j survives with probability w_j / sum.
        if(random->Uniform(0.0, 1.0) * sum < w)
            trials_->Keep();
    }
    record.interaction_int_parameters["batch_trials"] = batch_size_;
    record.interaction_int_parameters["batch_emitted"] = emitted;
    if(emitted == 0) {
        ZeroSecondaries(record);
        record.SetSamplerFactor(0.0);
        return;
    }
    trials_->Fill(record);
    record.SetSamplerFactor((sum / batch_size_) * envelope / sigma);
}

std::string BatchResampler::SamplerKey() const {
    return trials_->SamplerKey() + "|batch=" + std::to_string(batch_size_);
}

/////////////////////////////////////////

ExternalGeneratorCrossSection::ExternalGeneratorCrossSection(std::shared_ptr<ExternalGenerator> generator)
    : generator_(std::move(generator))
{
    if(!generator_)
        throw std::invalid_argument("ExternalGeneratorCrossSection: no generator");
    channels_ = generator_->GetChannels();
    if(channels_.empty())
        throw std::invalid_argument("ExternalGeneratorCrossSection: the generator declares no channels");
    if(generator_->SamplerKey().empty())
        throw std::invalid_argument("ExternalGeneratorCrossSection: the generator must name its sampler "
            "(SamplerKey) so that weights can check it");
}

bool ExternalGeneratorCrossSection::equal(CrossSection const & other) const {
    auto const * x = dynamic_cast<ExternalGeneratorCrossSection const *>(&other);
    return x && x->generator_ == generator_;
}

bool ExternalGeneratorCrossSection::HasChannel(siren::dataclasses::InteractionSignature const & signature) const {
    return std::find(channels_.begin(), channels_.end(), signature) != channels_.end();
}

double ExternalGeneratorCrossSection::TotalCrossSection(dataclasses::InteractionRecord const & record) const {
    if(!HasChannel(record.signature)) return 0.0;
    double const energy = record.primary_momentum[0];
    if(energy < generator_->Threshold(record.signature)) return 0.0;
    return generator_->TotalCrossSection(record.signature, energy);
}

double ExternalGeneratorCrossSection::DifferentialCrossSection(dataclasses::InteractionRecord const & record) const {
    // There is no differential cross section; FinalStateProbability is the
    // opaque placeholder that the weighter checks cancels.
    return TotalCrossSection(record);
}

double ExternalGeneratorCrossSection::InteractionThreshold(dataclasses::InteractionRecord const & record) const {
    return generator_->Threshold(record.signature);
}

void ExternalGeneratorCrossSection::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record,
                                                     std::shared_ptr<siren::utilities::SIREN_random> random) const {
    generator_->GenerateFinalState(record, random);
}

std::vector<siren::dataclasses::ParticleType> ExternalGeneratorCrossSection::GetPossibleTargets() const {
    std::set<siren::dataclasses::ParticleType> targets;
    for(auto const & c : channels_) targets.insert(c.target_type);
    return {targets.begin(), targets.end()};
}

std::vector<siren::dataclasses::ParticleType> ExternalGeneratorCrossSection::GetPossibleTargetsFromPrimary(
        siren::dataclasses::ParticleType primary_type) const {
    std::set<siren::dataclasses::ParticleType> targets;
    for(auto const & c : channels_) if(c.primary_type == primary_type) targets.insert(c.target_type);
    return {targets.begin(), targets.end()};
}

std::vector<siren::dataclasses::ParticleType> ExternalGeneratorCrossSection::GetPossiblePrimaries() const {
    std::set<siren::dataclasses::ParticleType> primaries;
    for(auto const & c : channels_) primaries.insert(c.primary_type);
    return {primaries.begin(), primaries.end()};
}

std::vector<dataclasses::InteractionSignature> ExternalGeneratorCrossSection::GetPossibleSignatures() const {
    return channels_;
}

std::vector<dataclasses::InteractionSignature> ExternalGeneratorCrossSection::GetPossibleSignaturesFromParents(
        siren::dataclasses::ParticleType primary_type,
        siren::dataclasses::ParticleType target_type) const {
    std::vector<dataclasses::InteractionSignature> out;
    for(auto const & c : channels_)
        if(c.primary_type == primary_type && c.target_type == target_type) out.push_back(c);
    return out;
}

double ExternalGeneratorCrossSection::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
    return TotalCrossSection(record) > 0.0 ? 1.0 : 0.0;
}

std::vector<std::string> ExternalGeneratorCrossSection::DensityVariables() const {
    return {};
}

bool ExternalGeneratorCrossSection::FinalStateDensityEvaluable(
        siren::dataclasses::InteractionSignature const &) const {
    return false;
}

std::string ExternalGeneratorCrossSection::FinalStateSamplerKey(
        siren::dataclasses::InteractionSignature const &) const {
    return generator_->SamplerKey();
}

} // namespace interactions
} // namespace siren
