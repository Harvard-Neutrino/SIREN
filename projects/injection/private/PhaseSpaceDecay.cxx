#include "SIREN/injection/PhaseSpaceDecay.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Errors.h"
#include <cmath>
#include <sstream>

namespace siren { namespace injection {
namespace {
std::string ChannelState(std::shared_ptr<PhaseSpaceChannel> const & channel) {
    std::ostringstream stream;
    { cereal::BinaryOutputArchive archive(stream); archive(channel); }
    return stream.str();
}
}

PhaseSpaceDecay::PhaseSpaceDecay(siren::dataclasses::InteractionSignature signature,
    std::vector<double> masses, double partial_width, double total_width,
    std::shared_ptr<PhaseSpaceChannel> physical_channel)
    : signature_(std::move(signature)), masses_(std::move(masses)),
      partial_width_(partial_width), total_width_(total_width),
      physical_channel_(std::move(physical_channel)) { Validate(); }

void PhaseSpaceDecay::Validate() const {
    using PT = siren::dataclasses::ParticleType;
    if (!physical_channel_ || signature_.primary_type == PT::unknown
        || signature_.target_type != PT::Decay || masses_.size() < 2
        || masses_.size() != signature_.secondary_types.size()
        || !std::isfinite(partial_width_) || partial_width_ <= 0
        || !std::isfinite(total_width_) || total_width_ < partial_width_)
        throw std::invalid_argument("Invalid PhaseSpaceDecay signature, masses, channel or widths");
    for (double mass : masses_)
        if (!std::isfinite(mass) || mass < 0)
            throw std::invalid_argument("PhaseSpaceDecay masses must be finite and nonnegative");
    auto expected = masses_.size() == 2 ? PhaseSpaceTopology::Decay2Body
        : masses_.size() == 3 ? PhaseSpaceTopology::Decay3Body : PhaseSpaceTopology::DecayNBody;
    if (physical_channel_->Topology() != expected
        || physical_channel_->Measure().type == PhaseSpaceMeasure::Type::Unspecified)
        throw std::invalid_argument("PhaseSpaceDecay needs a compatible, explicit physical measure");
    auto measure = physical_channel_->Measure();
    if (measure.type == PhaseSpaceMeasure::Type::OnShellCascade) {
        if (masses_.size() != 3 || measure.spectator != 0 || measure.pair_first != 1
            || measure.pair_second != 2 || masses_[1] != masses_[2]
            || !(measure.pair_mass > 2 * masses_[1]))
            throw std::invalid_argument("PhaseSpaceDecay OnShellCascade requires spectator at 0 "
                "and exactly equal pair masses at (1, 2) below pair_mass/2");
    }
}

bool PhaseSpaceDecay::equal(siren::interactions::Decay const & other) const {
    auto p = dynamic_cast<PhaseSpaceDecay const *>(&other);
    return p && signature_ == p->signature_ && masses_ == p->masses_
        && partial_width_ == p->partial_width_ && total_width_ == p->total_width_
        && (physical_channel_ == p->physical_channel_
            || ChannelState(physical_channel_) == ChannelState(p->physical_channel_));
}
double PhaseSpaceDecay::ParentDecayWidth(siren::dataclasses::InteractionRecord const & r) const {
    return r.signature.primary_type == signature_.primary_type ? total_width_ : 0;
}
double PhaseSpaceDecay::TotalDecayWidthAllFinalStates(siren::dataclasses::InteractionRecord const & r) const {
    return TotalDecayWidth(r.signature.primary_type);
}
double PhaseSpaceDecay::TotalDecayWidth(siren::dataclasses::ParticleType p) const {
    return p == signature_.primary_type ? partial_width_ : 0;
}
double PhaseSpaceDecay::TotalDecayWidth(siren::dataclasses::InteractionRecord const & r) const {
    return r.signature == signature_ ? partial_width_ : 0;
}
double PhaseSpaceDecay::DifferentialDecayWidth(siren::dataclasses::InteractionRecord const & r) const {
    return TotalDecayWidth(r) * FinalStateProbability(r);
}
double PhaseSpaceDecay::FinalStateProbability(siren::dataclasses::InteractionRecord const & r) const {
    if (r.signature != signature_) return 0;
    return physical_channel_->Density(nullptr, r);
}
void PhaseSpaceDecay::SampleFinalState(siren::dataclasses::CrossSectionDistributionRecord & output,
    std::shared_ptr<siren::utilities::SIREN_random> random) const
{
    if (output.signature != signature_)
        throw std::invalid_argument("PhaseSpaceDecay cannot sample a different signature");
    auto record = output.record;
    // Include caller updates made on the mutable distribution record after its
    // construction, before the channel adds/updates its own sampling metadata.
    record.interaction_parameters = output.GetInteractionParameters();
    record.secondary_masses = masses_;
    record.secondary_momenta.resize(masses_.size());
    record.secondary_helicities.assign(masses_.size(), 0);
    physical_channel_->Sample(random, nullptr, record);
    for (std::size_t i = 0; i < masses_.size(); ++i) {
        auto & daughter = output.GetSecondaryParticleRecord(i);
        daughter.SetMass(masses_[i]);
        daughter.SetFourMomentum(record.secondary_momenta[i]);
        daughter.SetHelicity(record.secondary_helicities[i]);
    }
    output.SetInteractionParameters(record.interaction_parameters);
}
std::vector<double> PhaseSpaceDecay::SecondaryMasses(
    std::vector<siren::dataclasses::ParticleType> const & types) const {
    if (types != signature_.secondary_types)
        throw std::invalid_argument("PhaseSpaceDecay secondary types do not match its signature");
    return masses_;
}
std::vector<siren::dataclasses::InteractionSignature> PhaseSpaceDecay::GetPossibleSignatures() const {
    return {signature_};
}
std::vector<siren::dataclasses::InteractionSignature> PhaseSpaceDecay::GetPossibleSignaturesFromParent(
    siren::dataclasses::ParticleType p) const {
    return p == signature_.primary_type ? GetPossibleSignatures()
        : std::vector<siren::dataclasses::InteractionSignature>{};
}
std::vector<std::string> PhaseSpaceDecay::DensityVariables() const {
    return {siren::dataclasses::PhaseSpaceMeasureName(Measure())};
}
PhaseSpaceTopology PhaseSpaceDecay::Topology() const { return physical_channel_->Topology(); }
PhaseSpaceMeasure PhaseSpaceDecay::Measure() const { return physical_channel_->Measure(); }
}}
