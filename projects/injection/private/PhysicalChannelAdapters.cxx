#include "SIREN/injection/PhysicalChannelAdapters.h"

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/utilities/Random.h"

#include <stdexcept>

namespace siren {
namespace injection {

namespace {

// Infer topology from a Decay's signatures.
PhaseSpaceTopology InferDecayTopology(
    std::shared_ptr<siren::interactions::Decay> const & decay)
{
    auto sigs = decay->GetPossibleSignatures();
    if (sigs.empty()) return PhaseSpaceTopology::Unspecified;
    size_t n = sigs.front().secondary_types.size();
    for (auto const & sig : sigs) {
        if (sig.secondary_types.size() != n) return PhaseSpaceTopology::Unspecified;
    }
    if (n == 2) return PhaseSpaceTopology::Decay2Body;
    if (n == 3) return PhaseSpaceTopology::Decay3Body;
    if (n > 3)  return PhaseSpaceTopology::DecayNBody;
    return PhaseSpaceTopology::Unspecified;
}

// Infer topology from a CrossSection's signatures.
PhaseSpaceTopology InferCrossSectionTopology(
    std::shared_ptr<siren::interactions::CrossSection> const & xs)
{
    auto sigs = xs->GetPossibleSignatures();
    if (sigs.empty()) return PhaseSpaceTopology::Unspecified;
    size_t n = sigs.front().secondary_types.size();
    for (auto const & sig : sigs) {
        if (sig.secondary_types.size() != n) return PhaseSpaceTopology::Unspecified;
    }
    if (n == 2) return PhaseSpaceTopology::Scatter2to2;
    if (n == 3) return PhaseSpaceTopology::Scatter2to3;
    return PhaseSpaceTopology::Unspecified;
}

} // anonymous namespace

// ================================================================ //
//  PhysicalDecayChannel                                              //
// ================================================================ //

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    topology_ = InferDecayTopology(decay_);
    measure_ = siren::dataclasses::MeasureFromConvention(decay_->Convention());
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    siren::dataclasses::InteractionSignature const & signature)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    size_t n = signature.secondary_types.size();
    if (n == 2) topology_ = PhaseSpaceTopology::Decay2Body;
    else if (n == 3) topology_ = PhaseSpaceTopology::Decay3Body;
    else if (n > 3) topology_ = PhaseSpaceTopology::DecayNBody;
    else topology_ = PhaseSpaceTopology::Unspecified;

    measure_ = siren::dataclasses::MeasureFromConvention(decay_->Convention());
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    PhaseSpaceConvention convention)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    topology_ = InferDecayTopology(decay_);
    measure_ = siren::dataclasses::MeasureFromConvention(convention);
    // Override topology from convention if it carries topology info
    auto conv_topo = siren::dataclasses::TopologyFromConvention(convention, 0);
    if (conv_topo != PhaseSpaceTopology::Unspecified) {
        topology_ = conv_topo;
    }
}

void PhysicalDecayChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    siren::dataclasses::CrossSectionDistributionRecord csdr(record);
    decay_->SampleFinalState(csdr, random);
    csdr.Finalize(record);
}

double PhysicalDecayChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    return decay_->FinalStateProbability(record);
}

std::string PhysicalDecayChannel::Name() const {
    return "PhysicalDecay";
}

PhaseSpaceTopology PhysicalDecayChannel::Topology() const {
    return topology_;
}

PhaseSpaceMeasure PhysicalDecayChannel::Measure() const {
    return measure_;
}

PhaseSpaceConvention PhysicalDecayChannel::Convention() const {
    // Reconstruct legacy convention from topology + measure
    return PhaseSpaceChannel::Convention();
}

std::shared_ptr<siren::interactions::Decay>
PhysicalDecayChannel::GetDecay() const {
    return decay_;
}

// ================================================================ //
//  PhysicalCrossSectionChannel                                       //
// ================================================================ //

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section)
    : cross_section_(std::move(cross_section))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    topology_ = InferCrossSectionTopology(cross_section_);
    measure_ = siren::dataclasses::MeasureFromConvention(
        cross_section_->Convention());
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    siren::dataclasses::InteractionSignature const & signature)
    : cross_section_(std::move(cross_section))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    size_t n = signature.secondary_types.size();
    if (n == 2) topology_ = PhaseSpaceTopology::Scatter2to2;
    else if (n == 3) topology_ = PhaseSpaceTopology::Scatter2to3;
    else topology_ = PhaseSpaceTopology::Unspecified;

    measure_ = siren::dataclasses::MeasureFromConvention(
        cross_section_->Convention());
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    PhaseSpaceConvention convention)
    : cross_section_(std::move(cross_section))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    topology_ = InferCrossSectionTopology(cross_section_);
    measure_ = siren::dataclasses::MeasureFromConvention(convention);
}

void PhysicalCrossSectionChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    siren::dataclasses::CrossSectionDistributionRecord csdr(record);
    cross_section_->SampleFinalState(csdr, random);
    csdr.Finalize(record);
}

double PhysicalCrossSectionChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    return cross_section_->FinalStateProbability(record);
}

std::string PhysicalCrossSectionChannel::Name() const {
    return "PhysicalCrossSection";
}

PhaseSpaceTopology PhysicalCrossSectionChannel::Topology() const {
    return topology_;
}

PhaseSpaceMeasure PhysicalCrossSectionChannel::Measure() const {
    return measure_;
}

PhaseSpaceConvention PhysicalCrossSectionChannel::Convention() const {
    return PhaseSpaceChannel::Convention();
}

std::shared_ptr<siren::interactions::CrossSection>
PhysicalCrossSectionChannel::GetCrossSection() const {
    return cross_section_;
}

} // namespace injection
} // namespace siren
