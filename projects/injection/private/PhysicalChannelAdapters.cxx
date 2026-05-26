#include "SIREN/injection/PhysicalChannelAdapters.h"

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/utilities/Random.h"

#include <stdexcept>

namespace siren {
namespace injection {

// ================================================================ //
//  PhysicalDecayChannel                                              //
// ================================================================ //

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay)
    : decay_(std::move(decay))
    , convention_(PhaseSpaceConvention::Custom)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    convention_ = decay_->Convention();
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    siren::dataclasses::InteractionSignature const & signature)
    : decay_(std::move(decay))
    , convention_(PhaseSpaceConvention::Custom)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    convention_ = decay_->Convention();
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    PhaseSpaceConvention convention)
    : decay_(std::move(decay))
    , convention_(convention)
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
}

void PhysicalDecayChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    // Create the CrossSectionDistributionRecord wrapper that
    // Decay::SampleFinalState expects, delegate, then finalize
    // back into the InteractionRecord.
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

PhaseSpaceConvention PhysicalDecayChannel::Convention() const {
    return convention_;
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
    , convention_(PhaseSpaceConvention::Custom)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    convention_ = cross_section_->Convention();
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    siren::dataclasses::InteractionSignature const & signature)
    : cross_section_(std::move(cross_section))
    , convention_(PhaseSpaceConvention::Custom)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    convention_ = cross_section_->Convention();
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    PhaseSpaceConvention convention)
    : cross_section_(std::move(cross_section))
    , convention_(convention)
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
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

PhaseSpaceConvention PhysicalCrossSectionChannel::Convention() const {
    return convention_;
}

std::shared_ptr<siren::interactions::CrossSection>
PhysicalCrossSectionChannel::GetCrossSection() const {
    return cross_section_;
}

} // namespace injection
} // namespace siren
