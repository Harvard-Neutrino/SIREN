#include "SIREN/injection/PhysicalChannelAdapters.h"

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/utilities/Random.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace siren {
namespace injection {

namespace {

std::string Lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

bool ContainsAny(std::vector<std::string> const & values,
                 std::vector<std::string> const & needles) {
    for (auto const & value : values) {
        std::string lower = Lower(value);
        for (auto const & needle : needles) {
            if (lower.find(needle) != std::string::npos) return true;
        }
    }
    return false;
}

bool EqualsAny(std::vector<std::string> const & values,
               std::vector<std::string> const & needles) {
    for (auto const & value : values) {
        std::string lower = Lower(value);
        for (auto const & needle : needles) {
            if (lower == needle) return true;
        }
    }
    return false;
}

PhaseSpaceConvention InferDecayConvention(
    std::shared_ptr<siren::interactions::Decay> const & decay,
    size_t n_secondaries)
{
    if (n_secondaries == 2) {
        auto variables = decay->DensityVariables();
        if (ContainsAny(variables, {"energy", "e_", "lab", "cone", "biased"})) {
            return PhaseSpaceConvention::Custom;
        }
        return PhaseSpaceConvention::RestFrameSolidAngle;
    }

    if (n_secondaries == 3) {
        auto variables = decay->DensityVariables();
        if (ContainsAny(variables, {"energy", "e_", "lab"})) {
            return PhaseSpaceConvention::Custom;
        }
        if (ContainsAny(variables, {"dalitz", "s12", "s13", "s23", "s_12", "s_13", "s_23"})) {
            return PhaseSpaceConvention::Dalitz;
        }
        if (ContainsAny(variables, {"theta", "cos"})) {
            return PhaseSpaceConvention::HelicityAngles;
        }
        return PhaseSpaceConvention::Custom;
    }

    return PhaseSpaceConvention::Custom;
}

PhaseSpaceConvention InferDecayConvention(
    std::shared_ptr<siren::interactions::Decay> const & decay)
{
    auto signatures = decay->GetPossibleSignatures();
    if (signatures.empty()) return PhaseSpaceConvention::Custom;

    size_t n_secondaries = signatures.front().secondary_types.size();
    for (auto const & sig : signatures) {
        if (sig.secondary_types.size() != n_secondaries) {
            return PhaseSpaceConvention::Custom;
        }
    }

    return InferDecayConvention(decay, n_secondaries);
}

PhaseSpaceConvention InferDecayConvention(
    std::shared_ptr<siren::interactions::Decay> const & decay,
    siren::dataclasses::InteractionSignature const & signature)
{
    return InferDecayConvention(decay, signature.secondary_types.size());
}

PhaseSpaceConvention InferCrossSectionConvention(
    std::shared_ptr<siren::interactions::CrossSection> const & cross_section)
{
    auto variables = cross_section->DensityVariables();
    bool has_bjorken_x =
        ContainsAny(variables, {"bjorken x", "bjorken_x", "log10(x)", "log x"}) ||
        EqualsAny(variables, {"x"});
    bool has_bjorken_y =
        ContainsAny(variables, {"bjorken y", "bjorken_y", "log10(y)", "log y"}) ||
        EqualsAny(variables, {"y"});
    bool has_q2 =
        ContainsAny(variables, {"q^2", "mandelstam"}) ||
        EqualsAny(variables, {"q2", "q^2", "t", "-t"});

    if (has_bjorken_x && has_bjorken_y) {
        return PhaseSpaceConvention::BjorkenXY;
    }
    if (has_q2) {
        return PhaseSpaceConvention::MandelstamST;
    }
    if (has_bjorken_y) {
        return PhaseSpaceConvention::BjorkenXY;
    }
    return PhaseSpaceConvention::Custom;
}

PhaseSpaceConvention InferCrossSectionConvention(
    std::shared_ptr<siren::interactions::CrossSection> const & cross_section,
    siren::dataclasses::InteractionSignature const &)
{
    return InferCrossSectionConvention(cross_section);
}

} // anonymous namespace

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
    convention_ = InferDecayConvention(decay_);
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
    convention_ = InferDecayConvention(decay_, signature);
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
    convention_ = InferCrossSectionConvention(cross_section_);
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
    convention_ = InferCrossSectionConvention(cross_section_, signature);
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
