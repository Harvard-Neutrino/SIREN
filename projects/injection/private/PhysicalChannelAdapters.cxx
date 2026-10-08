#include "SIREN/injection/PhysicalChannelAdapters.h"

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/PhaseSpaceConvention.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include <stdexcept>

namespace siren {
namespace injection {

namespace {

void ValidateConventionOverride(
    PhaseSpaceTopology signature_topology,
    PhaseSpaceConvention const & convention,
    char const * adapter_name)
{
    if (convention.topology != signature_topology) {
        throw siren::utilities::ConfigurationError(
            std::string(adapter_name) +
            " convention topology does not match the selected signature "
            "[siren-docs: errors#measure-compat]");
    }
    if (convention.measure.type == PhaseSpaceMeasure::Type::Unspecified) {
        throw siren::utilities::ConfigurationError(
            std::string(adapter_name) +
            " convention override requires an explicit measure "
            "[siren-docs: errors#measure-compat]");
    }
}

} // anonymous namespace

// ================================================================ //
//  PhysicalDecayChannel                                              //
// ================================================================ //

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    topology_ = decay_->Topology();
    measure_ = decay_->Measure();
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    siren::dataclasses::InteractionSignature const & signature)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    topology_ = decay_->TopologyForSignature(signature);
    measure_ = decay_->MeasureForSignature(signature);
}

PhysicalDecayChannel::PhysicalDecayChannel(
    std::shared_ptr<siren::interactions::Decay> decay,
    siren::dataclasses::InteractionSignature const & signature,
    PhaseSpaceConvention const & convention)
    : decay_(std::move(decay))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!decay_) {
        throw std::runtime_error("PhysicalDecayChannel requires a non-null Decay");
    }
    PhaseSpaceTopology signature_topology =
        decay_->TopologyForSignature(signature);
    ValidateConventionOverride(
        signature_topology, convention, "PhysicalDecayChannel");
    topology_ = convention.topology;
    measure_ = convention.measure;
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
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    topology_ = cross_section_->Topology();
    measure_ = cross_section_->Measure();
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    siren::dataclasses::InteractionSignature const & signature)
    : cross_section_(std::move(cross_section))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    topology_ = cross_section_->TopologyForSignature(signature);
    measure_ = cross_section_->MeasureForSignature(signature);
}

PhysicalCrossSectionChannel::PhysicalCrossSectionChannel(
    std::shared_ptr<siren::interactions::CrossSection> cross_section,
    siren::dataclasses::InteractionSignature const & signature,
    PhaseSpaceConvention const & convention)
    : cross_section_(std::move(cross_section))
    , topology_(PhaseSpaceTopology::Unspecified)
    , measure_(PhaseSpaceMeasure::Unspecified())
{
    if (!cross_section_) {
        throw std::runtime_error(
            "PhysicalCrossSectionChannel requires a non-null CrossSection");
    }
    PhaseSpaceTopology signature_topology =
        cross_section_->TopologyForSignature(signature);
    ValidateConventionOverride(
        signature_topology, convention, "PhysicalCrossSectionChannel");
    topology_ = convention.topology;
    measure_ = convention.measure;
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

std::shared_ptr<siren::interactions::CrossSection>
PhysicalCrossSectionChannel::GetCrossSection() const {
    return cross_section_;
}

} // namespace injection
} // namespace siren
