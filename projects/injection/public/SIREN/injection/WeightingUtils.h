#pragma once
#ifndef SIREN_WeightingUtils_H
#define SIREN_WeightingUtils_H

#include <memory>                 // for shared_ptr

#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace injection {

struct MultiChannelPhaseSpace;
class PhaseSpaceChannel;

using PhaseSpaceConvention = siren::dataclasses::PhaseSpaceConvention;

// Return the natural convention (topology and measure) of the interaction
// model matching the record's signature. When multiple models match, elect a
// common convention that every model density can reach pointwise.
PhaseSpaceConvention SelectedFinalStateConvention(
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

// Elect one convention in which both densities can be evaluated pointwise.
// The result is directional: if one side has an explicit azimuth and the other
// declares a uniform integrated azimuth, the explicit convention wins because
// the marginal can be lifted while the joint density cannot be marginalized at
// one point. Throws when the topologies or measure families have no common
// pointwise convention.
PhaseSpaceConvention ResolveCommonFinalStateConvention(
    PhaseSpaceConvention const & first,
    PhaseSpaceConvention const & second);

// Compute the channel selection probability: the fraction of the
// total interaction rate that belongs to the observed signature.
//
//   P(channel) = rate_selected / rate_total
//
// where rate = density * cross_section (for scattering) or
// 1/decay_length (for decays).  This is independent of the
// final-state kinematics.
double ChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

// Compute the full final-state probability:
//   P(event) = P(channel) * FinalStateProbability
//
// This overload evaluates the density in the model's own convention; the overload below converts into a requested convention.
// Uses the Decay/CrossSection model's own FinalStateProbability, in the
// model's own convention.
double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

// CrossSectionProbability with the final-state density converted into the
// requested convention. The weighter uses this to evaluate both sides in the
// common convention elected for the weight ratio.
double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention);

// Compute the full final-state probability using a phase space
// channel density instead of FinalStateProbability:
//   P(event) = P(channel) * phase_space.Density(record)
//
// Used on the generation side when biased phase space sampling
// is active.
double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    MultiChannelPhaseSpace const & phase_space);

// CrossSectionProbabilityWithPhaseSpace with the mixture density evaluated in
// the requested convention (DensityIn), including the topology check.
double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    MultiChannelPhaseSpace const & phase_space,
    PhaseSpaceConvention const & convention);

// Compute FinalStateProbability conditional on the selected signature. A
// single matching model takes the exact direct-density path; when several
// concrete models share the signature, return their rate-weighted mixture.
// Used for Fixed vertices where the path-interaction factor is suppressed.
double SelectedFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

// SelectedFinalStateProbability with the density converted into the
// requested convention.
double SelectedFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention);

// Channel-selection probability for a Fixed vertex.
//
// The injector rate-selects the channel (Injector::SelectChannel) even at a
// Fixed vertex, so the generation/physical density must carry the same
// selected_rate/total_rate factor the sampler applied. This helper enumerates
// candidate signatures exactly like ChannelSelectionProbability, but:
//   - returns exactly 1.0 (an exact float no-op) when at most one candidate
//     signature competes for the record's parent (and target where
//     applicable), so single-channel Fixed configs are unaffected;
//   - otherwise returns selected_rate/total_rate.
// It fails loud (throws siren::utilities::WeightCalculationError) rather than
// silently returning a wrong factor if total_rate <= 0 or the result is
// non-finite.
double FixedVertexChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

} // namespace injection
} // namespace siren

#endif // SIREN_WeightingUtils_H
