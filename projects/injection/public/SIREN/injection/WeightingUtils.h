#pragma once
#ifndef SIREN_WeightingUtils_H
#define SIREN_WeightingUtils_H

#include <memory>                 // for shared_ptr

namespace siren { namespace interactions { class InteractionCollection; } }
namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }

namespace siren {
namespace injection {

struct MultiChannelPhaseSpace;
class PhaseSpaceChannel;

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
// This is the original CrossSectionProbability function.
// Uses the Decay/CrossSection model's own FinalStateProbability.
double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

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

// Compute FinalStateProbability for the matched interaction only,
// without rate weighting. Used for Fixed vertices where there is no
// rate competition between channels.
double SelectedFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

} // namespace injection
} // namespace siren

#endif // SIREN_WeightingUtils_H
