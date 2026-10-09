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

PhaseSpaceConvention ResolveCommonFinalStateConvention(
    PhaseSpaceConvention const & first,
    PhaseSpaceConvention const & second);

double ChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

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

double FixedVertexChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

} // namespace injection
} // namespace siren

#endif // SIREN_WeightingUtils_H
