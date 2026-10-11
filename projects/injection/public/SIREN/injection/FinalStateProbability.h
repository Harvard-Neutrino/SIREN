#pragma once
#ifndef SIREN_FinalStateProbability_H
#define SIREN_FinalStateProbability_H

#include <memory>

#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace interactions { class InteractionCollection; } }

namespace siren {
namespace injection {

class PhysicalProcess;
class MultiChannelPhaseSpace;

using PhaseSpaceConvention = siren::dataclasses::PhaseSpaceConvention;

// Select a common convention in which both densities can be compared.
PhaseSpaceConvention ResolveCommonFinalStateConvention(
    PhaseSpaceConvention const & first,
    PhaseSpaceConvention const & second);

// The convention of a process's final-state density at this record: its
// registered proposal's, or else the matching interaction model's.
PhaseSpaceConvention ProcessFinalStateConvention(
    PhysicalProcess const & process,
    siren::dataclasses::InteractionRecord const & record);

// The final-state factor a process contributes in the requested convention.
// Uses a rate-weighted probability for a propagated vertex or the channel
// selection probability times the selected density for a fixed vertex. A
// registered proposal replaces the interaction model's own density.
double ProcessFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    PhysicalProcess const & process,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention, bool for_generation = false);

// Final-state factor for the injector's sampled interaction, in the process's
// own convention. Includes the interaction selection used during generation.
double ProcessGenerationProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    PhysicalProcess const & process,
    siren::dataclasses::InteractionRecord const & record);

// Collection-level diagnostics for physical, unconditioned probabilities.
double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    MultiChannelPhaseSpace const & phase_space);

double ChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

double FixedVertexChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

} // namespace injection
} // namespace siren

#endif // SIREN_FinalStateProbability_H
