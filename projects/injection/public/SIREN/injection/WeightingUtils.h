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

// Compute the probability of the observed final state given an
// interaction at this vertex.  Uses the interaction collection's
// FinalStateProbability methods (physical or injection).
double CrossSectionProbability(std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &);

// Same as above, but replaces the FinalStateProbability with the
// multi-channel phase space density for the matched interaction.
// Used on the generation side when biased phase space sampling
// is active.
double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const>,
    std::shared_ptr<siren::interactions::InteractionCollection const>,
    siren::dataclasses::InteractionRecord const &,
    MultiChannelPhaseSpace const &);

} // namespace injection
} // namespace siren

#endif // SIREN_WeightingUtils_H
