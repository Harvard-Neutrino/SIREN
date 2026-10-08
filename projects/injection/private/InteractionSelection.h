#pragma once

#include <memory>
#include <vector>

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace detector { class DetectorModel; } }
namespace siren { namespace interactions {
class Interaction;
class InteractionCollection;
} }

namespace siren {
namespace injection {
namespace detail {

// One entry in the rate competition performed at an interaction vertex.
// The concrete interaction pointer is deliberately transient: records retain
// the observable signature, while generation keeps this pointer only long
// enough to dispatch final-state metadata and sampling to the model that won.
struct InteractionCandidate {
    siren::dataclasses::InteractionSignature signature;
    double target_mass;
    double rate;
    std::shared_ptr<siren::interactions::Interaction> interaction;
};

std::vector<InteractionCandidate> EnumerateInteractionCandidates(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record);

} // namespace detail
} // namespace injection
} // namespace siren
