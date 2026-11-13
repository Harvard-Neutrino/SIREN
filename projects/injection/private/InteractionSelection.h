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
