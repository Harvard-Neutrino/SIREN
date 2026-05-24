#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/utilities/Random.h"

#include <numeric>
#include <stdexcept>

namespace siren {
namespace injection {

int MultiChannelPhaseSpace::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    if (channels.empty()) {
        throw std::runtime_error("MultiChannelPhaseSpace has no channels");
    }

    // Pick a channel according to the weights
    double r = random->Uniform(0, 1);
    double cumulative = 0.0;
    int selected = static_cast<int>(channels.size()) - 1;
    for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
        cumulative += weights[i];
        if (r < cumulative) {
            selected = i;
            break;
        }
    }

    channels[selected]->Sample(random, detector_model, record);
    return selected;
}

double MultiChannelPhaseSpace::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    double density = 0.0;
    for (size_t i = 0; i < channels.size(); ++i) {
        density += weights[i] * channels[i]->Density(detector_model, record);
    }
    return density;
}

} // namespace injection
} // namespace siren
