#pragma once

namespace siren {
namespace injection {

class MultiChannelPhaseSpace;

// Transfer finalized coefficients between a pilot mixture and its production clone.
void CopyTunedWeights(MultiChannelPhaseSpace & destination,
                      MultiChannelPhaseSpace const & source);

} // namespace injection
} // namespace siren
