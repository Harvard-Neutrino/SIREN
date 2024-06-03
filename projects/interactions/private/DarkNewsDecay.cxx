#include "SIREN/interactions/DarkNewsDecay.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "SIREN/interactions/Decay.h"               // for Decay
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/utilities/Random.h"                  // for SIREN_random
#include "SIREN/utilities/Errors.h"                  // for PythonImplementationError


namespace siren {
namespace interactions {

DarkNewsDecay::DarkNewsDecay() {}

bool DarkNewsDecay::equal(Decay const & other) const {
    const DarkNewsDecay* x = dynamic_cast<const DarkNewsDecay*>(&other);

    if(!x)
        return false;
    else
        return true;
}

double DarkNewsDecay::TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const {
    return TotalDecayWidth(interaction.signature.primary_type);
}

double DarkNewsDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow TotalDecayWidth to call
    throw(siren::utilities::PythonImplementationError("DarkNewsDecay::TotalDecayWidth should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const {
     // Should be implemented on the python side
    // Not pure virtual in order to allow FinalStateProbability to call
    throw(siren::utilities::PythonImplementationError("DarkNewsDecay::TotalDecayWidthForFinalState should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow FinalStateProbability to call
    throw(siren::utilities::PythonImplementationError("DarkNewsDecay::DifferentialDecayWidth should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}

void DarkNewsDecay::SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord & interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(siren::utilities::PythonImplementationError("DarkNewsDecay::SampleRecordFromDarkNews should be implemented in Python!"));
}

void DarkNewsDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const {
    SampleRecordFromDarkNews(interaction, random);
}

} // namespace interactions
} // namespace siren
