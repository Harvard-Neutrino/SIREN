#include "LeptonInjector/crosssections/DarkNewsDecay.h"

#include <array>                                              // for array
#include <cmath>                                              // for sqrt, M_PI
#include <string>                                             // for basic_s...
#include <vector>                                             // for vector
#include <stddef.h>                                           // for size_t

#include <rk/geom3.hh>                                        // for Vector3
#include <rk/rk.hh>                                           // for P4, Boost

#include "LeptonInjector/crosssections/Decay.h"               // for Decay
#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/utilities/Random.h"                  // for LI_random
#include "LeptonInjector/utilities/Errors.h"                  // for PythonImplementationError


namespace LI {
namespace crosssections {

DarkNewsDecay::DarkNewsDecay() {}

pybind11::object DarkNewsDecay::get_self() {
    return pybind11::cast<pybind11::none>(Py_None);
}

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

double DarkNewsDecay::TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow TotalDecayWidth to call
    throw(LI::utilities::PythonImplementationError("DarkNewsDecay::TotalDecayWidth should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const {
     // Should be implemented on the python side
    // Not pure virtual in order to allow FinalStateProbability to call
    throw(LI::utilities::PythonImplementationError("DarkNewsDecay::TotalDecayWidthForFinalState should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow FinalStateProbability to call
    throw(LI::utilities::PythonImplementationError("DarkNewsDecay::DifferentialDecayWidth should be implemented in Python!"));
    return 0;
}

double DarkNewsDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
  double dd = DifferentialDecayWidth(record);
  double td = TotalDecayWidthForFinalState(record);
  if (dd == 0) return 0.;
  else if (td == 0) return 0.;
  else return dd/td;
}

dataclasses::InteractionRecord DarkNewsDecay::SampleRecordFromDarkNews(dataclasses::InteractionRecord & interaction,  std::shared_ptr<LI::utilities::LI_random> random) const {
    // Should be implemented on the python side
    // Not pure virtual in order to allow SampleFinalState to call
    throw(LI::utilities::PythonImplementationError("DarkNewsDecay::SampleRecordFromDarkNews should be implemented in Python!"));
    return dataclasses::InteractionRecord();
}

void DarkNewsDecay::SampleFinalState(dataclasses::InteractionRecord & interaction,  std::shared_ptr<LI::utilities::LI_random> random) const {
    interaction = SampleRecordFromDarkNews(interaction,random);
}


} // namespace crosssections
} // namespace LI