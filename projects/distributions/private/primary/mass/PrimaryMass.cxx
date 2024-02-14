#include "LeptonInjector/distributions/primary/mass/PrimaryMass.h"

#include <tuple>                                              // for tie
#include <string>                                             // for basic_s...
#include <iostream>                                           // for operator<<
#include <stdlib.h>                                           // for abs

#include "LeptonInjector/dataclasses/InteractionRecord.h"     // for Interac...
#include "LeptonInjector/dataclasses/InteractionSignature.h"  // for Interac...
#include "LeptonInjector/dataclasses/Particle.h"              // for Particle
#include "LeptonInjector/distributions/Distributions.h"       // for Injecti...

namespace LI {
namespace distributions {

//---------------
// class PrimaryMass : PrimaryInjectionDistribution
//---------------

PrimaryMass::PrimaryMass(double primary_mass) :
    primary_mass(primary_mass)
{}

double PrimaryMass::GetPrimaryMass() const {
    return primary_mass;
}

void PrimaryMass::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::PrimaryDistributionRecord & record) const {
    record.SetMass(primary_mass);
}

double PrimaryMass::GenerationProbability(std::shared_ptr<LI::detector::DetectorModel const> detector_model, std::shared_ptr<LI::interactions::InteractionCollection const> interactions, LI::dataclasses::InteractionRecord const & record) const {
    if(2.0 * abs(record.primary_mass - primary_mass) / (record.primary_mass + primary_mass) > 1e-9) {
        std::cerr << "Event primary mass does not match injector primary mass!" << std::endl;
        std::cerr << "Event primary_mass: " << record.primary_mass << std::endl;
        std::cerr << "Injector primary_mass: " << primary_mass << std::endl;
        std::cerr << "Particle mass definitions should be consistent." << std::endl;
        std::cerr << "Are you using the wrong simulation?" << std::endl;
        return 0.0;
    }
    return 1.0;
}

std::vector<std::string> PrimaryMass::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryMass::Name() const {
    return "PrimaryMass";
}

std::shared_ptr<PrimaryInjectionDistribution> PrimaryMass::clone() const {
    return std::shared_ptr<PrimaryInjectionDistribution>(new PrimaryMass(*this));
}

bool PrimaryMass::equal(WeightableDistribution const & other) const {
    const PrimaryMass* x = dynamic_cast<const PrimaryMass*>(&other);

    if(!x)
        return false;
    else
        return primary_mass == x->primary_mass;
}

bool PrimaryMass::less(WeightableDistribution const & other) const {
    const PrimaryMass* x = dynamic_cast<const PrimaryMass*>(&other);
    return primary_mass == x->primary_mass;
}

} // namespace distributions
} // namespace LeptonInjector
