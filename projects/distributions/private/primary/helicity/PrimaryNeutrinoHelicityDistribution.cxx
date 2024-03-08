#include "SIREN/distributions/primary/helicity/PrimaryNeutrinoHelicityDistribution.h"

#include <array>                                              // for array
#include <string>                                             // for basic_s...
#include <stdlib.h>                                           // for abs

#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle
#include "SIREN/distributions/Distributions.h"       // for Injecti...
#include "SIREN/math/Vector3D.h"                     // for Vector3D

namespace SI {
namespace distributions {

//---------------
// class PrimaryNeutrinoHelicityDistribution : PrimaryInjectionDistribution
//---------------
void PrimaryNeutrinoHelicityDistribution::Sample(std::shared_ptr<SI::utilities::LI_random> rand, std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::PrimaryDistributionRecord & record) const {
    SI::dataclasses::ParticleType const & t = record.type;
    if(static_cast<int32_t>(t) > 0) // Particles are left handed, anti-particles are right handed
        record.SetHelicity(-0.5);
    else
        record.SetHelicity(0.5);
}

double PrimaryNeutrinoHelicityDistribution::GenerationProbability(std::shared_ptr<SI::detector::DetectorModel const> detector_model, std::shared_ptr<SI::interactions::InteractionCollection const> interactions, SI::dataclasses::InteractionRecord const & record) const {
    std::array<double, 4> const & mom = record.primary_momentum;
    SI::math::Vector3D dir(mom[1], mom[2], mom[3]);
    dir.normalize();

    if(abs(0.5 - abs(record.primary_helicity)) > 1e-9) // Helicity magnitude must be 0.5
        return 0.0;

    SI::dataclasses::ParticleType const & t = record.signature.primary_type;
    // Particles are left handed, anti-particles are right handed
    if(static_cast<int32_t>(t) > 0) {
        if(record.primary_helicity < 0) // expect opposite direction
            return 1.0;
        else
            return 0.0;
    } else {
        if(record.primary_helicity > 0) // expect same direction
            return 1.0;
        else
            return 0.0;
    }
}

PrimaryNeutrinoHelicityDistribution::PrimaryNeutrinoHelicityDistribution() {}

std::vector<std::string> PrimaryNeutrinoHelicityDistribution::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryNeutrinoHelicityDistribution::Name() const {
    return "PrimaryNeutrinoHelicityDistribution";
}

std::shared_ptr<PrimaryInjectionDistribution> PrimaryNeutrinoHelicityDistribution::clone() const {
    return std::shared_ptr<PrimaryNeutrinoHelicityDistribution>(new PrimaryNeutrinoHelicityDistribution(*this));
}

bool PrimaryNeutrinoHelicityDistribution::equal(WeightableDistribution const & other) const {
    const PrimaryNeutrinoHelicityDistribution* x = dynamic_cast<const PrimaryNeutrinoHelicityDistribution*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool PrimaryNeutrinoHelicityDistribution::less(WeightableDistribution const & other) const {
    return false;
}

} // namespace distributions
} // namespace SIREN
