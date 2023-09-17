#include "LeptonInjector/distributions/primary/helicity/PrimaryNeutrinoHelicityDistribution.h"

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace distributions {

//---------------
// class PrimaryNeutrinoHelicityDistribution : InjectionDistribution
//---------------
void PrimaryNeutrinoHelicityDistribution::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    LI::dataclasses::Particle::ParticleType & t = record.signature.primary_type;
    if(t > 0) // Particles are left handed, anti-particles are right handed
        record.primary_helicity = -0.5;
    else
        record.primary_helicity = 0.5;
}

double PrimaryNeutrinoHelicityDistribution::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    std::array<double, 4> const & mom = record.primary_momentum;
    LI::math::Vector3D dir(mom[1], mom[2], mom[3]);
    dir.normalize();

    if(abs(0.5 - abs(record.primary_helicity)) > 1e-9) // Helicity magnitude must be 0.5
        return 0.0;

    LI::dataclasses::Particle::ParticleType const & t = record.signature.primary_type;
    // Particles are left handed, anti-particles are right handed
    if(t > 0) {
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

std::shared_ptr<InjectionDistribution> PrimaryNeutrinoHelicityDistribution::clone() const {
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
} // namespace LeptonInjector
