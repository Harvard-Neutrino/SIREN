#include "LeptonInjector/distributions/primary/type/PrimaryInjector.h"

#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/distributions/Distributions.h"

namespace LI {
namespace distributions {

//---------------
// class PrimaryInjector : InjectionDistribution
//---------------

PrimaryInjector::PrimaryInjector(LI::dataclasses::Particle::ParticleType primary_type, double primary_mass) :
    primary_type(primary_type),
    primary_mass(primary_mass)
{}

LI::dataclasses::Particle::ParticleType PrimaryInjector::PrimaryType() const {
    return primary_type;
}

double PrimaryInjector::PrimaryMass() const {
    return primary_mass;
}

void PrimaryInjector::Sample(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord & record) const {
    record.signature.primary_type = primary_type;
    record.primary_mass = primary_mass;
}
double PrimaryInjector::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    if(record.signature.primary_type != primary_type)
        return 0.0;
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

std::vector<std::string> PrimaryInjector::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryInjector::Name() const {
    return "PrimaryInjector";
}

std::shared_ptr<InjectionDistribution> PrimaryInjector::clone() const {
    return std::shared_ptr<InjectionDistribution>(new PrimaryInjector(*this));
}

bool PrimaryInjector::equal(WeightableDistribution const & other) const {
    const PrimaryInjector* x = dynamic_cast<const PrimaryInjector*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(primary_type, primary_mass)
            ==
            std::tie(x->primary_type, x->primary_mass);
}

bool PrimaryInjector::less(WeightableDistribution const & other) const {
    const PrimaryInjector* x = dynamic_cast<const PrimaryInjector*>(&other);
    return
        std::tie(primary_type, primary_mass)
        <
        std::tie(x->primary_type, x->primary_mass);
}

} // namespace distributions
} // namespace LeptonInjector
