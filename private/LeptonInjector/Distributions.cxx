#include <tuple>

#include "LeptonInjector/Distributions.h"
#include "earthmodel-service/EarthModelCalculator.h"

namespace LeptonInjector {

namespace {
    double one_minus_exp_of_negative(double x) {
        if(x < 1e-1) {
            return std::exp(std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0);
        } else {
            return 1.0 - std::exp(-x);
        }
    }
}

//---------------
// class WeightableDistribution
//---------------

std::vector<std::string> WeightableDistribution::DensityVariables() const {
    return {};
}

bool WeightableDistribution::operator==(WeightableDistribution const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool WeightableDistribution::operator<(WeightableDistribution const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class InjectionDistribution
//---------------

void InjectionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
}


//---------------
// class PrimaryInjector : InjectionDistribution
//---------------

PrimaryInjector::PrimaryInjector(LeptonInjector::Particle::ParticleType primary_type, double primary_mass) :
    primary_type(primary_type),
    primary_mass(primary_mass)
{}

LeptonInjector::Particle::ParticleType PrimaryInjector::PrimaryType() const {
    return primary_type;
}

double PrimaryInjector::PrimaryMass() const {
    return primary_mass;
}

void PrimaryInjector::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    record.signature.primary_type = primary_type;
    record.primary_mass = primary_mass;
}
double PrimaryInjector::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
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


//---------------
// class TargetMomentumDistribution : InjectionDistribution
//---------------

void TargetMomentumDistribution::Sample(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        CrossSectionCollection const & cross_sections,
        InteractionRecord & record) const {
    record.target_momentum = SampleMomentum(rand, earth_model, cross_sections, record);
}

std::vector<std::string> TargetMomentumDistribution::DensityVariables() const {
    return std::vector<std::string>{"TargetMomentum"};
}

//---------------
// class TargetAtRest : TargetMomentumDistribution : InjectionDistribution
//---------------
std::array<double, 4> TargetAtRest::SampleMomentum(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        CrossSectionCollection const & cross_sections,
        InteractionRecord const & record) const {
    return std::array<double, 4>{record.target_mass, 0, 0, 0};
}

double TargetAtRest::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    return 1.0;
}

std::vector<std::string> TargetAtRest::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> TargetAtRest::clone() const {
    return std::shared_ptr<InjectionDistribution>(new TargetAtRest(*this));
}

std::string TargetAtRest::Name() const {
    return "TargetAtRest";
}

bool TargetAtRest::equal(WeightableDistribution const & other) const {
    const TargetAtRest* x = dynamic_cast<const TargetAtRest*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool TargetAtRest::less(WeightableDistribution const & other) const {
    return false;
}

//---------------
// class PrimaryEnergyDistribution : InjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        CrossSectionCollection const & cross_sections,
        InteractionRecord & record) const {
    record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

//---------------
// class PowerLaw : PrimaryEnergyDistribution
//---------------
PowerLaw::PowerLaw(double powerLawIndex, double energyMin, double energyMax)
    : powerLawIndex(powerLawIndex)
    , energyMin(energyMin)
    , energyMax(energyMax)
{}

double PowerLaw::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    if(energyMin == energyMax)
        return energyMin; //return the only allowed energy

    if(powerLawIndex == 1.0) //sample uniformly in log space
        return pow(10.0, rand->Uniform(log10(energyMin), log10(energyMax)));
    else {
        double u = rand->Uniform();
        double energyP = (1 - u) * pow(energyMin, 1 - powerLawIndex) + u * pow(energyMax, 1 - powerLawIndex);
        return pow(energyP, 1 / (1 - powerLawIndex));
    }
}

double PowerLaw::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    if(energyMin == energyMax)
        return 1.0; // only one allowed energy

    if(powerLawIndex == 1.0)
        return 1.0 / (record.primary_momentum[0] * log(energyMax / energyMin));
    else {
        return pow(record.primary_momentum[0], -powerLawIndex) * (powerLawIndex - 1.0) * (pow(energyMin, powerLawIndex - 1.0) - pow(energyMax, powerLawIndex - 1.0));
    }
}

std::string PowerLaw::Name() const {
    return "PowerLaw";
}

std::shared_ptr<InjectionDistribution> PowerLaw::clone() const {
    return std::shared_ptr<InjectionDistribution>(new PowerLaw(*this));
}

bool PowerLaw::equal(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, powerLawIndex)
            ==
            std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

bool PowerLaw::less(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);
    return
        std::tie(energyMin, energyMax, powerLawIndex)
        <
        std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

//---------------
// class ModifiedMoyalPlusExponentialEnergyDistribution : PrimaryEnergyDistribution
//---------------

double ModifiedMoyalPlusExponentialEnergyDistribution::unnormed_pdf(double energy) const {
    double x = (energy - mu) / sigma;
    double moyal = (A / sigma) * std::exp(-(x + std::exp(-x)/2)) / std::sqrt(2.0 * M_PI);
    double exponential = (B / l) * std::exp(-energy / l);
    return moyal + exponential;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

ModifiedMoyalPlusExponentialEnergyDistribution::ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , mu(mu)
    , sigma(sigma)
    , A(A)
    , l(l)
    , B(B)
{
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = earthmodel::Integration::rombergIntegrate(integrand, energyMin, energyMax);
}

double ModifiedMoyalPlusExponentialEnergyDistribution::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    // Metropolis-Hastings algorithm to sample from PDF.
    // Pass in a function pointer for the PDF

    double energy, density, test_energy, test_density, odds;
    bool accept;

    // sample an initial point uniformly
    energy = rand->Uniform(energyMin, energyMax);
    density = pdf(energy);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        test_energy = rand->Uniform(energyMin, energyMax);
        test_density = pdf(test_energy);
        odds = test_density / density;
        accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
        if(accept) {
            energy = test_energy;
            density = test_density;
        }
    }

    return energy;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string ModifiedMoyalPlusExponentialEnergyDistribution::Name() const {
    return "ModifiedMoyalPlusExponentialEnergyDistribution";
}

std::shared_ptr<InjectionDistribution> ModifiedMoyalPlusExponentialEnergyDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ModifiedMoyalPlusExponentialEnergyDistribution(*this));
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::equal(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, mu, sigma, A, l, B)
            ==
            std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::less(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, mu, sigma, A, l, B)
        <
        std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

//---------------
// class PrimaryDirectionDistribution : InjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir = SampleDirection(rand, earth_model, cross_sections, record);
    double energy = record.primary_momentum[0];
    double mass = record.primary_mass;
    double momentum = std::sqrt(energy*energy - mass*mass);
    record.primary_momentum[1] = momentum * dir.GetX();
    record.primary_momentum[2] = momentum * dir.GetY();
    record.primary_momentum[3] = momentum * dir.GetZ();
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

//---------------
// class IsotropicDirection : PrimaryDirectionDistribution
//---------------
earthmodel::Vector3D IsotropicDirection::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    double nx = rand->Uniform(0, 1);
    double ny = rand->Uniform(0, 1);
    double nz = rand->Uniform(0, 1);
    earthmodel::Vector3D res(nx, ny, nz);
    res.normalize();
    return res;
}

double IsotropicDirection::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    return 1.0 / (4.0 * M_PI);
}

std::shared_ptr<InjectionDistribution> IsotropicDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new IsotropicDirection(*this));
}

std::string IsotropicDirection::Name() const {
    return "IsotropicDirection";
}

bool IsotropicDirection::equal(WeightableDistribution const & other) const {
    const IsotropicDirection* x = dynamic_cast<const IsotropicDirection*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool IsotropicDirection::less(WeightableDistribution const & other) const {
    return false;
}

//---------------
// class FixedDirection : PrimaryDirectionDistribution
//---------------
earthmodel::Vector3D FixedDirection::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    return dir;
}

double FixedDirection::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    if(abs(1.0 - earthmodel::scalar_product(dir, event_dir)) < 1e-9)
        return 1.0;
    else
        return 0.0;
}

std::vector<std::string> FixedDirection::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> FixedDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new FixedDirection(*this));
}

std::string FixedDirection::Name() const {
    return "FixedDirection";
}

bool FixedDirection::equal(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9);
}

bool FixedDirection::less(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);
    if(abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        double X = dir.GetX();
        double Y = dir.GetY();
        double Z = dir.GetZ();
        double other_X = dir.GetX();
        double other_Y = dir.GetY();
        double other_Z = dir.GetZ();
        return
            std::tie(X, Y, Z)
            <
            std::tie(other_X, other_Y, other_Z);
    }
}

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(earthmodel::Vector3D dir, double opening_angle) : dir(dir), opening_angle(opening_angle) {
    this->dir.normalize();
    if(this->dir == earthmodel::Vector3D(0,0,1)) {
        rotation = earthmodel::Quaternion(0,0,0,1);
    } else {
        earthmodel::Vector3D r = cross_product(earthmodel::Vector3D(0, 0, 1), dir);
        r.normalize();
        rotation = earthmodel::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
    }
}

earthmodel::Vector3D Cone::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const{
    double theta = cos(rand->Uniform(acos(opening_angle), 1));
    double phi = rand->Uniform(0, 2.0 * M_PI);
    earthmodel::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(earthmodel::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double theta = acos(earthmodel::scalar_product(dir, event_dir));
    if(theta < opening_angle)
        return 1.0 / (2.0 * M_PI * (1.0 - cos(opening_angle)));
    else
        return 0.0;
}

std::shared_ptr<InjectionDistribution> Cone::clone() const {
    return std::shared_ptr<InjectionDistribution>(new Cone(*this));
}

std::string Cone::Name() const {
    return "Cone";
}

bool Cone::equal(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9
            and opening_angle == x->opening_angle);
}

bool Cone::less(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);
    if(abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        return opening_angle < x->opening_angle;
    }
}

//---------------
// class VertexPositionDistribution : InjectionDistribution
//---------------
void VertexPositionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
    record.interaction_vertex[0] = pos.GetX();
    record.interaction_vertex[1] = pos.GetY();
    record.interaction_vertex[2] = pos.GetZ();
}

std::vector<std::string> VertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"InteractionVertexPosition"};
}

//---------------
// class CylinderVolumePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D CylinderVolumePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), z);
    return cylinder.LocalToGlobalPosition(pos);
}

double CylinderVolumePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D pos(record.interaction_vertex);
    double z = pos.GetZ();
    double r = sqrt(pos.GetX() * pos.GetX() + pos.GetY() * pos.GetY());
    if(abs(z) >= 0.5 * cylinder.GetZ()
            or r <= cylinder.GetInnerRadius()
            or r >= cylinder.GetRadius()) {
        return 0.0;
    } else {
        return 1.0 / ((cylinder.GetRadius() * cylinder.GetRadius() - cylinder.GetInnerRadius() * cylinder.GetInnerRadius()) * cylinder.GetZ());
    }
}


CylinderVolumePositionDistribution::CylinderVolumePositionDistribution(earthmodel::Cylinder) : cylinder(cylinder) {}

std::string CylinderVolumePositionDistribution::Name() const {
    return "CylinderVolumePositionDistribution";
}

std::shared_ptr<InjectionDistribution> CylinderVolumePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new CylinderVolumePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> CylinderVolumePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & interaction) const {
    earthmodel::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pos(interaction.interaction_vertex);
    std::vector<earthmodel::Geometry::Intersection> intersections = cylinder.Intersections(pos, dir);
    earthmodel::EarthModel::SortIntersections(intersections);
    if(intersections.size() == 0) {
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));
    } else if(intersections.size() >= 2) {
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(intersections.front().position, intersections.back().position);
    } else {
        throw std::runtime_error("Only found one cylinder intersection!");
    }
}

bool CylinderVolumePositionDistribution::equal(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (cylinder == x->cylinder);
}

bool CylinderVolumePositionDistribution::less(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);
    return cylinder < x->cylinder;
}

//---------------
// class DepthFunction
//---------------
DepthFunction::DepthFunction() {}

double DepthFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool DepthFunction::operator==(DepthFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool DepthFunction::operator<(DepthFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class RangeFunction
//---------------
RangeFunction::RangeFunction() {}

double RangeFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool RangeFunction::operator==(RangeFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool RangeFunction::operator<(RangeFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class DecayRangeFunction
//---------------
//
//
double DecayRangeFunction::DecayLength(double particle_mass, double decay_width, double energy) {
    double beta = sqrt(energy*energy - particle_mass*particle_mass) / energy;
    double gamma = energy / particle_mass;
    double time_in_rest_frame = 1.0 / decay_width; // inverse GeV
    double time_in_lab_frame = time_in_rest_frame * gamma; // inverse GeV
    constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per inverse GeV
    double length = time_in_lab_frame * beta * iGeV_in_m; // meters = ((inverse GeV * dimensionless) * (meters per inverse GeV))
    return length; // meters
}

double DecayRangeFunction::DecayLength(InteractionSignature const & signature, double energy) const {
    return DecayRangeFunction::DecayLength(particle_mass, decay_width, energy);
}

double DecayRangeFunction::Range(InteractionSignature const & signature, double energy) const {
    return DecayLength(signature, energy) * multiplier;
}

double DecayRangeFunction::operator()(InteractionSignature const & signature, double energy) const {
    return Range(signature, energy);
}

double DecayRangeFunction::Multiplier() const {
    return multiplier;
}

double DecayRangeFunction::ParticleMass() const {
    return particle_mass;
}

double DecayRangeFunction::DecayWidth() const {
    return decay_width;
}

DecayRangeFunction::DecayRangeFunction(double particle_mass, double decay_width, double multiplier) : particle_mass(particle_mass), decay_width(decay_width), multiplier(multiplier) {}

bool DecayRangeFunction::equal(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(particle_mass, decay_width, multiplier)
            ==
            std::tie(x->particle_mass, x->decay_width, x->multiplier);
}

bool DecayRangeFunction::less(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    return
        std::tie(particle_mass, decay_width, multiplier)
        <
        std::tie(x->particle_mass, x->decay_width, x->multiplier);
}

//---------------
// class ColumnDepthPositionDistribution : VertexPositionDistribution
//---------------
earthmodel::Vector3D ColumnDepthPositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D ColumnDepthPositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    std::set<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections.GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double totalInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);
    double expTotalInteractionDepth = exp(totalInteractionDepth);

    double y = rand->Uniform();
    double traversedInteractionDepth = totalInteractionDepth - log(y + totalInteractionDepth * (1-y));

    double dist = path.GetDistanceFromStartAlongPath(traversedInteractionDepth, targets, total_cross_sections);
    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double ColumnDepthPositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    std::set<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections.GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double totalInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    path.SetPoints(path.GetFirstPoint(), vertex);

    double traversedInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    double prob_density = exp(-traversedInteractionDepth) / one_minus_exp_of_negative(totalInteractionDepth);
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

ColumnDepthPositionDistribution::ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), depth_function(depth_function), target_types(target_types) {}

std::string ColumnDepthPositionDistribution::Name() const {
    return "ColumnDepthPositionDistribution";
}

std::shared_ptr<InjectionDistribution> ColumnDepthPositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ColumnDepthPositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> ColumnDepthPositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool ColumnDepthPositionDistribution::equal(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (depth_function and x->depth_function and *depth_function == *x->depth_function)
                    or (!depth_function and !x->depth_function)
                )
            and target_types == x->target_types);
}

bool ColumnDepthPositionDistribution::less(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);
    bool depth_less =
        (!depth_function and x->depth_function) // this->NULL and other->(not NULL)
        or (depth_function and x->depth_function // both not NULL
                and *depth_function < *x->depth_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, depth_less, x->target_types);
}

//---------------
// class RangePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D RangePositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D RangePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    std::set<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections.GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double totalInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);
    double expTotalInteractionDepth = exp(totalInteractionDepth);

    double y = rand->Uniform();
    double traversedInteractionDepth = totalInteractionDepth - log(y + totalInteractionDepth * (1-y));

    double dist = path.GetDistanceFromStartAlongPath(traversedInteractionDepth, targets, total_cross_sections);
    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double RangePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    std::set<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections.GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double totalInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    path.SetPoints(path.GetFirstPoint(), vertex);

    double traversedInteractionDepth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    double prob_density = exp(-traversedInteractionDepth) / one_minus_exp_of_negative(totalInteractionDepth);
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

RangePositionDistribution::RangePositionDistribution() {}

RangePositionDistribution::RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {}

std::string RangePositionDistribution::Name() const {
    return "RangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> RangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new RangePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> RangePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));
    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool RangePositionDistribution::equal(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            and target_types == x->target_types);
}

bool RangePositionDistribution::less(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, range_less, x->target_types);
}

//---------------
// class DecayRangePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D DecayRangePositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D DecayRangePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    double y = rand->Uniform();
    double total_distance = path.GetDistance();
    double dist = -decay_length * log(y * (exp(-total_distance/decay_length) - 1) + 1);

    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double DecayRangePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    double total_distance = path.GetDistance();
    double dist = earthmodel::scalar_product(path.GetDirection(), vertex - path.GetFirstPoint());

    double prob_density = exp(-dist / decay_length) / (decay_length * (1.0 - exp(-total_distance / decay_length))); // m^-1
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3
    return prob_density;
}

DecayRangePositionDistribution::DecayRangePositionDistribution() {}

DecayRangePositionDistribution::DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {}

std::string DecayRangePositionDistribution::Name() const {
    return "DecayRangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> DecayRangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new DecayRangePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> DecayRangePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool DecayRangePositionDistribution::equal(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            and target_types == x->target_types);
}

bool DecayRangePositionDistribution::less(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, range_less, x->target_types);
}

//---------------
// class PrimaryNeutrinoHelicityDistribution : InjectionDistribution
//---------------
void PrimaryNeutrinoHelicityDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    std::array<double, 4> & mom = record.primary_momentum;
    double momentum = sqrt(mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3]);
    double factor = 0.5;
    Particle::ParticleType & t = record.signature.primary_type;
    if(t > 0) // Particles are left handed, anti-particles are right handed
        record.primary_helicity = -0.5;
    else
        record.primary_helicity = 0.5;
}

double PrimaryNeutrinoHelicityDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    std::array<double, 4> const & mom = record.primary_momentum;
    earthmodel::Vector3D dir(mom[1], mom[2], mom[3]);
    dir.normalize();

    if(abs(0.5 - abs(record.primary_helicity)) > 1e-9) // Helicity magnitude must be 0.5
        return 0.0;

    Particle::ParticleType const & t = record.signature.primary_type;
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
    const PrimaryNeutrinoHelicityDistribution* x = dynamic_cast<const PrimaryNeutrinoHelicityDistribution*>(&other);
    return false;
}

} // namespace LeptonInjector
