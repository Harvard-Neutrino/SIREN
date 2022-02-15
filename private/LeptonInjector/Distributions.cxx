#include "LeptonInjector/Distributions.h"
#include "earthmodel-service/EarthModelCalculator.h"

namespace LeptonInjector {

//---------------
// class WeightableDistribution
//---------------

std::vector<std::string> WeightableDistribution::DensityVariables() const {
    return {};
}

//---------------
// class InjectionDistribution
//---------------

void InjectionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
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

//---------------
// class ArbPDF : PrimaryEnergyDistribution
//---------------
ArbPDF::ArbPDF(double minE_, double maxE_, std::vector<double> params_, double (*PDF_)(double, std::vector<double>))
    : PDF(PDF_)
    , minE(minE_)
    , maxE(maxE_)
    , params(params_)
{
    PDF = PDF_;
    minE = minE_;
    maxE = maxE_;
    params = params_;
    std::function<double(double)> integrand = [&] (double x) -> double {
        return (*PDF)(x, params);
    };
    integral = earthmodel::Integration::rombergIntegrate(integrand, minE, maxE);
}

double ArbPDF::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    //Nick TODO: implement Metropolis-Hastings algorithm to sample from PDF. 
    // Pass in a function pointer for the PDF
    
    double E, testE, odds;
    bool accept;
    
    // sample an initial point uniformly
    E = rand->Uniform(minE,maxE);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        accept = true;
        testE = rand->Uniform(minE,maxE);
        odds = (*PDF)(testE,params)/(*PDF)(E,params);
        accept = (odds > 1.) || rand->Uniform(0,1)<odds;
        if(accept) E = testE;
    }

    return E;
}

double ArbPDF::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {

    double density = (*PDF)(record.primary_momentum[0], params);
    return density / integral;
}

std::string ArbPDF::Name() const {
    return "ArbPDF";
}

std::shared_ptr<InjectionDistribution> ArbPDF::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ArbPDF(*this));
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

//---------------
// class DepthFunction
//---------------
DepthFunction::DepthFunction() {}

double DepthFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

//---------------
// class RangeFunction
//---------------
RangeFunction::RangeFunction() {}

double RangeFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

//---------------
// class DecayRangeFunction
//---------------
//
//
double DecayRangeFunction::DecayLength(InteractionSignature const & signature, double energy) const {
    std::array<double, 4> lab_momentum{energy, 0.0, 0.0, sqrt(energy*energy - particle_mass*particle_mass)}; // GeV
    double beta = sqrt((lab_momentum[1]*lab_momentum[1] + lab_momentum[2]*lab_momentum[2] + lab_momentum[3]*lab_momentum[3]) / (lab_momentum[0]*lab_momentum[0])); // dimensionless
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    double time_in_rest_frame = 1.0 / decay_width; // inverse GeV
    double time_in_lab_frame = time_in_rest_frame * gamma; // inverse GeV
    constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per inverse GeV
    double length = time_in_lab_frame * beta * iGeV_in_m; // meters = ((inverse GeV * dimensionless) * (meters per inverse GeV))
    return length; // meters
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
    path.ExtendFromStartByColumnDepth(lepton_depth, target_types);
    path.ClipToOuterBounds();

    double totalColumnDepth = path.GetColumnDepthInBounds(target_types);

    double traversedColumnDepth = totalColumnDepth * rand->Uniform();
    double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth);
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
    path.ExtendFromStartByColumnDepth(lepton_depth, target_types);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    double totalColumnDepth = path.GetColumnDepthInBounds(target_types); // g/cm^2
    double density = earth_model->GetDensity(vertex, target_types); // g/cm^3
    double prob_density = density / totalColumnDepth * 100; // (cm^-1 * cm/m) -> m^-1
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

    double totalColumnDepth = path.GetColumnDepthInBounds(target_types);

    double traversedColumnDepth = totalColumnDepth * rand->Uniform();
    double dist = path.GetDistanceFromStartAlongPath(traversedColumnDepth, target_types);
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

    double totalColumnDepth = path.GetColumnDepthInBounds(target_types); // g/cm^2
    double density = earth_model->GetDensity(vertex, target_types); // g/cm^3
    double prob_density = density / totalColumnDepth * 100; // (cm^-1 * cm/m) -> m^-1
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

//---------------
// class PrimaryNeutrinoSpinDistribution : InjectionDistribution
//---------------
void PrimaryNeutrinoSpinDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord & record) const {
    std::array<double, 4> & mom = record.primary_momentum;
    double momentum = sqrt(mom[1]*mom[1] + mom[2]*mom[2] + mom[3]*mom[3]);
    double factor = 0.5 / momentum;
    Particle::ParticleType & t = record.signature.primary_type;
    if(t > 0) // Particles are left handed, anti-particles are right handed
        factor = -factor;
    record.primary_spin[0] = factor * mom[1];
    record.primary_spin[1] = factor * mom[2];
    record.primary_spin[2] = factor * mom[3];
}

double PrimaryNeutrinoSpinDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel> earth_model, CrossSectionCollection const & cross_sections, InteractionRecord const & record) const {
    std::array<double, 4> const & mom = record.primary_momentum;
    earthmodel::Vector3D dir(mom[1], mom[2], mom[3]);
    dir.normalize();
    earthmodel::Vector3D spin_dir(record.primary_spin);
    double spin_magnitude = spin_dir.magnitude();
    spin_dir.normalize();

    double dot = earthmodel::scalar_product(dir, spin_dir);

    if(abs(0.5 - spin_magnitude) > 1e-9) // Spin magnitude must be 0.5
        return 0.0;

    Particle::ParticleType const & t = record.signature.primary_type;
    // Particles are left handed, anti-particles are right handed
    if(t > 0) {
        if(abs(1.0 + dot) > 1e-9) // expect opposite direction
            return 0.0;
    } else {
        if(abs(1.0 - dot) > 1e-9) // expect same direction
            return 0.0;
    }
    return 1.0;
}

PrimaryNeutrinoSpinDistribution::PrimaryNeutrinoSpinDistribution() {}

std::vector<std::string> PrimaryNeutrinoSpinDistribution::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryNeutrinoSpinDistribution::Name() const {
    return "PrimaryNeutrinoSpinDistribution";
}

std::shared_ptr<InjectionDistribution> PrimaryNeutrinoSpinDistribution::clone() const {
    return std::shared_ptr<PrimaryNeutrinoSpinDistribution>(new PrimaryNeutrinoSpinDistribution(*this));
}

} // namespace LeptonInjector
