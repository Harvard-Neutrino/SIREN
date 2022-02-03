#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"

#include <rk/rk.hh>

// namespace constants = boost::math::constants;

namespace LeptonInjector{

//---------------
// class InjectorBase
//---------------

InjectorBase::InjectorBase(
        unsigned int events_to_inject,
        Particle::ParticleType primary_type,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::vector<std::shared_ptr<InjectionDistribution>> distributions,
        std::shared_ptr<LI_random> random) :
    events_to_inject(events_to_inject),
    primary_type(primary_type),
    cross_sections(primary_type, cross_sections),
    earth_model(earth_model),
    distributions(distributions),
    random(random)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        Particle::ParticleType primary_type,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random) :
    events_to_inject(events_to_inject),
    primary_type(primary_type),
    cross_sections(primary_type, cross_sections),
    earth_model(earth_model),
    random(random)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        CrossSectionCollection cross_sections) :
    events_to_inject(events_to_inject),
    cross_sections(cross_sections)
{}

InteractionRecord InjectorBase::NewRecord() const {
    InteractionRecord record;
    record.signature.primary_type = primary_type;
    record.primary_mass = Particle(primary_type).GetMass();
    return record;
}

void InjectorBase::SetRandom(std::shared_ptr<LI_random> random) {
    this->random = random;
}

void InjectorBase::SampleCrossSection(InteractionRecord & record) const {
    std::vector<Particle::ParticleType> const & possible_targets = cross_sections.TargetTypes();
    std::vector<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
    std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    earthmodel::Vector3D intVertex(record.interaction_vertex[0],record.interaction_vertex[1],record.interaction_vertex[2]);

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(intVertex, earthmodel::Vector3D());

    double total_prob = 0.0;
    std::vector<Particle::ParticleType> single;
    std::vector<double> probs;
    std::vector<Particle::ParticleType> matching_targets;
    std::vector<InteractionSignature> matching_signatures;
    std::vector<std::shared_ptr<CrossSection>> matching_cross_sections;
    for(auto const target : possible_targets) {
        if(available_targets.find(target) != available_targets.end()) {
            // Get target density
            single.push_back(target);
            double target_density = earth_model->GetDensity(intVertex, single);
            single.clear();
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections.GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<InteractionSignature> signatures = cross_section->GetPossibleSignatures();
                for(auto const & signature : signatures) {
                    record.signature = signature;
                    // Add total cross section times density to the total prob
                    total_prob += target_density * cross_section->TotalCrossSection(record);
                    // Add total prob to probs
                    probs.push_back(total_prob);
                    // Add target and cross section pointer to the lists
                    matching_targets.push_back(target);
                    matching_cross_sections.push_back(cross_section);
                    matching_signatures.push_back(signature);
                }
            }
        }
    }
    // Throw a random number
    double r = random->Uniform(0, total_prob);
    // Choose the target and cross section
    unsigned int index = 0;
    for(; (index < probs.size()-1) and (r > probs[index]); ++index) {}
    record.signature.target_type = matching_targets[index];
    record.signature = matching_signatures[index];
    record.target_mass = earth_model->GetTargetMass(record.signature.target_type);
    record.target_momentum = {record.target_mass,0,0,0};
    matching_cross_sections[index]->SampleFinalState(record, random);
}

void InjectorBase::SampleSecondaryDecay(InteractionRecord & record) const {
    // This function takes an interaction record containing an HNL and simulates the decay to a photon
    // Currently assumes Majorana HNL
    // Final state photon added to secondary particle vectors in InteractionRecord

    // Find the HNL in the secondary particle vector and save its momentum/cartesian direction
    unsigned int lepton_index = (record.signature.secondary_types[0] == Particle::ParticleType::NuF4 or record.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
    std::array<double, 4> hnl_momentum = record.secondary_momenta[lepton_index];
    // stga3::FourVector<double> pHNL_lab{hnl_momentum[0], hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]};
    // double hnl_mass = std::sqrt(pHNL_lab | pHNL_lab);
    double hnl_mass = record.secondary_masses[lepton_index];
    rk::P4 pHNL_lab(geom3::Vector3(hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]), hnl_mass);
    earthmodel::Vector3D hnl_dir(hnl_momentum[1],hnl_momentum[2],hnl_momentum[3]);
    hnl_dir.normalize();

    // Calculate the decay location of the HNL
    double decay_loc = -1 * record.decay_length * std::log(random->Uniform(0,1));
    record.decay_vertex = {record.interaction_vertex[0] + decay_loc*hnl_dir.GetX(),
                           record.interaction_vertex[1] + decay_loc*hnl_dir.GetY(),
                           record.interaction_vertex[2] + decay_loc*hnl_dir.GetZ()};

    // Sample decay angles
    // Majorana Case: Isotropic Decay
    double costh = random->Uniform(-1,1);
    double theta = std::acos(costh);
    double phi = random->Uniform(0,2*Constants::pi);
    rk::P4 pGamma_HNLrest(
            geom3::Vector3(
                hnl_mass/2.0*std::cos(phi)*std::sin(theta),
                hnl_mass/2.0*std::sin(phi)*std::sin(theta),
                hnl_mass/2.0*costh),
            0.0);

    // Boost gamma to lab frame
    rk::Boost boost_to_lab = pHNL_lab.labBoost();
    rk::P4 pGamma_lab = boost_to_lab * pGamma_HNLrest;

    std::array<double,4> gamma_momentum;
    gamma_momentum[0] = pGamma_lab.e();
    gamma_momentum[1] = pGamma_lab.px();
    gamma_momentum[2] = pGamma_lab.py();
    gamma_momentum[3] = pGamma_lab.pz();
    record.secondary_momenta.push_back(gamma_momentum);
    record.secondary_masses.push_back(0);
    record.signature.secondary_types.push_back(Particle::ParticleType::Gamma);
}

void InjectorBase::SamplePairProduction(InteractionRecord & record) {
    // function for simulating the pair production of the photon created in HNL decay
    // considers the different radiation lengths of materials in the detector
    // Nick TODO: comment more

    earthmodel::MaterialModel mat_model = earth_model->GetMaterials();
    earthmodel::Vector3D decay_vtx(record.decay_vertex[0],
                                   record.decay_vertex[1],
                                   record.decay_vertex[2]);
    unsigned int gamma_index = record.secondary_momenta.size() - 1;
    earthmodel::Vector3D decay_dir(record.secondary_momenta[gamma_index][1],
                                   record.secondary_momenta[gamma_index][2],
                                   record.secondary_momenta[gamma_index][3]);
    decay_dir.normalize();
    if(std::isnan(decay_dir.magnitude())){
        for(int j = 0; j < record.secondary_momenta.size(); ++j)
        {
        std::cout << j << " ";
        for(int k = 0; k < 4; ++k) std::cout << record.secondary_momenta[j][k] << " ";
        std::cout << std::endl;
        }
        return;
    }
    earthmodel::Path path(earth_model, decay_vtx, decay_dir, 0);
    path.ComputeIntersections();
    std::vector<double> X0;
    std::vector<double> P;
    std::vector<double> D;
    double x0; double p; double density;
    D.push_back(0.);
    double N = 0;
    double lnP_nopp = 0; // for calculating the probability that no pair production occurs
    earthmodel::Geometry::IntersectionList ilist = path.GetIntersections();
    earthmodel::Vector3D density_point = decay_vtx;
    int j = 0;
    for(auto& i : ilist.intersections){
        if(i.distance<0 || std::isinf(i.distance)) continue;
        D.push_back(i.distance);
        x0 = (9./7.)*mat_model.GetMaterialRadLength(i.matID); // in g/cm^2
        density_point += 0.5*(i.position - density_point);
        density = earth_model->GetDensity(density_point);
        x0 *= 0.01/density; // in m
        X0.push_back(x0);
        p = std::exp(-D[j]/x0) - std::exp(-D[j+1]/x0);
        P.push_back(p);
        N += p;
        lnP_nopp += -(D[j+1] - D[j])/x0;
        density_point = i.position;
        ++j;
    }
    record.prob_nopairprod = std::exp(lnP_nopp);


    // sample the PDF by inverting the CDF
    double X = random->Uniform(0,1);
    double C = 0;
    if(P.size() > 0) {
        for(j = 0; j < P.size(); ++j){
            C += P[j]/N;
            if(C>X) {C -= P[j]/N; break;}
        }
        double pairprod_dist = -X0[j]*std::log(X - C + std::exp(-D[j]/X0[j]));
        record.pairprod_vertex = {record.decay_vertex[0] + pairprod_dist*decay_dir.GetX(),
                                  record.decay_vertex[1] + pairprod_dist*decay_dir.GetY(),
                                  record.decay_vertex[2] + pairprod_dist*decay_dir.GetZ()};
    }

}

InteractionRecord InjectorBase::GenerateEvent() {
    InteractionRecord record = this->NewRecord();
    for(auto & distribution : distributions) {
        distribution->Sample(random, earth_model, cross_sections, record);
    }
    SampleCrossSection(record);
    injected_events += 1;
    return record;
}

std::string InjectorBase::Name() const {
    return("InjectorBase");
}

InjectorBase::operator bool() const {
    return injected_events < events_to_inject;
}

//---------------
// class RangedLeptonInjector : InjectorBase
//---------------
RangedLeptonInjector::RangedLeptonInjector(
        unsigned int events_to_inject,
        Particle::ParticleType primary_type,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        std::shared_ptr<RangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    disk_radius(disk_radius),
    endcap_length(endcap_length),
    InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random)
{
    std::vector<Particle::ParticleType> target_types = this->cross_sections.TargetTypes();
    position_distribution = std::make_shared<RangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
}

InteractionRecord RangedLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);

    // Sample decay angle of photon
    SampleSecondaryDecay(event);

    // Sample pair production location
    SamplePairProduction(event);

    injected_events += 1;
    return event;
}

std::string RangedLeptonInjector::Name() const {
    return("RangedInjector");
}

//---------------
// class DecayRangeLeptonInjector : InjectorBase
//---------------
DecayRangeLeptonInjector::DecayRangeLeptonInjector(
        unsigned int events_to_inject,
        Particle::ParticleType primary_type,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        std::shared_ptr<DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    disk_radius(disk_radius),
    endcap_length(endcap_length),
    InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random)
{
    std::vector<Particle::ParticleType> target_types = this->cross_sections.TargetTypes();
    position_distribution = std::make_shared<DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
}

InteractionRecord DecayRangeLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);

    // Sample decay angle of photon
    SampleSecondaryDecay(event);

    // Sample pair production location
    SamplePairProduction(event);

    injected_events += 1;
    return event;
}

std::string DecayRangeLeptonInjector::Name() const {
    return("DecayRangeInjector");
}

//---------------
// class VolumeLeptonInjector : InjectorBase
//---------------
VolumeLeptonInjector::VolumeLeptonInjector(
        unsigned int events_to_inject,
        Particle::ParticleType primary_type,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        earthmodel::Cylinder cylinder) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    position_distribution(std::make_shared<CylinderVolumePositionDistribution>(cylinder)),
    InjectorBase(events_to_inject, primary_type, cross_sections, earth_model, random)
{}

InteractionRecord VolumeLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);
    injected_events += 1;
    return event;
}

std::string VolumeLeptonInjector::Name() const {
    return("VolumeInjector");
}

} // namespace LeptonInjector

