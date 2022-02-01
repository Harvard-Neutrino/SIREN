#include <cassert>
#include <fstream>
#include <algorithm>

#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/EventProps.h"

#include "phys-services/CrossSection.h"

#include "earthmodel-service/Path.h"

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

    // Find the HNL in the secondar particle vector and save its momentum/cartesian direction
    unsigned int lepton_index = (record.signature.secondary_types[0] == Particle::ParticleType::NuF4 or record.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
    std::array<double, 4> hnl_momentum = record.secondary_momenta[lepton_index];
    stga3::FourVector<double> pHNL_lab{hnl_momentum[0], hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]};
    double hnl_mass = std::sqrt(pHNL_lab | pHNL_lab);
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
    stga3::FourVector<double> pGamma_HNLrest{hnl_mass/2.0,
                                             hnl_mass/2.0*std::cos(phi)*std::sin(theta),
                                             hnl_mass/2.0*std::sin(phi)*std::sin(theta),
                                             hnl_mass/2.0*costh};

    // Boost gamma to lab frame
    stga3::Beta<double> beta_to_hnl_rest = stga3::beta_to_rest_frame_of(pHNL_lab);
    stga3::Boost<double> boost_to_lab = stga3::boost_from_beta(-beta_to_hnl_rest);
    stga3::FourVector<double> pGamma_lab = stga3::apply_boost(boost_to_lab,pGamma_HNLrest);

    std::array<double,4> gamma_momentum;
    gamma_momentum[0] = pGamma_lab.e0();
    gamma_momentum[1] = pGamma_lab.e1();
    gamma_momentum[2] = pGamma_lab.e2();
    gamma_momentum[3] = pGamma_lab.e3();
    record.secondary_momenta.push_back(gamma_momentum);
    record.signature.secondary_types.push_back(Particle::ParticleType::Gamma);
}

void InjectorBase::SamplePairProduction(InteractionRecord & record) {
    // Nick TODO: finish implementing this function which samples the photon pair produciton location
    earthmodel::Vector3D decay_vtx(record.decay_vertex[0],
                                   record.decay_vertex[1],
                                   record.decay_vertex[2]);
    unsigned int gamma_index = record.secondary_momenta.size() - 1;
    earthmodel::Vector3D decay_dir(record.secondary_momenta[gamma_index][1],
                                   record.secondary_momenta[gamma_index][2],
                                   record.secondary_momenta[gamma_index][3]);
    decay_dir.normalize();
    earthmodel::Path path(earth_model, decay_vtx, decay_dir, 0);
    path.ComputeIntersections();
    std::vector<double> P;
    double N = 0;
    earthmodel::Geometry::IntersectionList ilist = path.GetIntersections();
    std::cout << "HNL DECAY VERTEX: " << decay_vtx.GetX() << " " << decay_vtx.GetY() << " " << decay_vtx.GetZ() << std::endl;
    std::cout << "HNL DECAY DIRECTION: " << decay_dir.GetX() << " " << decay_dir.GetY() << " " << decay_dir.GetZ() << std::endl;
    std::cout << "INTERSECTIONS:\n";
    for(auto& i : ilist.intersections){
        std::cout << i.hierarchy << " " << i.distance << " " << i.position.GetX() << " " << i.position.GetY() << " " << i.position.GetZ() << std::endl;
    }
    std::cout << "\n";

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

