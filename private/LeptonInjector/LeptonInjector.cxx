#include <cassert>
#include <fstream>
#include <algorithm>

#include <rk/rk.hh>

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Distributions.h"

#include "LeptonInjector/LeptonInjector.h"


// namespace constants = boost::math::constants;

namespace LeptonInjector {

//---------------
// class InjectorBase
//---------------

InjectorBase::InjectorBase() {}

InjectorBase::InjectorBase(
        unsigned int events_to_inject,
        std::shared_ptr<PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::vector<std::shared_ptr<InjectionDistribution>> distributions,
        std::shared_ptr<LI_random> random) :
    events_to_inject(events_to_inject),
    primary_injector(primary_injector),
    cross_sections(std::make_shared<CrossSectionCollection>(primary_injector->PrimaryType(), cross_sections)),
    earth_model(earth_model),
    distributions(distributions),
    random(random)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        std::shared_ptr<PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random) :
    events_to_inject(events_to_inject),
    primary_injector(primary_injector),
    cross_sections(std::make_shared<CrossSectionCollection>(primary_injector->PrimaryType(), cross_sections)),
    earth_model(earth_model),
    random(random)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        std::shared_ptr<CrossSectionCollection> cross_sections) :
    events_to_inject(events_to_inject),
    cross_sections(cross_sections)
{}

InteractionRecord InjectorBase::NewRecord() const {
    InteractionRecord record;
    primary_injector->Sample(random, earth_model, cross_sections, record);
    return record;
}

void InjectorBase::SetRandom(std::shared_ptr<LI_random> random) {
    this->random = random;
}

void InjectorBase::SampleCrossSection(InteractionRecord & record) const {
    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();
    std::set<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
    std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);

    double total_prob = 0.0;
    std::vector<double> probs;
    std::vector<Particle::ParticleType> matching_targets;
    std::vector<InteractionSignature> matching_signatures;
    std::vector<std::shared_ptr<CrossSection>> matching_cross_sections;
    for(auto const target : available_targets) {
        if(possible_targets.find(target) != possible_targets.end()) {
            // Get target density
            double target_density = earth_model->GetParticleDensity(intersections, interaction_vertex, target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections->GetCrossSectionsForTarget(target);
            unsigned int xs_i = 0;
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                unsigned int sig_i = 0;
                for(auto const & signature : signatures) {
                    record.signature = signature;
                    record.target_mass = earth_model->GetTargetMass(target);
                    record.target_momentum = {record.target_mass,0,0,0};
                    // Add total cross section times density to the total prob
                    total_prob += target_density * cross_section->TotalCrossSection(record);
                    // Add total prob to probs
                    probs.push_back(total_prob);
                    // Add target and cross section pointer to the lists
                    matching_targets.push_back(target);
                    matching_cross_sections.push_back(cross_section);
                    matching_signatures.push_back(signature);
                    sig_i += 1;
                }
                xs_i += 1;
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

double InjectorBase::CrossSectionProbability(InteractionRecord const & record) const {
    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();
    std::set<Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
    std::set<Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    earthmodel::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    earthmodel::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);

    double total_prob = 0.0;
    double selected_prob = 0.0;
    for(auto const target : available_targets) {
        if(possible_targets.find(target) != possible_targets.end()) {
            // Get target density
            double target_density = earth_model->GetParticleDensity(intersections, interaction_vertex, target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<CrossSection>> const & target_cross_sections = cross_sections->GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<InteractionSignature> signatures = cross_section->GetPossibleSignatures();
                for(auto const & signature : signatures) {
                    // Add total cross section times density to the total prob
                    double target_prob = target_density * cross_section->TotalCrossSection(record);
                    total_prob += target_prob;
                    // Add up total cross section times density times final state prob for matching signatures
                    if(signature == record.signature) {
                        selected_prob += target_prob * cross_section->FinalStateProbability(record);
                    }
                }
            }
        }
    }
    return selected_prob / total_prob;
}

void InjectorBase::SampleSecondaryDecay(InteractionRecord const & interaction, DecayRecord & decay, double decay_width) const {
    // This function takes an interaction record containing an HNL and simulates the decay to a photon
    // Currently assumes Majorana HNL
    // Final state photon added to secondary particle vectors in InteractionRecord

    // Find the HNL in the secondary particle vector and save its momentum/cartesian direction
    unsigned int lepton_index = (interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == Particle::ParticleType::NuF4Bar) ? 0 : 1;
    Particle::ParticleType hnl_type = interaction.signature.secondary_types[lepton_index];
    double hnl_mass = interaction.secondary_masses[lepton_index];
    std::array<double, 4> hnl_momentum = interaction.secondary_momenta[lepton_index];
    double hnl_helicity = interaction.secondary_helicity[lepton_index];

    // Store the HNL as the primary for the decay
    decay.signature.primary_type = hnl_type;
    decay.primary_mass = hnl_mass;
    decay.primary_momentum = hnl_momentum;
    decay.primary_helicity = hnl_helicity;

    rk::P4 pHNL_lab(geom3::Vector3(hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]), hnl_mass);
    earthmodel::Vector3D hnl_dir(hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]);
    hnl_dir.normalize();

    // Calculate the decay location of the HNL
    double decay_length = LeptonInjector::DecayRangeFunction::DecayLength(hnl_mass, decay_width, hnl_momentum[0]);
    double decay_loc = -1 * decay_length * std::log(random->Uniform(0,1));
    decay.decay_vertex = earthmodel::Vector3D(interaction.interaction_vertex) + decay_loc * hnl_dir;

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

    decay.signature.secondary_types.resize(1);
    decay.secondary_masses.resize(1);
    decay.secondary_momenta.resize(1);
    decay.secondary_helicity.resize(1);

    decay.signature.secondary_types[0] = Particle::ParticleType::Gamma;
    decay.secondary_masses[0] = 0;
    decay.secondary_momenta[0][0] = pGamma_lab.e();
    decay.secondary_momenta[0][1] = pGamma_lab.px();
    decay.secondary_momenta[0][2] = pGamma_lab.py();
    decay.secondary_momenta[0][3] = pGamma_lab.pz();
    decay.secondary_helicity[0] = std::copysign(1.0, decay.primary_helicity);

    decay.decay_parameters.resize(1);
    decay.decay_parameters[0] = decay_length;
}

void InjectorBase::SamplePairProduction(DecayRecord const & decay, InteractionRecord & interaction) const {
    // function for simulating the pair production of the photon created in HNL decay
    // considers the different radiation lengths of materials in the detector
    // Nick TODO: comment more

    earthmodel::MaterialModel const & mat_model = earth_model->GetMaterials();
    earthmodel::Vector3D decay_vtx(decay.decay_vertex);
    unsigned int gamma_index = 0;
    earthmodel::Vector3D decay_dir(decay.secondary_momenta[gamma_index][1],
                                   decay.secondary_momenta[gamma_index][2],
                                   decay.secondary_momenta[gamma_index][3]);

    interaction.signature.primary_type = decay.signature.secondary_types[gamma_index];
    interaction.primary_mass = decay.secondary_masses[gamma_index];
    interaction.primary_momentum = decay.secondary_momenta[gamma_index];
    interaction.primary_helicity = decay.secondary_helicity[gamma_index];

    decay_dir.normalize();

    earthmodel::Path path(earth_model, decay_vtx, decay_dir, 0);
    path.ComputeIntersections();

    std::vector<double> X0;
    std::vector<double> P;
    std::vector<double> D;
    double x0; double p; double density;
    D.push_back(0.);
    double N = 0;
    double lnP_nopp = 0; // for calculating the probability that no pair production occurs
    earthmodel::Geometry::IntersectionList const & ilist = path.GetIntersections();
    earthmodel::Vector3D density_point = decay_vtx;
    int j = 0;
    for(auto const & i : ilist.intersections) {
        if(i.distance<0 || std::isinf(i.distance))
            continue;
        D.push_back(i.distance);
        x0 = (9./7.)*mat_model.GetMaterialRadiationLength(i.matID); // in g/cm^2
        density_point += 0.5*(i.position - density_point);
        density = earth_model->GetMassDensity(density_point);
        x0 *= 0.01/density; // in m
        X0.push_back(x0);
        p = std::exp(-D[j]/x0) - std::exp(-D[j+1]/x0);
        P.push_back(p);
        N += p;
        lnP_nopp += -(D[j+1] - D[j])/x0;
        density_point = i.position;
        ++j;
    }

    interaction.interaction_parameters.resize(1);
    interaction.interaction_parameters[0] = std::exp(lnP_nopp);

    // sample the PDF by inverting the CDF
    double X = random->Uniform(0, 1);
    double C = 0;
    if(P.size() > 0) {
        for(j = 0; j < P.size(); ++j){
            C += P[j]/N;
            if(C>X) {C -= P[j]/N; break;}
        }
        double pairprod_dist = -X0[j]*std::log(X - C + std::exp(-D[j]/X0[j]));
        interaction.interaction_vertex = earthmodel::Vector3D(decay.decay_vertex) + pairprod_dist * decay_dir;
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

double InjectorBase::GenerationProbability(InteractionRecord const & record) const {
    double probability = 1.0;
    for(auto const & dist : distributions) {
        probability *= dist->GenerationProbability(earth_model, cross_sections, record);
    }
    probability *= CrossSectionProbability(record);
    probability *= events_to_inject;
    return probability;
}

std::set<std::vector<std::string>> InjectorBase::DensityVariables() const {
    std::set<std::vector<std::string>> variable_sets;
    std::vector<std::string> variables;
    for(auto const & dist : distributions) {
        std::vector<std::string> new_variables = dist->DensityVariables();
        variables.reserve(variables.size() + new_variables.size());
        variables.insert(variables.end(), new_variables.begin(), new_variables.end());
    }
    std::vector<std::shared_ptr<CrossSection>> xs_vec = cross_sections->GetCrossSections();
    for(auto const & xs : xs_vec) {
        std::vector<std::string> new_variables = xs->DensityVariables();
        std::vector<std::string> variable_list;
        variable_list.reserve(variables.size() + new_variables.size());
        variable_list.insert(variable_list.end(), variables.begin(), variables.end());
        variable_list.insert(variable_list.end(), new_variables.begin(), new_variables.end());
        variable_sets.insert(variable_list);
    }
    return variable_sets;
}

std::string InjectorBase::Name() const {
    return("InjectorBase");
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> InjectorBase::InjectionBounds(InteractionRecord const & interaction) const {
    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));
}

std::vector<std::shared_ptr<InjectionDistribution>> InjectorBase::GetInjectionDistributions() const {
    return distributions;
}

std::shared_ptr<earthmodel::EarthModel> InjectorBase::GetEarthModel() const {
    return earth_model;
}

std::shared_ptr<CrossSectionCollection> InjectorBase::GetCrossSections() const {
    return cross_sections;
}

unsigned int InjectorBase::InjectedEvents() const {
    return injected_events;
}

unsigned int InjectorBase::EventsToInject() const {
    return events_to_inject;
}

InjectorBase::operator bool() const {
    return injected_events < events_to_inject;
}

//---------------
// class RangedLeptonInjector : InjectorBase
//---------------
RangedLeptonInjector::RangedLeptonInjector() {}

RangedLeptonInjector::RangedLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        std::shared_ptr<RangeFunction> range_func,
        double disk_radius,
        double endcap_length,
        std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    disk_radius(disk_radius),
    endcap_length(endcap_length),
    helicity_distribution(helicity_distribution),
    InjectorBase(events_to_inject, primary_injector, cross_sections, earth_model, random)
{
    std::set<Particle::ParticleType> target_types = this->cross_sections->TargetTypes();
    position_distribution = std::make_shared<RangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
    distributions = {target_momentum_distribution, energy_distribution, helicity_distribution, direction_distribution, position_distribution};
}

InteractionRecord RangedLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose the helicity
    helicity_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);

    injected_events += 1;
    return event;
}

std::string RangedLeptonInjector::Name() const {
    return("RangedInjector");
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> RangedLeptonInjector::InjectionBounds(InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

//---------------
// class DecayRangeLeptonInjector : InjectorBase
//---------------
DecayRangeLeptonInjector::DecayRangeLeptonInjector() {}

DecayRangeLeptonInjector::DecayRangeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        std::shared_ptr<DecayRangeFunction> range_func,
        double disk_radius,
        double endcap_length,
        std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    disk_radius(disk_radius),
    endcap_length(endcap_length),
    helicity_distribution(helicity_distribution),
    InjectorBase(events_to_inject, primary_injector, cross_sections, earth_model, random)
{
    std::set<Particle::ParticleType> target_types = this->cross_sections->TargetTypes();
    position_distribution = std::make_shared<DecayRangePositionDistribution>(disk_radius, endcap_length, range_func, target_types);
    distributions = {target_momentum_distribution, energy_distribution, helicity_distribution, direction_distribution, position_distribution};
}

InteractionRecord DecayRangeLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose the helicity
    helicity_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a direction on the sphere
    direction_distribution->Sample(random, earth_model, cross_sections, event);

    // Pick a position for the vertex
    position_distribution->Sample(random, earth_model, cross_sections, event);

    // Sample the cross section and final state
    SampleCrossSection(event);

    injected_events += 1;
    return event;
}

std::string DecayRangeLeptonInjector::Name() const {
    return("DecayRangeInjector");
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> DecayRangeLeptonInjector::InjectionBounds(InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

//---------------
// class VolumeLeptonInjector : InjectorBase
//---------------
VolumeLeptonInjector::VolumeLeptonInjector() {}

VolumeLeptonInjector::VolumeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<CrossSection>> cross_sections,
        std::shared_ptr<earthmodel::EarthModel> earth_model,
        std::shared_ptr<LI_random> random,
        std::shared_ptr<PrimaryEnergyDistribution> edist,
        std::shared_ptr<PrimaryDirectionDistribution> ddist,
        std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution,
        earthmodel::Cylinder cylinder,
        std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution) :
    energy_distribution(edist),
    direction_distribution(ddist),
    target_momentum_distribution(target_momentum_distribution),
    position_distribution(std::make_shared<CylinderVolumePositionDistribution>(cylinder)),
    helicity_distribution(helicity_distribution),
    InjectorBase(events_to_inject, primary_injector, cross_sections, earth_model, random) {
    distributions = {target_momentum_distribution, energy_distribution, helicity_distribution, direction_distribution, position_distribution};
}

InteractionRecord VolumeLeptonInjector::GenerateEvent() {
    InteractionRecord event = NewRecord();

    // Choose a target momentum
    target_momentum_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose an energy
    energy_distribution->Sample(random, earth_model, cross_sections, event);

    // Choose the helicity
    helicity_distribution->Sample(random, earth_model, cross_sections, event);

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

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> VolumeLeptonInjector::InjectionBounds(InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

} // namespace LeptonInjector

