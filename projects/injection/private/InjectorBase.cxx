#include <cassert>
#include <fstream>
#include <algorithm>

#include <rk/rk.hh>

#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/injection/Weighter.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"

// For CrossSectionProbability
#include "LeptonInjector/injection/WeightingUtils.h"

#include "LeptonInjector/injection/InjectorBase.h"

#include "LeptonInjector/utilities/Errors.h"

namespace LI {
namespace injection {

//---------------
// class InjectorBase
//---------------

InjectorBase::InjectorBase() {}

InjectorBase::InjectorBase(
        unsigned int events_to_inject,
        std::shared_ptr<LI::distributions::PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> distributions,
        std::shared_ptr<LI::utilities::LI_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    primary_injector(primary_injector),
    cross_sections(std::make_shared<LI::crosssections::CrossSectionCollection>(primary_injector->PrimaryType(), cross_sections)),
    earth_model(earth_model),
    distributions(distributions)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        std::shared_ptr<LI::distributions::PrimaryInjector> primary_injector,
        std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::shared_ptr<LI::utilities::LI_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    primary_injector(primary_injector),
    cross_sections(std::make_shared<LI::crosssections::CrossSectionCollection>(primary_injector->PrimaryType(), cross_sections)),
    earth_model(earth_model)
{}

InjectorBase::InjectorBase(unsigned int events_to_inject,
        std::shared_ptr<LI::crosssections::CrossSectionCollection> cross_sections) :
    events_to_inject(events_to_inject),
    cross_sections(cross_sections)
{}

LI::dataclasses::InteractionRecord InjectorBase::NewRecord() const {
    LI::dataclasses::InteractionRecord record;
    primary_injector->Sample(random, earth_model, cross_sections, record);
    return record;
}

void InjectorBase::SetRandom(std::shared_ptr<LI::utilities::LI_random> random) {
    this->random = random;
}

void InjectorBase::SampleCrossSection(LI::dataclasses::InteractionRecord & record) const {

    // Make sure the particle has interacted
    if(std::isnan(record.interaction_vertex[0]) ||
       std::isnan(record.interaction_vertex[1]) ||
       std::isnan(record.interaction_vertex[2])) {
	    throw(LI::utilities::InjectionFailure("No particle interaction!"));
    }

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();
    std::set<LI::dataclasses::Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(record.interaction_vertex);
    std::set<LI::dataclasses::Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();


    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(interaction_vertex, primary_direction);

    double total_prob = 0.0;
    std::vector<double> probs;
    std::vector<LI::dataclasses::Particle::ParticleType> matching_targets;
    std::vector<LI::dataclasses::InteractionSignature> matching_signatures;
    std::vector<std::shared_ptr<LI::crosssections::CrossSection>> matching_cross_sections;
    LI::dataclasses::InteractionRecord fake_record = record;
    for(auto const target : available_targets) {
        if(possible_targets.find(target) != possible_targets.end()) {
            // Get target density
            double target_density = earth_model->GetParticleDensity(intersections, interaction_vertex, target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<LI::crosssections::CrossSection>> const & target_cross_sections = cross_sections->GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<LI::dataclasses::InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                for(auto const & signature : signatures) {
                    fake_record.signature = signature;
                    fake_record.target_mass = earth_model->GetTargetMass(target);
                    fake_record.target_momentum = {fake_record.target_mass,0,0,0};
                    // Add total cross section times density to the total prob
                    total_prob += target_density * cross_section->TotalCrossSection(fake_record);
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
    double selected_prob = 0.0;
    for(unsigned int i=0; i<probs.size(); ++i) {
        if(matching_signatures[index] == matching_signatures[i]) {
            selected_prob += (i > 0 ? probs[i] - probs[i - 1] : probs[i]);
        }
    }
    if(total_prob == 0 or selected_prob == 0)
        throw(LI::utilities::InjectionFailure("No valid interactions for this event!"));
    record.target_mass = earth_model->GetTargetMass(record.signature.target_type);
    record.target_momentum = {record.target_mass,0,0,0};
    matching_cross_sections[index]->SampleFinalState(record, random);
}

void InjectorBase::SampleSecondaryDecay(LI::dataclasses::InteractionRecord const & interaction, LI::dataclasses::DecayRecord & decay, double decay_width, double alpha_gen, double alpha_phys, LI::geometry::Geometry *fiducial = nullptr, double buffer = 0) const {
    // This function takes an interaction record containing an HNL and simulates the decay to a photon
    // Samples according to (1 + alpha * cos(theta))/2 and returns physical weight
    // Final state photon added to secondary particle vectors in LI::dataclasses::InteractionRecord

    // Find the HNL in the secondary particle vector and save its momentum/cartesian direction
    unsigned int lepton_index = (interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4Bar) ? 0 : 1;
    LI::dataclasses::Particle::ParticleType hnl_type = interaction.signature.secondary_types[lepton_index];
    double hnl_mass = interaction.secondary_masses[lepton_index];
    std::array<double, 4> hnl_momentum = interaction.secondary_momenta[lepton_index];
    double hnl_helicity = interaction.secondary_helicity[lepton_index];

    // Store the HNL as the primary for the decay
    decay.signature.primary_type = hnl_type;
    decay.primary_mass = hnl_mass;
    decay.primary_momentum = hnl_momentum;
    decay.primary_helicity = hnl_helicity;

    rk::P4 pHNL_lab(geom3::Vector3(hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]), hnl_mass);
    LI::math::Vector3D hnl_dir(hnl_momentum[1], hnl_momentum[2], hnl_momentum[3]);
    hnl_dir.normalize();

    // Calculate the decay location of the HNL
    // Require the decay to happen within a fid vol if possible
    double decay_length = LI::distributions::DecayRangeFunction::DecayLength(hnl_mass, decay_width, hnl_momentum[0]);
    double decay_weight = 1.0;
    double a=0,b=0;
    double C = random->Uniform(0,1);
    if(fiducial) {
				std::vector<LI::geometry::Geometry::Intersection> ints = fiducial->Intersections(interaction.interaction_vertex,hnl_dir);
				if(ints.size()!=0 && ints[ints.size()-1].distance > 0) {
						a = std::max(0.,ints[0].distance - buffer);
						b = ints[ints.size()-1].distance;
						C*=(1-std::exp(-(b-a)/decay_length));
						decay_weight = std::exp(-a/decay_length) - std::exp(-b/decay_length);
				}

    }

    double decay_loc = a + -1 * decay_length * std::log(1-C);
    decay.decay_vertex = LI::math::Vector3D(interaction.interaction_vertex) + decay_loc * hnl_dir;

    // Sample decay angles
    double X = random->Uniform(0,1);
    double costh;
    // Majorana Case
    if(alpha_gen==0) {
        costh = 2*X - 1;
    }
    // Dirac case (alpha = 1,-1 based for L/R handed HNLs)
    else {
        costh = -1./alpha_gen + sqrt(1./std::pow(alpha_gen,2) + (4*X - 2)/alpha_gen + 1);
    }

    double theta = std::acos(costh);
    double phi = random->Uniform(0,2*LI::utilities::Constants::pi);
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
    decay.secondary_momenta.resize(2);
    decay.secondary_helicity.resize(1);

    decay.signature.secondary_types[0] = LI::dataclasses::Particle::ParticleType::Gamma;
    decay.secondary_masses[0] = 0;
    decay.secondary_momenta[0][0] = pGamma_lab.e();
    decay.secondary_momenta[0][1] = pGamma_lab.px();
    decay.secondary_momenta[0][2] = pGamma_lab.py();
    decay.secondary_momenta[0][3] = pGamma_lab.pz();
    decay.secondary_momenta[1][0] = pGamma_HNLrest.e();
    decay.secondary_momenta[1][1] = pGamma_HNLrest.px();
    decay.secondary_momenta[1][2] = pGamma_HNLrest.py();
    decay.secondary_momenta[1][3] = pGamma_HNLrest.pz();
    decay.secondary_helicity[0] = std::copysign(1.0, decay.primary_helicity);

    decay.decay_parameters.resize(5);
    decay.decay_parameters[0] = decay_length;
    decay.decay_parameters[1] = decay_weight;
    decay.decay_parameters[2] = (1+alpha_phys*costh)/(1+alpha_gen*costh);
    decay.decay_parameters[3] = a;
    decay.decay_parameters[4] = b;
}

void InjectorBase::SamplePairProduction(LI::dataclasses::DecayRecord const & decay, LI::dataclasses::InteractionRecord & interaction) const {
    // function for simulating the pair production of the photon created in HNL decay
    // considers the different radiation lengths of materials in the detector
    // Nick TODO: comment more

    LI::detector::MaterialModel const & mat_model = earth_model->GetMaterials();
    LI::math::Vector3D decay_vtx(decay.decay_vertex);
    unsigned int gamma_index = 0;
    LI::math::Vector3D decay_dir(decay.secondary_momenta[gamma_index][1],
                                   decay.secondary_momenta[gamma_index][2],
                                   decay.secondary_momenta[gamma_index][3]);

    interaction.signature.primary_type = decay.signature.secondary_types[gamma_index];
    interaction.primary_mass = decay.secondary_masses[gamma_index];
    interaction.primary_momentum = decay.secondary_momenta[gamma_index];
    interaction.primary_helicity = decay.secondary_helicity[gamma_index];

    decay_dir.normalize();

    LI::detector::Path path(earth_model, decay_vtx, decay_dir, 0);
    path.ComputeIntersections();

    std::vector<double> X0;
    std::vector<double> P;
    std::vector<double> D;
    double x0; double p; double density;
    D.push_back(0.);
    double N = 0;
    double lnP_nopp = 0; // for calculating the probability that no pair production occurs
    LI::geometry::Geometry::IntersectionList const & ilist = path.GetIntersections();
    LI::math::Vector3D density_point = decay_vtx;
    unsigned int i = 0;
    for(unsigned int j=0; j < ilist.intersections.size(); ++j) {
        auto const & intersection = ilist.intersections[j];
        if(intersection.distance<0 || std::isinf(intersection.distance))
            continue;
        D.push_back(intersection.distance);
        x0 = (9./7.)*mat_model.GetMaterialRadiationLength(intersection.matID); // in g/cm^2
        density_point += 0.5*(intersection.position - density_point);
        density = earth_model->GetMassDensity(density_point);
        x0 *= 0.01/density; // in m
        X0.push_back(x0);
        p = std::exp(-D[i]/x0) - std::exp(-D[i+1]/x0);
        P.push_back(p);
        N += p;
        lnP_nopp += -(D[i+1] - D[i])/x0;
        density_point = intersection.position;
        ++i;
    }

    interaction.interaction_parameters.resize(1);
    interaction.interaction_parameters[0] = std::exp(lnP_nopp);

    // sample the PDF by inverting the CDF
    double X = random->Uniform(0, 1);
    double C = 0;
    if(P.size() > 0) {
        unsigned int j=0;
        for(; j < P.size(); ++j){
            C += P[j]/N;
            if(C>X) {C -= P[j]/N; break;}
        }
        double pairprod_dist = -X0[j]*std::log(-N*(X - C) + std::exp(-D[j]/X0[j]));
        interaction.interaction_vertex = LI::math::Vector3D(decay.decay_vertex) + pairprod_dist * decay_dir;
    }
}

LI::dataclasses::InteractionRecord InjectorBase::GenerateEvent() {
    LI::dataclasses::InteractionRecord record;
    while(true) {
        try {
            record = this->NewRecord();
            for(auto & distribution : distributions) {
                distribution->Sample(random, earth_model, cross_sections, record);
            }
            SampleCrossSection(record);
            break;
        } catch(LI::utilities::InjectionFailure const & e) {
            continue;
        }
    }
    injected_events += 1;
    return record;
}

double InjectorBase::GenerationProbability(LI::dataclasses::InteractionRecord const & record) const {
    double probability = 1.0;
    for(auto const & dist : distributions) {
        double prob = dist->GenerationProbability(earth_model, cross_sections, record);
        probability *= prob;
    }
    double prob = LI::injection::CrossSectionProbability(earth_model, cross_sections, record);
    probability *= prob;
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
    std::vector<std::shared_ptr<LI::crosssections::CrossSection>> xs_vec = cross_sections->GetCrossSections();
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

std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectorBase::InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return std::pair<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
}

std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> InjectorBase::GetInjectionDistributions() const {
    return distributions;
}

std::shared_ptr<LI::detector::EarthModel> InjectorBase::GetEarthModel() const {
    return earth_model;
}

std::shared_ptr<LI::crosssections::CrossSectionCollection> InjectorBase::GetCrossSections() const {
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

} // namespace injection
} // namespace LI

