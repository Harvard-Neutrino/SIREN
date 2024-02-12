#include "LeptonInjector/injection/Injector.h"

#include <array>
#include <cmath>
#include <string>
#include <algorithm>

#include <rk/rk.hh>

#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/interactions/Decay.h"
#include "LeptonInjector/dataclasses/DecayRecord.h"
#include "LeptonInjector/dataclasses/DecaySignature.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/detector/DetectorModel.h"
#include "LeptonInjector/detector/MaterialModel.h"
#include "LeptonInjector/detector/Path.h"
#include "LeptonInjector/detector/Coordinates.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"
#include "LeptonInjector/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/injection/Process.h"
#include "LeptonInjector/injection/WeightingUtils.h"
#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/utilities/Constants.h"
#include "LeptonInjector/utilities/Errors.h"
#include "LeptonInjector/utilities/Random.h"

namespace LI {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class Injector
//---------------

Injector::Injector() {}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> detector_model,
        std::shared_ptr<LI::utilities::LI_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::InjectionProcess> primary_process,
        std::shared_ptr<LI::utilities::LI_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{
    SetPrimaryProcess(primary_process);
}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::InjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{
    SetPrimaryProcess(primary_process);
    for(auto secondary_process : secondary_processes) {
        AddSecondaryProcess(secondary_process);
    }
}


std::shared_ptr<distributions::VertexPositionDistribution> Injector::FindPrimaryVertexDistribution(std::shared_ptr<LI::injection::InjectionProcess> process) {
    for(auto distribution : process->GetInjectionDistributions()) {
        distributions::VertexPositionDistribution * raw_ptr = dynamic_cast<distributions::VertexPositionDistribution*>(distribution.get());
        if(raw_ptr)
            return std::dynamic_pointer_cast<distributions::VertexPositionDistribution>(distribution);
    }
    throw(LI::utilities::AddProcessFailure("No primary vertex distribution specified!"));
}

std::shared_ptr<distributions::SecondaryVertexPositionDistribution> Injector::FindSecondaryVertexDistribution(std::shared_ptr<LI::injection::SecondaryInjectionProcess> process) {
    for(std::shared_ptr<distributions::SecondaryInjectionDistribution> distribution : process->GetSecondaryInjectionDistributions()) {
        distributions::SecondaryVertexPositionDistribution * raw_ptr = dynamic_cast<distributions::SecondaryVertexPositionDistribution*>(distribution.get());
        if(raw_ptr)
            return std::dynamic_pointer_cast<distributions::SecondaryVertexPositionDistribution>(distribution);
    }
    throw(LI::utilities::AddProcessFailure("No secondary vertex distribution specified!"));
}

void Injector::SetPrimaryProcess(std::shared_ptr<LI::injection::InjectionProcess> primary) {
    std::shared_ptr<distributions::VertexPositionDistribution> vtx_dist;
    try {
        vtx_dist = FindPrimaryVertexDistribution(primary);
    } catch(LI::utilities::AddProcessFailure const & e) {
        std::cerr << e.what() << std::endl;
        exit(0);
    }
    primary_process = primary;
    primary_position_distribution = vtx_dist;
}

void Injector::AddSecondaryProcess(std::shared_ptr<LI::injection::SecondaryInjectionProcess> secondary) {
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> vtx_dist;
    try {
        vtx_dist = FindSecondaryVertexDistribution(secondary);
    } catch(LI::utilities::AddProcessFailure const & e) {
        std::cerr << e.what() << std::endl;
        exit(0);
    }
    secondary_processes.push_back(secondary);
    secondary_position_distributions.push_back(vtx_dist);
    secondary_process_map.insert({secondary->GetPrimaryType(), secondary});
    secondary_position_distribution_map.insert({secondary->GetPrimaryType(), vtx_dist});
}

LI::dataclasses::InteractionRecord Injector::NewRecord() const {
    LI::dataclasses::InteractionRecord record;
    record.signature.primary_type = primary_process->GetPrimaryType();
    return record;
}

void Injector::SetRandom(std::shared_ptr<LI::utilities::LI_random> random) {
    this->random = random;
}

void Injector::SampleCrossSection(LI::dataclasses::InteractionRecord & record) const {
    SampleCrossSection(record, primary_process->GetInteractions());
}

void Injector::SampleCrossSection(LI::dataclasses::InteractionRecord & record, std::shared_ptr<LI::interactions::InteractionCollection> interactions) const {

    // Make sure the particle has interacted
    if(std::isnan(record.interaction_vertex[0]) ||
            std::isnan(record.interaction_vertex[1]) ||
            std::isnan(record.interaction_vertex[2])) {
        throw(LI::utilities::InjectionFailure("No particle interaction!"));
    }

    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = interactions->TargetTypes();

    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::set<LI::dataclasses::Particle::ParticleType> available_targets = detector_model->GetAvailableTargets(intersections, DetectorPosition(record.interaction_vertex));

    double total_prob = 0.0;
    double xsec_prob = 0.0;
    std::vector<double> probs;
    std::vector<LI::dataclasses::Particle::ParticleType> matching_targets;
    std::vector<LI::dataclasses::InteractionSignature> matching_signatures;
    std::vector<std::shared_ptr<LI::interactions::CrossSection>> matching_cross_sections;
    std::vector<std::shared_ptr<LI::interactions::Decay>> matching_decays;
    LI::dataclasses::InteractionRecord fake_record = record;
    double fake_prob;
    if (interactions->HasCrossSections()) {
        for(auto const target : available_targets) {
            if(possible_targets.find(target) != possible_targets.end()) {
                // Get target density
                double target_density = detector_model->GetParticleDensity(intersections, DetectorPosition(interaction_vertex), target);
                // Loop over cross sections that have this target
                std::vector<std::shared_ptr<LI::interactions::CrossSection>> const & target_cross_sections = interactions->GetCrossSectionsForTarget(target);
                for(auto const & cross_section : target_cross_sections) {
                    // Loop over cross section signatures with the same target
                    std::vector<LI::dataclasses::InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                    for(auto const & signature : signatures) {
                        fake_record.signature = signature;
                        fake_record.target_mass = detector_model->GetTargetMass(target);
                        // Add total cross section times density to the total prob
                        fake_prob = target_density * cross_section->TotalCrossSection(fake_record);
                        total_prob += fake_prob;
                        xsec_prob += fake_prob;
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
    }
    if (interactions->HasDecays()) {
        for(auto const & decay : interactions->GetDecays() ) {
            for(auto const & signature : decay->GetPossibleSignaturesFromParent(record.signature.primary_type)) {
                fake_record.signature = signature;
                // fake_prob has units of 1/cm to match cross section probabilities
                fake_prob = 1./(decay->TotalDecayLengthForFinalState(fake_record)/LI::utilities::Constants::cm);
                total_prob += fake_prob;
                // Add total prob to probs
                probs.push_back(total_prob);
                // Add target and decay pointer to the lists
                matching_targets.push_back(LI::dataclasses::Particle::ParticleType::Decay);
                matching_decays.push_back(decay);
                matching_signatures.push_back(signature);
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
    record.target_mass = detector_model->GetTargetMass(record.signature.target_type);
    LI::dataclasses::CrossSectionDistributionRecord xsec_record(record);
    if(r <= xsec_prob)
        matching_cross_sections[index]->SampleFinalState(xsec_record, random);
    else
        matching_decays[index - matching_cross_sections.size()]->SampleFinalState(xsec_record, random);
    xsec_record.Finalize(record);
}

void Injector::SampleNeutrissimoDecay(LI::dataclasses::InteractionRecord const & interaction, LI::dataclasses::DecayRecord & decay, double decay_width, double alpha_gen, double alpha_phys, LI::geometry::Geometry *fiducial = nullptr, double buffer = 0) const {
    // This function takes an interaction record containing an HNL and simulates the decay to a photon
    // Samples according to (1 + alpha * cos(theta))/2 and returns physical weight
    // Final state photon added to secondary particle vectors in LI::dataclasses::InteractionRecord

    // Find the HNL in the secondary particle vector and save its momentum/cartesian direction
    unsigned int lepton_index = (interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4 or interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuF4Bar) ? 0 : 1;
    LI::dataclasses::Particle::ParticleType hnl_type = interaction.signature.secondary_types[lepton_index];
    double hnl_mass = interaction.secondary_masses[lepton_index];
    std::array<double, 4> hnl_momentum = interaction.secondary_momenta[lepton_index];
    double hnl_helicity = interaction.secondary_helicities[lepton_index];

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
    decay.secondary_helicities.resize(1);

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
    decay.secondary_helicities[0] = std::copysign(1.0, decay.primary_helicity);

    decay.decay_parameters.clear();
    decay.decay_parameters["decay_length"] = decay_length;
    decay.decay_parameters["decay_weight"] = decay_weight;
    decay.decay_parameters["angular_weight"] = (1+alpha_phys*costh)/(1+alpha_gen*costh);
    decay.decay_parameters["min_distance"] = a;
    decay.decay_parameters["max_distance"] = b;
}

// Function to sample secondary processes
//
// Modifies an InteractionRecord with the new event
LI::dataclasses::InteractionRecord Injector::SampleSecondaryProcess(LI::dataclasses::SecondaryDistributionRecord & secondary_record) const {
    std::shared_ptr<LI::injection::SecondaryInjectionProcess> secondary_process = secondary_process_map.at(secondary_record.type);
    std::shared_ptr<LI::interactions::InteractionCollection> secondary_interactions = secondary_process->GetInteractions();
    std::vector<std::shared_ptr<LI::distributions::SecondaryInjectionDistribution>> secondary_distributions = secondary_process->GetSecondaryInjectionDistributions();

    size_t max_tries = 1000;
    size_t tries = 0;
    size_t failed_tries = 0;
    std::cout << "Trying to generate secondary process!" << std::endl;
    while(true) {
        try {
            for(auto & distribution : secondary_distributions) {
                distribution->Sample(random, detector_model, secondary_process->GetInteractions(), secondary_record);
            }
            LI::dataclasses::InteractionRecord record;
            secondary_record.Finalize(record);
            SampleCrossSection(record, secondary_interactions);
            return record;
        } catch(LI::utilities::InjectionFailure const & e) {
            failed_tries += 1;
            std::cerr << e.what() << std::endl;
            if(tries > max_tries) {
                std::cerr << "Failed to generate secondary process!" << std::endl;
                std::cout << "Tries: " << tries << " Failed Tries: " << failed_tries << std::endl;
                throw(LI::utilities::InjectionFailure("Failed to generate secondary process!"));
                break;
            }
            continue;
        }
        if(tries > max_tries) {
            std::cerr << "Failed to generate secondary process!" << std::endl;
            std::cout << "Tries: " << tries << " Failed Tries: " << failed_tries << std::endl;
            throw(LI::utilities::InjectionFailure("Failed to generate secondary process!"));
            break;
        }
    }
    return LI::dataclasses::InteractionRecord();
}

LI::dataclasses::InteractionTree Injector::GenerateEvent() {
    LI::dataclasses::InteractionRecord record;
    size_t max_tries = 1000;
    size_t tries = 0;
    size_t failed_tries = 0;
    // Initial Process
    std::cout << "Trying to generate primary process!" << std::endl;
    while(true) {
        tries += 1;
        try {
            std::cout << "Try: " << tries << std::endl;
            record = this->NewRecord();
            for(auto & distribution : primary_process->GetInjectionDistributions()) {
                distribution->Sample(random, detector_model, primary_process->GetInteractions(), record);
            }
            SampleCrossSection(record);
            break;
        } catch(LI::utilities::InjectionFailure const & e) {
            failed_tries += 1;
            std::cerr << e.what() << std::endl;
            if(tries > max_tries) {
                std::cerr << "Failed to generate primary process!" << std::endl;
                std::cout << "Tries: " << tries << " Failed Tries: " << failed_tries << std::endl;
                throw(LI::utilities::InjectionFailure("Failed to generate primary process!"));
                break;
            }
            continue;
        }
        if(tries > max_tries) {
            std::cerr << "Failed to generate primary process!" << std::endl;
            std::cout << "Tries: " << tries << " Failed Tries: " << failed_tries << std::endl;
            throw(LI::utilities::InjectionFailure("Failed to generate primary process!"));
            break;
        }
    }
    std::cout << "Finished primary sampling! Tries: " << tries << " Failed Tries: " << failed_tries << std::endl;
    LI::dataclasses::InteractionTree tree;
    std::shared_ptr<LI::dataclasses::InteractionTreeDatum> parent = tree.add_entry(record);

    // Secondary Processes
    std::deque<std::tuple<std::shared_ptr<LI::dataclasses::InteractionTreeDatum>, std::shared_ptr<LI::dataclasses::SecondaryDistributionRecord>>> secondaries;
    std::function<void(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> add_secondaries = [&](std::shared_ptr<LI::dataclasses::InteractionTreeDatum> parent) {
        for(size_t i=0; parent->record.signature.secondary_types.size(); ++i) {
            LI::dataclasses::ParticleType const & type = parent->record.signature.secondary_types[i];
            std::map<LI::dataclasses::Particle::ParticleType, std::shared_ptr<LI::injection::SecondaryInjectionProcess>>::iterator it = secondary_process_map.find(type);
            if(it == secondary_process_map.end())
                continue;
            if(stopping_condition(parent, i))
                continue;
            secondaries.emplace_back(
                parent,
                std::make_shared<LI::dataclasses::SecondaryDistributionRecord>(parent->record, i)
            );
        }
    };

    add_secondaries(parent);
    std::cout << "secondaries.size() = " << secondaries.size() << std::endl;
    while(secondaries.size() > 0) {
        std::cout << "In while loop" << std::endl;
        std::cout << "Iterating over secondaries" << std::endl;
        for(int i = secondaries.size() - 1; i >= 0; --i) {
            std::shared_ptr<LI::dataclasses::InteractionTreeDatum> parent = std::get<0>(secondaries[i]);
            std::shared_ptr<LI::dataclasses::SecondaryDistributionRecord> secondary_dist = std::get<1>(secondaries[i]);
            secondaries.erase(secondaries.begin() + i);

            std::cout << "Sampling secondary process" << std::endl;
            LI::dataclasses::InteractionRecord secondary_record = SampleSecondaryProcess(*secondary_dist);
            std::cout << "Adding new parent" << std::endl;
            std::shared_ptr<LI::dataclasses::InteractionTreeDatum> secondary_datum = tree.add_entry(secondary_record);
            add_secondaries(secondary_datum);
            std::cout << "Appending new parent to list" << std::endl;
        }
    }
    injected_events += 1;
    std::cout << "Returning tree!" << std::endl;
    return tree;
}

double Injector::SecondaryGenerationProbability(std::shared_ptr<LI::dataclasses::InteractionTreeDatum> const & datum) const {
    return SecondaryGenerationProbability(datum, secondary_process_map.at(datum->record.signature.primary_type));
}

double Injector::SecondaryGenerationProbability(std::shared_ptr<LI::dataclasses::InteractionTreeDatum> const & datum,
        std::shared_ptr<LI::injection::SecondaryInjectionProcess> process) const {
    double probability = 1.0;
    for(auto const & dist : process->GetSecondaryInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), datum->record);
        probability *= prob;
    }
    double prob = LI::injection::CrossSectionProbability(detector_model, process->GetInteractions(), datum->record);
    probability *= prob;
    return probability;
}

double Injector::GenerationProbability(LI::dataclasses::InteractionTree const & tree) const {
    double probability = 1.0;
    std::set<std::shared_ptr<LI::dataclasses::InteractionTreeDatum>>::const_iterator it = tree.tree.cbegin();
    while(it != tree.tree.cend()) {
        if((*it)->depth()==0) probability *= GenerationProbability((*it));
        else probability *= SecondaryGenerationProbability((*it));
        ++it;
    }
    return probability;
}

double Injector::GenerationProbability(std::shared_ptr<LI::dataclasses::InteractionTreeDatum> const & datum,
        std::shared_ptr<LI::injection::InjectionProcess> process) const {
    double probability = 1.0;
    if(!process) { // assume we are dealing with the primary process
        process = primary_process;
        probability *= events_to_inject; // only do this for the primary process
    }
    for(auto const & dist : process->GetInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), datum->record);
        probability *= prob;
    }
    double prob = LI::injection::CrossSectionProbability(detector_model, process->GetInteractions(), datum->record);
    probability *= prob;
    return probability;
}

double Injector::GenerationProbability(LI::dataclasses::InteractionRecord const & record,
        std::shared_ptr<LI::injection::InjectionProcess> process) const {
    double probability = 1.0;
    if(!process) { // assume we are dealing with the primary process
        process = primary_process;
        probability *= events_to_inject; // only do this for the primary process
    }
    for(auto const & dist : process->GetInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), record);
        probability *= prob;
    }
    double prob = LI::injection::CrossSectionProbability(detector_model, process->GetInteractions(), record);
    probability *= prob;
    return probability;
}

// TODO: do we need to save secondary process variables here?
std::set<std::vector<std::string>> Injector::DensityVariables() const {
    std::set<std::vector<std::string>> variable_sets;
    std::vector<std::string> variables;
    for(auto const & dist : primary_process->GetInjectionDistributions()) {
        std::vector<std::string> new_variables = dist->DensityVariables();
        variables.reserve(variables.size() + new_variables.size());
        variables.insert(variables.end(), new_variables.begin(), new_variables.end());
    }
    std::vector<std::shared_ptr<LI::interactions::CrossSection>> xs_vec = primary_process->GetInteractions()->GetCrossSections();
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

std::string Injector::Name() const {
    return("Injector");
}

std::tuple<LI::math::Vector3D, LI::math::Vector3D> Injector::PrimaryInjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    if(!primary_position_distribution) {
        return std::tuple<LI::math::Vector3D, LI::math::Vector3D>(LI::math::Vector3D(0, 0, 0), LI::math::Vector3D(0, 0, 0));
    }
    return primary_position_distribution->InjectionBounds(detector_model, primary_process->GetInteractions(), interaction);
}

// Assumes there is a secondary process and position distribuiton for the provided particle type
std::tuple<LI::math::Vector3D, LI::math::Vector3D> Injector::SecondaryInjectionBounds(LI::dataclasses::InteractionRecord const & record) const {
    return secondary_position_distribution_map.at(record.signature.primary_type)->InjectionBounds(detector_model, secondary_process_map.at(record.signature.primary_type)->GetInteractions(), record);
}

std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> Injector::GetInjectionDistributions() const {
    return primary_process->GetInjectionDistributions();
}

std::shared_ptr<LI::detector::DetectorModel> Injector::GetDetectorModel() const {
    return detector_model;
}

std::shared_ptr<LI::interactions::InteractionCollection> Injector::GetInteractions() const {
    return primary_process->GetInteractions();
}

unsigned int Injector::InjectedEvents() const {
    return injected_events;
}

unsigned int Injector::EventsToInject() const {
    return events_to_inject;
}

Injector::operator bool() const {
    return injected_events < events_to_inject;
}

} // namespace injection
} // namespace LI

