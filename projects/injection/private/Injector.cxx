#include "SIREN/injection/Injector.h"

#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream> // Added for debug output

#include <rk/rk.hh>

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/DarkNewsCrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/dataclasses/DecaySignature.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/MaterialModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/distributions/Distributions.h"
#include "SIREN/distributions/primary/vertex/DecayRangeFunction.h"
#include "SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "SIREN/distributions/secondary/vertex/SecondaryVertexPositionDistribution.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/WeightingUtils.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

//---------------
// class Injector
//---------------

Injector::Injector() {}

Injector::Injector(
        unsigned int events_to_inject,
        std::string filename,
        std::shared_ptr<siren::utilities::SIREN_random> random) :
    events_to_inject(events_to_inject),
    random(random)
{
    LoadInjector(filename);
}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<siren::detector::DetectorModel> detector_model,
        std::shared_ptr<siren::utilities::SIREN_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<siren::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::shared_ptr<siren::utilities::SIREN_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{
    SetPrimaryProcess(primary_process);
}

Injector::Injector(
        unsigned int events_to_inject,
        std::shared_ptr<siren::detector::DetectorModel> detector_model,
        std::shared_ptr<injection::PrimaryInjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes,
        std::shared_ptr<siren::utilities::SIREN_random> random) :
    events_to_inject(events_to_inject),
    random(random),
    detector_model(detector_model)
{
    SetPrimaryProcess(primary_process);
    for(auto secondary_process : secondary_processes) {
        AddSecondaryProcess(secondary_process);
    }
}

std::shared_ptr<distributions::VertexPositionDistribution> Injector::FindPrimaryVertexDistribution(std::shared_ptr<siren::injection::PrimaryInjectionProcess> process) {
    for(auto distribution : process->GetPrimaryInjectionDistributions()) {
        distributions::VertexPositionDistribution * raw_ptr = dynamic_cast<distributions::VertexPositionDistribution*>(distribution.get());
        if(raw_ptr)
            return std::dynamic_pointer_cast<distributions::VertexPositionDistribution>(distribution);
    }
    throw(siren::utilities::AddProcessFailure("No primary vertex distribution specified!"));
}

std::shared_ptr<distributions::SecondaryVertexPositionDistribution> Injector::FindSecondaryVertexDistribution(std::shared_ptr<siren::injection::SecondaryInjectionProcess> process) {
    for(std::shared_ptr<distributions::SecondaryInjectionDistribution> distribution : process->GetSecondaryInjectionDistributions()) {
        distributions::SecondaryVertexPositionDistribution * raw_ptr = dynamic_cast<distributions::SecondaryVertexPositionDistribution*>(distribution.get());
        if(raw_ptr)
            return std::dynamic_pointer_cast<distributions::SecondaryVertexPositionDistribution>(distribution);
    }
    throw(siren::utilities::AddProcessFailure("No secondary vertex distribution specified!"));
}

void Injector::SetPrimaryProcess(std::shared_ptr<siren::injection::PrimaryInjectionProcess> primary) {
    std::shared_ptr<distributions::VertexPositionDistribution> vtx_dist;
    try {
        vtx_dist = FindPrimaryVertexDistribution(primary);
    } catch(siren::utilities::AddProcessFailure const & e) {
        std::cerr << e.what() << std::endl;
        exit(0);
    }
    primary_process = primary;
    primary_position_distribution = vtx_dist;
}

void Injector::AddSecondaryProcess(std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary) {
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> vtx_dist;
    try {
        vtx_dist = FindSecondaryVertexDistribution(secondary);
    } catch(siren::utilities::AddProcessFailure const & e) {
        std::cerr << e.what() << std::endl;
        exit(0);
    }
    secondary_processes.push_back(secondary);
    secondary_position_distributions.push_back(vtx_dist);
    secondary_process_map.insert({secondary->GetPrimaryType(), secondary});
    secondary_position_distribution_map.insert({secondary->GetPrimaryType(), vtx_dist});
}

void Injector::SetSecondaryProcesses(std::vector<std::shared_ptr<siren::injection::SecondaryInjectionProcess>> secondaries) {
    secondary_processes.clear();
    secondary_position_distributions.clear();
    secondary_process_map.clear();
    secondary_position_distribution_map.clear();
    for(auto secondary : secondaries) {
        AddSecondaryProcess(secondary);
    }
}

siren::dataclasses::InteractionRecord Injector::NewRecord() const {
    siren::dataclasses::InteractionRecord record;
    record.signature.primary_type = primary_process->GetPrimaryType();
    return record;
}

void Injector::SetRandom(std::shared_ptr<siren::utilities::SIREN_random> random) {
    this->random = random;
}

void Injector::SampleCrossSection(siren::dataclasses::InteractionRecord & record) const {
    SampleCrossSection(record, primary_process->GetInteractions());
}

void Injector::SampleCrossSection(siren::dataclasses::InteractionRecord & record, std::shared_ptr<siren::interactions::InteractionCollection> interactions) const {
    // Make sure the particle has interacted
    if(std::isnan(record.interaction_vertex[0]) ||
       std::isnan(record.interaction_vertex[1]) ||
       std::isnan(record.interaction_vertex[2])) {
        std::cout << "InjectionFailure: Interaction vertex contains NaN: [" 
                  << record.interaction_vertex[0] << ", " 
                  << record.interaction_vertex[1] << ", " 
                  << record.interaction_vertex[2] << "]" << std::endl;
        throw(siren::utilities::InjectionFailure("No particle interaction!"));
    }

    std::set<siren::dataclasses::ParticleType> const & possible_targets = interactions->TargetTypes();
    std::cout << "Possible targets: ";
    for (auto const & target : possible_targets) {
        std::cout << target << " ";
    }
    std::cout << std::endl;

    siren::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);
    std::cout << "Interaction vertex: [" 
              << interaction_vertex.GetX() << ", " 
              << interaction_vertex.GetY() << ", " 
              << interaction_vertex.GetZ() << "]" << std::endl;

    siren::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();
    std::cout << "Primary direction: [" 
              << primary_direction.GetX() << ", " 
              << primary_direction.GetY() << ", " 
              << primary_direction.GetZ() << "]" << std::endl;

    siren::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::cout << "Number of intersections: " << intersections.intersections.size() << std::endl;

    std::set<siren::dataclasses::ParticleType> available_targets = detector_model->GetAvailableTargets(intersections, DetectorPosition(record.interaction_vertex));
    std::cout << "Available targets: ";
    for (auto const & target : available_targets) {
        std::cout << target << " ";
    }
    std::cout << std::endl;

    double total_prob = 0.0;
    double xsec_prob = 0.0;
    std::vector<double> probs;
    std::vector<siren::dataclasses::ParticleType> matching_targets;
    std::vector<siren::dataclasses::InteractionSignature> matching_signatures;
    std::vector<std::shared_ptr<siren::interactions::CrossSection>> matching_cross_sections;
    std::vector<std::shared_ptr<siren::interactions::Decay>> matching_decays;
    siren::dataclasses::InteractionRecord fake_record = record;
    double fake_prob;
    if (interactions->HasCrossSections()) {
        for(auto const target : available_targets) {
            if(possible_targets.find(target) != possible_targets.end()) {
                // Get target density
                double target_density = detector_model->GetParticleDensity(intersections, DetectorPosition(interaction_vertex), target);
                std::cout << "Target: " << target << ", Density: " << target_density << std::endl;
                // Loop over cross sections that have this target
                std::vector<std::shared_ptr<siren::interactions::CrossSection>> const & target_cross_sections = interactions->GetCrossSectionsForTarget(target);
                for(auto const & cross_section : target_cross_sections) {
                    // Loop over cross section signatures with the same target
                    std::vector<siren::dataclasses::InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                    for(auto const & signature : signatures) {
                        fake_record.signature = signature;
                        fake_record.target_mass = detector_model->GetTargetMass(target);
                        // Add total cross section times density to the total prob
                        fake_prob = target_density * cross_section->TotalCrossSection(fake_record);
                        std::cout << "Signature: " << signature << ", Cross-section: " << cross_section->TotalCrossSection(fake_record) 
                                  << ", Density * Cross-section: " << fake_prob << std::endl;
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
                fake_prob = 1./(decay->TotalDecayLengthForFinalState(fake_record)/siren::utilities::Constants::cm);
                std::cout << "Decay signature: " << signature << ", Probability: " << fake_prob << std::endl;
                total_prob += fake_prob;
                // Add total prob to probs
                probs.push_back(total_prob);
                // Add target and decay pointer to the lists
                matching_targets.push_back(siren::dataclasses::ParticleType::Decay);
                matching_decays.push_back(decay);
                matching_signatures.push_back(signature);
            }
        }
    }

    std::cout << "Total probability: " << total_prob << std::endl;
    if(total_prob == 0) {
        std::cout << "InjectionFailure: No valid interactions for this event!" << std::endl;
        throw(siren::utilities::InjectionFailure("No valid interactions for this event!"));
    }
    // Throw a random number
    double r = random->Uniform(0, total_prob);
    std::cout << "Random number: " << r << std::endl;
    // Choose the target and cross section
    unsigned int index = 0;
    for(; (index+1 < probs.size()) and (r > probs[index]); ++index) {}
    record.signature.target_type = matching_targets[index];
    record.signature = matching_signatures[index];
    double selected_prob = 0.0;
    for(unsigned int i=0; i<probs.size(); ++i) {
        if(matching_signatures[index] == matching_signatures[i]) {
            selected_prob += (i > 0 ? probs[i] - probs[i - 1] : probs[i]);
        }
    }
    std::cout << "Selected target: " << record.signature.target_type 
              << ", Selected signature: " << record.signature 
              << ", Selected probability: " << selected_prob << std::endl;
    if(selected_prob == 0) {
        std::cout << "InjectionFailure: No valid interactions for this event!" << std::endl;
        throw(siren::utilities::InjectionFailure("No valid interactions for this event!"));
    }
    record.target_mass = detector_model->GetTargetMass(record.signature.target_type);
    siren::dataclasses::CrossSectionDistributionRecord xsec_record(record);
    if(r <= xsec_prob) {
        std::cout << "Sampling cross-section final state" << std::endl;
        matching_cross_sections[index]->SampleFinalState(xsec_record, random);
    } else {
        std::cout << "Sampling decay final state" << std::endl;
        matching_decays[index - matching_cross_sections.size()]->SampleFinalState(xsec_record, random);
    }
    xsec_record.Finalize(record);
}

// Function to sample secondary processes
//
// Modifies an InteractionRecord with the new event
siren::dataclasses::InteractionRecord Injector::SampleSecondaryProcess(siren::dataclasses::SecondaryDistributionRecord & secondary_record) const {
    std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary_process = secondary_process_map.at(secondary_record.type);
    std::shared_ptr<siren::interactions::InteractionCollection> secondary_interactions = secondary_process->GetInteractions();
    std::vector<std::shared_ptr<siren::distributions::SecondaryInjectionDistribution>> secondary_distributions = secondary_process->GetSecondaryInjectionDistributions();

    size_t max_tries = 100;
    size_t tries = 0;
    size_t failed_tries = 0;
    while(true) {
        try {
            for(auto & distribution : secondary_distributions) {
                distribution->Sample(random, detector_model, secondary_process->GetInteractions(), secondary_record);
            }
            siren::dataclasses::InteractionRecord record;
            secondary_record.Finalize(record);
            SampleCrossSection(record, secondary_interactions);
            return record;
        } catch(siren::utilities::InjectionFailure const & e) {
            failed_tries += 1;
            if(tries > max_tries) {
                std::cout << "InjectionFailure: Failed to generate secondary process after " << max_tries << " tries!" << std::endl;
                throw(siren::utilities::InjectionFailure("Failed to generate secondary process!"));
                break;
            }
            continue;
        }
        if(tries > max_tries) {
            std::cout << "InjectionFailure: Failed to generate secondary process after " << max_tries << " tries!" << std::endl;
            throw(siren::utilities::InjectionFailure("Failed to generate secondary process!"));
            break;
        }
    }
    return siren::dataclasses::InteractionRecord();
}

siren::dataclasses::InteractionTree Injector::GenerateEvent() {
    siren::dataclasses::InteractionRecord record;
    size_t max_tries = 100;
    size_t tries = 0;
    size_t failed_tries = 0;
    // Initial Process
    while(true) {
        tries += 1;
        std::cout << "GenerateEvent attempt: " << tries << std::endl;
        try {
            siren::dataclasses::PrimaryDistributionRecord primary_record(primary_process->GetPrimaryType());
            std::cout << "Sampling primary distributions for type: " << primary_process->GetPrimaryType() << std::endl;
            for(auto & distribution : primary_process->GetPrimaryInjectionDistributions()) {
                std::cout << "Sampling distribution: " << distribution->Name() << std::endl;
                distribution->Sample(random, detector_model, primary_process->GetInteractions(), primary_record);
            }
            std::cout << "Sampled energy: " << primary_record.GetEnergy() 
                      << ", position: [" << primary_record.GetInitialPosition()[0] 
                      << ", " << primary_record.GetInitialPosition()[1] 
                      << ", " << primary_record.GetInitialPosition()[2] 
                      << "], direction: [" << primary_record.GetDirection()[0] 
                      << ", " << primary_record.GetDirection()[1] 
                      << ", " << primary_record.GetDirection()[2] << "]" << std::endl;
            primary_record.Finalize(record);
            SampleCrossSection(record);
            break;
        } catch(siren::utilities::InjectionFailure const & e) {
            std::cout << "InjectionFailure: " << e.what() << ", Try: " << tries << std::endl;
            failed_tries += 1;
            if(tries > max_tries) {
                std::cout << "InjectionFailure: Failed to generate primary process after " << max_tries << " tries!" << std::endl;
                throw(siren::utilities::InjectionFailure("Failed to generate primary process!"));
                break;
            }
            continue;
        }
        if(tries > max_tries) {
            std::cout << "InjectionFailure: Failed to generate primary process after " << max_tries << " tries!" << std::endl;
            throw(siren::utilities::InjectionFailure("Failed to generate primary process!"));
            break;
        }
    }
    siren::dataclasses::InteractionTree tree;
    std::shared_ptr<siren::dataclasses::InteractionTreeDatum> parent = tree.add_entry(record);

    // Secondary Processes
    std::deque<std::tuple<std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, std::shared_ptr<siren::dataclasses::SecondaryDistributionRecord>>> secondaries;
    std::function<void(std::shared_ptr<siren::dataclasses::InteractionTreeDatum>)> add_secondaries = [&](std::shared_ptr<siren::dataclasses::InteractionTreeDatum> parent) {
        for(size_t i=0; i<parent->record.signature.secondary_types.size(); ++i) {
            siren::dataclasses::ParticleType const & type = parent->record.signature.secondary_types[i];
            std::map<siren::dataclasses::ParticleType, std::shared_ptr<siren::injection::SecondaryInjectionProcess>>::iterator it = secondary_process_map.find(type);
            if(it == secondary_process_map.end()) {
                continue;
            }
            if(stopping_condition(parent, i)) {
                continue;
            }
            secondaries.emplace_back(
                parent,
                std::make_shared<siren::dataclasses::SecondaryDistributionRecord>(parent->record, i)
            );
        }
    };

    add_secondaries(parent);
    while(secondaries.size() > 0) {
        for(int i = secondaries.size() - 1; i >= 0; --i) {
            std::shared_ptr<siren::dataclasses::InteractionTreeDatum> parent = std::get<0>(secondaries[i]);
            std::shared_ptr<siren::dataclasses::SecondaryDistributionRecord> secondary_dist = std::get<1>(secondaries[i]);
            secondaries.erase(secondaries.begin() + i);

            siren::dataclasses::InteractionRecord secondary_record = SampleSecondaryProcess(*secondary_dist);
            std::shared_ptr<siren::dataclasses::InteractionTreeDatum> secondary_datum = tree.add_entry(secondary_record, parent);
            add_secondaries(secondary_datum);
        }
    }
    injected_events += 1;
    return tree;
}

double Injector::SecondaryGenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum) const {
    return SecondaryGenerationProbability(datum, secondary_process_map.at(datum->record.signature.primary_type));
}

double Injector::SecondaryGenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum,
        std::shared_ptr<siren::injection::SecondaryInjectionProcess> process) const {
    double probability = 1.0;
    for(auto const & dist : process->GetSecondaryInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), datum->record);
        probability *= prob;
    }
    double prob = siren::injection::CrossSectionProbability(detector_model, process->GetInteractions(), datum->record);
    probability *= prob;
    return probability;
}

double Injector::GenerationProbability(siren::dataclasses::InteractionTree const & tree) const {
    double probability = 1.0;
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTreeDatum>>::const_iterator it = tree.tree.cbegin();
    while(it != tree.tree.cend()) {
        if((*it)->depth()==0) probability *= GenerationProbability((*it));
        else probability *= SecondaryGenerationProbability((*it));
        ++it;
    }
    return probability;
}

double Injector::GenerationProbability(std::shared_ptr<siren::dataclasses::InteractionTreeDatum> const & datum,
        std::shared_ptr<siren::injection::PrimaryInjectionProcess> process) const {
    double probability = 1.0;
    if(!process) { // assume we are dealing with the primary process
        process = primary_process;
        unsigned int stat_weight = (events_to_inject > 0) ? events_to_inject : 1;
        probability *= stat_weight; // only do this for the primary process
    }
    for(auto const & dist : process->GetPrimaryInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), datum->record);
        probability *= prob;
    }
    double prob = siren::injection::CrossSectionProbability(detector_model, process->GetInteractions(), datum->record);
    probability *= prob;
    return probability;
}

double Injector::GenerationProbability(siren::dataclasses::InteractionRecord const & record,
        std::shared_ptr<siren::injection::PrimaryInjectionProcess> process) const {
    double probability = 1.0;
    if(!process) { // assume we are dealing with the primary process
        process = primary_process;
        unsigned int stat_weight = (events_to_inject > 0) ? events_to_inject : 1;
        probability *= stat_weight; // only do this for the primary process
    }
    for(auto const & dist : process->GetPrimaryInjectionDistributions()) {
        double prob = dist->GenerationProbability(detector_model, process->GetInteractions(), record);
        probability *= prob;
    }
    double prob = siren::injection::CrossSectionProbability(detector_model, process->GetInteractions(), record);
    probability *= prob;
    return probability;
}

// TODO: do we need to save secondary process variables here?
std::set<std::vector<std::string>> Injector::DensityVariables() const {
    std::set<std::vector<std::string>> variable_sets;
    std::vector<std::string> variables;
    for(auto const & dist : primary_process->GetPrimaryInjectionDistributions()) {
        std::vector<std::string> new_variables = dist->DensityVariables();
        variables.reserve(variables.size() + new_variables.size());
        variables.insert(variables.end(), new_variables.begin(), new_variables.end());
    }
    std::vector<std::shared_ptr<siren::interactions::CrossSection>> xs_vec = primary_process->GetInteractions()->GetCrossSections();
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

std::tuple<siren::math::Vector3D, siren::math::Vector3D> Injector::PrimaryInjectionBounds(siren::dataclasses::InteractionRecord const & interaction) const {
    if(!primary_position_distribution) {
        return std::tuple<siren::math::Vector3D, siren::math::Vector3D>(siren::math::Vector3D(0, 0, 0), siren::math::Vector3D(0, 0, 0));
    }
    return primary_position_distribution->InjectionBounds(detector_model, primary_process->GetInteractions(), interaction);
}

// Assumes there is a secondary process and position distribution for the provided particle type
std::tuple<siren::math::Vector3D, siren::math::Vector3D> Injector::SecondaryInjectionBounds(siren::dataclasses::InteractionRecord const & record) const {
    return secondary_position_distribution_map.at(record.signature.primary_type)->InjectionBounds(detector_model, secondary_process_map.at(record.signature.primary_type)->GetInteractions(), record);
}

std::vector<std::shared_ptr<siren::distributions::PrimaryInjectionDistribution>> Injector::GetPrimaryInjectionDistributions() const {
    return primary_process->GetPrimaryInjectionDistributions();
}

std::shared_ptr<siren::detector::DetectorModel> Injector::GetDetectorModel() const {
    return detector_model;
}

void Injector::SetDetectorModel(std::shared_ptr<siren::detector::DetectorModel> detector_model) {
    this->detector_model = detector_model;
}

std::shared_ptr<siren::interactions::InteractionCollection> Injector::GetInteractions() const {
    return primary_process->GetInteractions();
}

unsigned int Injector::InjectedEvents() const {
    return injected_events;
}

unsigned int Injector::EventsToInject() const {
    return events_to_inject;
}

void Injector::ResetInjectedEvents() {
    injected_events = 0;
}

Injector::operator bool() const {
    return events_to_inject == 0 or injected_events < events_to_inject;
}

void Injector::SaveInjector(std::string const & filename) const {
    std::ofstream os(filename, std::ios::binary);
    ::cereal::BinaryOutputArchive archive(os);
    this->save(archive,0);
}

void Injector::LoadInjector(std::string const & filename) {
    std::ifstream is(filename, std::ios::binary);
    ::cereal::BinaryInputArchive archive(is);
    this->load(archive,0);
}

} // namespace injection
} // namespace siren
