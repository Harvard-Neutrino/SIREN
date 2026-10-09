#include "SIREN/injection/Injector.h"

#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>

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
#include "InteractionSelection.h"
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
    std::shared_ptr<distributions::VertexPositionDistribution> vertex_distribution;
    std::shared_ptr<distributions::VertexPositionDistribution> initial_position_only;
    for(auto distribution : process->GetPrimaryInjectionDistributions()) {
        std::set<distributions::DistributionVariable> variables = distribution->SetVariables();
        std::shared_ptr<distributions::VertexPositionDistribution> position_distribution =
            std::dynamic_pointer_cast<distributions::VertexPositionDistribution>(distribution);
        if(variables.count(distributions::DistributionVariable::InteractionVertex)) {
            vertex_distribution = position_distribution;
        } else if(position_distribution && variables.count(distributions::DistributionVariable::InitialPosition)) {
            // An external distribution carrying only the initial position
            // doubles as a fixed-vertex distribution when nothing downstream
            // places the vertex.
            initial_position_only = position_distribution;
        }
    }
    if(vertex_distribution)
        return vertex_distribution;
    if(initial_position_only)
        return initial_position_only;
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
    std::shared_ptr<distributions::VertexPositionDistribution> vtx_dist = FindPrimaryVertexDistribution(primary);
    primary_process = primary;
    primary_position_distribution = vtx_dist;
}

void Injector::AddSecondaryProcess(std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary) {
    std::shared_ptr<distributions::SecondaryVertexPositionDistribution> vtx_dist = FindSecondaryVertexDistribution(secondary);
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

std::shared_ptr<siren::utilities::SIREN_random> Injector::GetRandom() const {
    return random;
}

void Injector::SampleCrossSection(siren::dataclasses::InteractionRecord & record) const {
    SampleCrossSection(record, primary_process->GetInteractions());
}

void Injector::SampleCrossSection(siren::dataclasses::InteractionRecord & record, std::shared_ptr<siren::interactions::InteractionCollection> interactions) const {
    std::shared_ptr<siren::interactions::Interaction> selected =
        SelectChannel(record, interactions);
    SampleMatchingFinalState(record, selected);
}

std::shared_ptr<siren::interactions::Interaction> Injector::SelectChannel(
    siren::dataclasses::InteractionRecord & record,
    std::shared_ptr<siren::interactions::InteractionCollection> interactions) const
{
    if(!interactions) {
        throw siren::utilities::ConfigurationError(
            "Injector: the process has no InteractionCollection, so no "
            "interaction channel can be selected "
            "[siren-docs: errors#configuration]");
    }
    // Make sure the particle has interacted
    if(std::isnan(record.interaction_vertex[0]) ||
            std::isnan(record.interaction_vertex[1]) ||
            std::isnan(record.interaction_vertex[2])) {
        throw(siren::utilities::InjectionFailure("No particle interaction!"));
    }

    double total_prob = 0.0;
    std::vector<detail::InteractionCandidate> candidates =
        detail::EnumerateInteractionCandidates(
            detector_model, interactions, record);
    std::vector<double> cumulative_rates;
    cumulative_rates.reserve(candidates.size());
    for (detail::InteractionCandidate const & candidate : candidates) {
        total_prob += candidate.rate;
        cumulative_rates.push_back(total_prob);
    }

    if (!std::isfinite(total_prob))
        throw siren::utilities::ConfigurationError("Total interaction rate overflow");
    if(total_prob == 0)
        throw(siren::utilities::InjectionFailure("No valid interactions for this event!"));

    double r = random->Uniform(0, total_prob);
    std::size_t index = 0;
    for(; (index + 1 < cumulative_rates.size())
            && (r > cumulative_rates[index]); ++index) {}
    detail::InteractionCandidate const & selected = candidates[index];
    record.signature = selected.signature;
    record.target_mass = selected.target_mass;
    return selected.interaction;
}

void Injector::SampleMatchingFinalState(
    siren::dataclasses::InteractionRecord & record,
    std::shared_ptr<siren::interactions::Interaction> selected_interaction) const
{
    if(!selected_interaction) {
        throw siren::utilities::ConfigurationError(
            "Injector: no selected interaction is available for final-state "
            "sampling "
            "[siren-docs: errors#configuration]");
    }
    siren::dataclasses::CrossSectionDistributionRecord xsec_record(record);
    double proposed_time = xsec_record.GetInteractionTime();
    if (auto decay = std::dynamic_pointer_cast<siren::interactions::Decay>(
            selected_interaction)) {
        decay->SampleFinalState(xsec_record, random);
        proposed_time = decay->SampleDecayTime(xsec_record, random);
    } else if (auto cross_section =
            std::dynamic_pointer_cast<siren::interactions::CrossSection>(
                selected_interaction)) {
        cross_section->SampleFinalState(xsec_record, random);
        proposed_time = cross_section->SampleInteractionTime(xsec_record, random);
    } else {
        throw siren::utilities::ConfigurationError(
            "Injector: selected interaction is neither a CrossSection nor a "
            "Decay [siren-docs: errors#configuration]");
    }
    if (proposed_time != xsec_record.GetInteractionTime()) {
        xsec_record.SetInteractionTime(proposed_time);
    }
    xsec_record.Finalize(record);
}

namespace {

void PreparePhaseSpaceFinalState(
    siren::dataclasses::InteractionRecord & record,
    std::shared_ptr<siren::interactions::Interaction> selected_interaction)
{
    std::vector<double> secondary_masses;
    std::vector<double> secondary_helicities;
    if (auto decay = std::dynamic_pointer_cast<siren::interactions::Decay>(
            selected_interaction)) {
        secondary_masses = decay->SecondaryMasses(
            record.signature.secondary_types);
        secondary_helicities = decay->SecondaryHelicities(record);
    } else if (auto cross_section =
            std::dynamic_pointer_cast<siren::interactions::CrossSection>(
                selected_interaction)) {
        secondary_masses = cross_section->SecondaryMasses(record);
        secondary_helicities = cross_section->SecondaryHelicities(record);
    } else {
        throw siren::utilities::ConfigurationError(
            "Injector: selected interaction is neither a CrossSection nor a "
            "Decay [siren-docs: errors#configuration]");
    }

    size_t n_secondaries = record.signature.secondary_types.size();
    if (secondary_masses.size() != n_secondaries) {
        throw siren::utilities::ConfigurationError(
            "Injector: SecondaryMasses returned " + std::to_string(secondary_masses.size())
            + " masses for " + std::to_string(n_secondaries) + " secondaries "
            "[siren-docs: errors#configuration]");
    }
    if (secondary_helicities.size() != n_secondaries) {
        secondary_helicities.assign(n_secondaries, 0.0);
    }

    record.secondary_masses = secondary_masses;
    record.secondary_helicities = secondary_helicities;
    record.secondary_ids.resize(n_secondaries);
    for (auto & id : record.secondary_ids) {
        id = siren::dataclasses::ParticleID::GenerateID();
    }
    record.secondary_momenta.assign(n_secondaries, {0.0, 0.0, 0.0, 0.0});
    record.secondary_times.assign(n_secondaries, record.interaction_time);
}

void ApplySelectedInteractionTime(
    siren::dataclasses::InteractionRecord & record,
    std::shared_ptr<siren::interactions::Interaction> selected_interaction,
    std::shared_ptr<siren::utilities::SIREN_random> random)
{
    siren::dataclasses::CrossSectionDistributionRecord xsec_record(record);
    double proposed_time;
    if (auto decay = std::dynamic_pointer_cast<siren::interactions::Decay>(
            selected_interaction)) {
        proposed_time = decay->SampleDecayTime(xsec_record, random);
    } else if (auto cross_section =
            std::dynamic_pointer_cast<siren::interactions::CrossSection>(
                selected_interaction)) {
        proposed_time = cross_section->SampleInteractionTime(xsec_record, random);
    } else {
        throw siren::utilities::ConfigurationError(
            "Injector: selected interaction is neither a CrossSection nor a "
            "Decay [siren-docs: errors#configuration]");
    }
    if (proposed_time != record.interaction_time) {
        record.interaction_time = proposed_time;
    }
    record.secondary_times.assign(
        record.signature.secondary_types.size(), record.interaction_time);
}

} // anonymous namespace

// Function to sample secondary processes
//
// Modifies an InteractionRecord with the new event
siren::dataclasses::InteractionRecord Injector::SampleSecondaryProcess(siren::dataclasses::SecondaryDistributionRecord & secondary_record) const {
    std::shared_ptr<siren::injection::SecondaryInjectionProcess> secondary_process = secondary_process_map.at(secondary_record.type);
    std::shared_ptr<siren::interactions::InteractionCollection> secondary_interactions = secondary_process->GetInteractions();
    std::vector<std::shared_ptr<siren::distributions::SecondaryInjectionDistribution>> secondary_distributions = secondary_process->GetSecondaryInjectionDistributions();

    for(auto & distribution : secondary_distributions) {
        distribution->Sample(random, detector_model, secondary_process->GetInteractions(), secondary_record);
    }
    siren::dataclasses::InteractionRecord record;
    secondary_record.Finalize(record);

    // Select the concrete interaction and signature in one rate-weighted draw.
    std::shared_ptr<siren::interactions::Interaction> selected_interaction =
        SelectChannel(record, secondary_interactions);

    if (secondary_process->HasPhaseSpace(record.signature)) {
        PreparePhaseSpaceFinalState(record, selected_interaction);
        secondary_process->GetPhaseSpace(record.signature)->Sample(random, detector_model, record);
        ApplySelectedInteractionTime(record, selected_interaction, random);
    } else {
        SampleMatchingFinalState(record, selected_interaction);
    }
    return record;
}

siren::dataclasses::InteractionTree Injector::GenerateEvent() {
    if(injection_attempts >= events_to_inject) {
        throw(std::runtime_error("Injector has already made the maximum number of injection attempts!"));
    }
    injection_attempts += 1;
    siren::dataclasses::InteractionRecord record;
    // Initial Process
    try {
        siren::dataclasses::PrimaryDistributionRecord primary_record(primary_process->GetPrimaryType());
        for(auto & distribution : primary_process->GetPrimaryInjectionDistributions()) {
            distribution->Sample(random, detector_model, primary_process->GetInteractions(), primary_record);
        }
        primary_record.Finalize(record);

        // Select the concrete interaction and signature in one rate-weighted draw.
        std::shared_ptr<siren::interactions::Interaction> selected_interaction =
            SelectChannel(record, primary_process->GetInteractions());

        if (primary_process->HasPhaseSpace(record.signature)) {
            PreparePhaseSpaceFinalState(record, selected_interaction);
            primary_process->GetPhaseSpace(record.signature)->Sample(random, detector_model, record);
            ApplySelectedInteractionTime(record, selected_interaction, random);
        } else {
            SampleMatchingFinalState(record, selected_interaction);
        }
    } catch(siren::utilities::InjectionFailure const & e) {
        failed_events += 1;
        int primary_pdg = static_cast<int>(primary_process->GetPrimaryType());
        failure_ledger_.Record(0, primary_pdg, e.reason(), e.what());
        last_failure_message_ = e.what();
        // A vertex-distribution throw fires before Finalize sets the signature,
        // so stamp the primary type to keep the partial tree readable.
        if(record.signature.primary_type == siren::dataclasses::ParticleType::unknown) {
            record.signature.primary_type = primary_process->GetPrimaryType();
        }
        siren::dataclasses::InteractionTree partial_tree;
        partial_tree.add_entry(record);
        last_failed_tree_ = std::move(partial_tree);
        return siren::dataclasses::InteractionTree();
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
                unregistered_secondary_count_ += 1;
                continue;
            }
            if(!stopping_condition || stopping_condition(tree, parent, i)) {
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

            int current_secondary_pdg = static_cast<int>(secondary_dist->type);
            try {
                siren::dataclasses::InteractionRecord secondary_record = SampleSecondaryProcess(*secondary_dist);
                std::shared_ptr<siren::dataclasses::InteractionTreeDatum> secondary_datum = tree.add_entry(secondary_record, parent);
                // Daughter record is authoritative for its production time; keep the
                // parent's secondary_times slot in sync (single write point, after the
                // daughter override is finalized).
                size_t sidx = secondary_dist->GetSecondaryIndex();
                if(sidx < parent->record.secondary_times.size())
                    parent->record.secondary_times[sidx] = secondary_record.primary_initial_time;
                add_secondaries(secondary_datum);
            } catch(siren::utilities::InjectionFailure const & e) {
                failed_events += 1;
                int secondary_depth = static_cast<int>(parent->depth(tree)) + 1;
                int parent_pdg = static_cast<int>(parent->record.signature.primary_type);
                std::ostringstream oss;
                oss << e.what() << " (secondary pdg " << current_secondary_pdg << ")";
                failure_ledger_.Record(secondary_depth, parent_pdg,
                    e.reason(), oss.str());
                last_failure_message_ = e.what();
                last_failed_tree_ = std::move(tree);
                return siren::dataclasses::InteractionTree();
            }
        }
    }
    tree.header.event_number = injected_events;
    tree.header.provenance["generator"] = "SIREN";
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
    double prob;
    auto phase_space = process->GetPhaseSpace(datum->record.signature);
    if(phase_space) {
        prob = siren::injection::CrossSectionProbabilityWithPhaseSpace(
            detector_model, process->GetInteractions(), datum->record,
            *phase_space);
    } else {
        prob = siren::injection::CrossSectionProbability(
            detector_model, process->GetInteractions(), datum->record);
    }
    probability *= prob;
    return probability;
}

double Injector::GenerationProbability(siren::dataclasses::InteractionTree const & tree) const {
    double probability = 1.0;
    std::vector<std::shared_ptr<siren::dataclasses::InteractionTreeDatum>>::const_iterator it = tree.tree.cbegin();
    while(it != tree.tree.cend()) {
        if((*it)->is_root()) probability *= GenerationProbability((*it));
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
    double prob;
    auto phase_space = process->GetPhaseSpace(datum->record.signature);
    if(phase_space) {
        prob = siren::injection::CrossSectionProbabilityWithPhaseSpace(
            detector_model, process->GetInteractions(), datum->record,
            *phase_space);
    } else {
        prob = siren::injection::CrossSectionProbability(
            detector_model, process->GetInteractions(), datum->record);
    }
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
    auto phase_space = process->GetPhaseSpace(record.signature);
    double prob = phase_space
        ? siren::injection::CrossSectionProbabilityWithPhaseSpace(
            detector_model, process->GetInteractions(), record, *phase_space)
        : siren::injection::CrossSectionProbability(
            detector_model, process->GetInteractions(), record);
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

// Assumes there is a secondary process and position distribuiton for the provided particle type
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

unsigned int Injector::InjectionAttempts() const {
    return injection_attempts;
}

unsigned int Injector::EventsToInject() const {
    return events_to_inject;
}

unsigned int Injector::FailedEvents() const {
    return failed_events;
}

unsigned int Injector::UnregisteredSecondaryCount() const {
    return unregistered_secondary_count_;
}

std::string Injector::GetLastFailureMessage() const {
    return last_failure_message_;
}

siren::dataclasses::InteractionTree const & Injector::GetLastFailedTree() const {
    return last_failed_tree_;
}

FailureLedger const & Injector::GetFailureLedger() const {
    return failure_ledger_;
}

void Injector::ResetInjectedEvents(unsigned int events_to_inject) {
    this->events_to_inject = events_to_inject;
    injected_events = 0;
    injection_attempts = 0;
    failed_events = 0;
    unregistered_secondary_count_ = 0;
    failure_ledger_.Clear();
    last_failure_message_.clear();
    last_failed_tree_ = siren::dataclasses::InteractionTree();
}

void Injector::ResetInjectedEvents() {
    injected_events = 0;
    injection_attempts = 0;
    failed_events = 0;
    unregistered_secondary_count_ = 0;
    failure_ledger_.Clear();
    last_failure_message_.clear();
    last_failed_tree_ = siren::dataclasses::InteractionTree();
}

Injector::operator bool() const {
    return events_to_inject == 0 or injected_events < events_to_inject;
}

namespace {
// Header word marking a version-stamped injector archive; headerless
// archives begin directly with the EventsToInject payload.
constexpr std::uint32_t kInjectorArchiveMagic = 0x53494E4A; // "SINJ"
} // anonymous namespace

void Injector::SaveInjector(std::string const & filename) const {
    if (stopping_condition)
        throw std::runtime_error("Injector archives cannot preserve a stopping-condition callback");
    std::ofstream os(filename, std::ios::binary);
    if(!os) {
        throw std::runtime_error(
            "Failed to open injector archive '" + filename + "' for writing");
    }
    ::cereal::BinaryOutputArchive archive(os);
    std::uint32_t magic = kInjectorArchiveMagic;
    // ::cereal::detail::Version<Injector> is what CEREAL_CLASS_VERSION registered,
    // so the header version can never drift from the class version.
    std::uint32_t version = ::cereal::detail::Version<Injector>::version;
    archive(magic, version);
    this->save(archive, version);
}

void Injector::LoadInjector(std::string const & filename) {
    // A missing file otherwise surfaces as a cryptic cereal stream error; name it.
    {
        std::ifstream is(filename, std::ios::binary);
        if(!is) {
            throw std::runtime_error(
                "Failed to load injector archive '" + filename + "': cannot open file");
        }
    }
    auto preserved_random = random;
    auto preserved_stopping_condition = stopping_condition;
    {
        std::ifstream is(filename, std::ios::binary);
        ::cereal::BinaryInputArchive archive(is);
        std::uint32_t magic = 0;
        try {
            archive(magic);
        } catch(...) {
            // Too short to hold a magic word; fall through to the headerless path.
            magic = 0;
        }
        if(magic == kInjectorArchiveMagic) {
            try {
                std::uint32_t version = 0;
                archive(version);
                Injector temp;
                temp.load(archive, version);
                *this = std::move(temp);
            } catch(std::exception const & e) {
                throw std::runtime_error(
                    "Failed to load injector archive '" + filename
                    + "': the headered parse failed: " + e.what());
            }
            // Version 2 restored the RNG into `random`; keep it. Older archives
            // left it null, so fall back to the pre-load engine.
            if(!random) random = preserved_random;
            stopping_condition = preserved_stopping_condition;
            return;
        }
    }
    // Headerless (legacy version-0) archive: the first word was the
    // EventsToInject payload, not a magic word. Reparse the whole stream from
    // the start as the version-0 schema, again into a temporary.
    try {
        std::ifstream is(filename, std::ios::binary);
        ::cereal::BinaryInputArchive archive(is);
        Injector temp;
        temp.load(archive, 0);
        *this = std::move(temp);
    } catch(std::exception const & e) {
        throw std::runtime_error(
            "Failed to load injector archive '" + filename
            + "': not a headered archive and the headerless version-0 parse "
              "failed: " + e.what());
    }
    // Headerless archives are version 0 and never carry the RNG, so restore the
    // pre-load engine.
    if(!random) random = preserved_random;
    stopping_condition = preserved_stopping_condition;
}

} // namespace injection
} // namespace siren
