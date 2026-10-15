#include "SIREN/injection/WeightingUtils.h"
#include "SIREN/injection/PhaseSpaceChannel.h"
#include "SIREN/injection/Process.h"
#include "InteractionSelection.h"

#include <vector>                                                 // for vector
#include <cmath>                                                  // for isfinite
#include <cstddef>                                                // for size_t
#include <string>                                                 // for to_string
#include <sstream>

#include "SIREN/utilities/Errors.h"                     // for WeightCalculationError
#include "SIREN/interactions/CrossSection.h"            // for Cro...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/interactions/Decay.h"                   // for Decay
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/dataclasses/InteractionSignature.h"      // for Int...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...

namespace siren {
namespace injection {

namespace {

double ConvertFinalStateDensity(
    double density,
    siren::dataclasses::PhaseSpaceTopology source_topology,
    siren::dataclasses::PhaseSpaceMeasure const & source_measure,
    PhaseSpaceConvention const & target,
    siren::dataclasses::InteractionRecord const & record)
{
    if (source_topology != target.topology) {
        throw siren::utilities::MeasureCompatibilityError(
            "Final-state topology does not match the requested weighting "
            "convention [siren-docs: errors#measure-compat]");
    }
    if (source_measure == target.measure) return density;
    return ConvertDensity(
        density, source_measure, target.measure, source_topology, record);
}

void MergeConvention(
    PhaseSpaceConvention & result,
    bool & found,
    siren::dataclasses::PhaseSpaceTopology topology,
    siren::dataclasses::PhaseSpaceMeasure const & measure)
{
    if (!found) {
        result.topology = topology;
        result.measure = measure;
        found = true;
        return;
    }
    if (result.topology == topology && result.measure == measure) return;

    if (result.topology == topology &&
        PhaseSpaceDensityConvertible(
            topology, result.measure, measure)) {
        result.measure = measure;
        return;
    }
    if (result.topology == topology &&
        PhaseSpaceDensityConvertible(
            topology, measure, result.measure)) {
        return;
    }

    std::ostringstream oss;
    oss << "Multiple interaction models match one signature but declare "
        << "different final-state conventions ("
        << siren::dataclasses::PhaseSpaceTopologyName(result.topology) << "/"
        << siren::dataclasses::PhaseSpaceMeasureName(result.measure) << " and "
        << siren::dataclasses::PhaseSpaceTopologyName(topology) << "/"
        << siren::dataclasses::PhaseSpaceMeasureName(measure)
        << "); there is no common pointwise density measure "
        << "[siren-docs: errors#measure-compat]";
    throw siren::utilities::MeasureCompatibilityError(oss.str());
}

// The interaction candidates at a vertex, with the total rate and the rate of
// those matching the record's signature. Enumerating candidates traverses the
// geometry, so each probability enumerates them once.
struct CandidateRates {
    std::vector<detail::InteractionCandidate> candidates;
    double total = 0.0;
    double selected = 0.0;
};

CandidateRates AccumulateRates(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    CandidateRates rates;
    rates.candidates = detail::EnumerateInteractionCandidates(
        detector_model, interactions, record);
    for (detail::InteractionCandidate const & candidate : rates.candidates) {
        rates.total += candidate.rate;
        if (candidate.signature == record.signature) {
            rates.selected += candidate.rate;
        }
    }
    return rates;
}

// The rate-weighted FinalStateProbability of the candidates matching the
// record's signature, in `convention`.
double RateWeightedFinalStateProbability(
    std::vector<detail::InteractionCandidate> const & candidates,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention)
{
    double selected_final_state = 0.0;
    for (detail::InteractionCandidate const & candidate : candidates) {
        if (candidate.signature != record.signature) continue;
        double density;
        if (auto decay = std::dynamic_pointer_cast<siren::interactions::Decay>(
                candidate.interaction)) {
            density = ConvertFinalStateDensity(
                decay->FinalStateProbability(record),
                decay->TopologyForSignature(candidate.signature),
                decay->MeasureForSignature(candidate.signature),
                convention, record);
        } else if (auto cross_section =
                std::dynamic_pointer_cast<siren::interactions::CrossSection>(
                    candidate.interaction)) {
            density = ConvertFinalStateDensity(
                cross_section->FinalStateProbability(record),
                cross_section->TopologyForSignature(candidate.signature),
                cross_section->MeasureForSignature(candidate.signature),
                convention, record);
        } else {
            throw siren::utilities::ConfigurationError(
                "Interaction candidate is neither a CrossSection nor a Decay "
                "[siren-docs: errors#configuration]");
        }
        selected_final_state += candidate.rate * density;
    }

    return selected_final_state;
}

} // anonymous namespace

PhaseSpaceConvention SelectedFinalStateConvention(
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    PhaseSpaceConvention result;
    bool found = false;

    for (auto const & decay : interactions->GetDecays()) {
        for (auto const & signature :
             decay->GetPossibleSignaturesFromParent(
                 record.signature.primary_type)) {
            if (signature != record.signature) continue;
            MergeConvention(
                result, found, decay->TopologyForSignature(signature),
                decay->MeasureForSignature(signature));
        }
    }
    for (auto const & target_xs : interactions->GetCrossSectionsByTarget()) {
        for (auto const & cross_section : target_xs.second) {
            for (auto const & signature :
                 cross_section->GetPossibleSignaturesFromParents(
                     record.signature.primary_type, target_xs.first)) {
                if (signature != record.signature) continue;
                MergeConvention(
                    result, found,
                    cross_section->TopologyForSignature(signature),
                    cross_section->MeasureForSignature(signature));
            }
        }
    }
    return result;
}

PhaseSpaceConvention ResolveCommonFinalStateConvention(
    PhaseSpaceConvention const & first,
    PhaseSpaceConvention const & second)
{
    if (first == second) return first;

    if (first.topology != second.topology) {
        std::ostringstream oss;
        oss << "Final-state densities have different topologies ("
            << siren::dataclasses::PhaseSpaceTopologyName(first.topology)
            << " and "
            << siren::dataclasses::PhaseSpaceTopologyName(second.topology)
            << "); there is no common pointwise convention "
            << "[siren-docs: errors#measure-compat]";
        throw siren::utilities::MeasureCompatibilityError(oss.str());
    }

    if (PhaseSpaceDensityConvertible(
            first.topology, first.measure, second.measure)) {
        return second;
    }
    if (PhaseSpaceDensityConvertible(
            first.topology, second.measure, first.measure)) {
        return first;
    }

    std::ostringstream oss;
    oss << "Final-state densities use incompatible measures ("
        << siren::dataclasses::PhaseSpaceMeasureName(first.measure)
        << " and "
        << siren::dataclasses::PhaseSpaceMeasureName(second.measure)
        << ") within "
        << siren::dataclasses::PhaseSpaceTopologyName(first.topology)
        << "; there is no common pointwise convention "
        << "[siren-docs: errors#measure-compat]";
    throw siren::utilities::MeasureCompatibilityError(oss.str());
}

double SelectedFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    return SelectedFinalStateProbability(
        detector_model, interactions, record,
        SelectedFinalStateConvention(interactions, record));
}

double SelectedFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention)
{
    // For one matching model, preserve the historical exact direct-density
    // path. For several models sharing one signature, condition on that
    // signature and return their rate-weighted density mixture.
    std::size_t matching_count = 0;
    double direct_density = 0.0;
    if (interactions->HasDecays()) {
        for (auto const & decay : interactions->GetDecays()) {
            for (auto const & sig :
                 decay->GetPossibleSignaturesFromParent(
                     record.signature.primary_type)) {
                if (sig == record.signature) {
                    ++matching_count;
                    direct_density = ConvertFinalStateDensity(
                        decay->FinalStateProbability(record),
                        decay->TopologyForSignature(sig),
                        decay->MeasureForSignature(sig),
                        convention, record);
                }
            }
        }
    }
    if (interactions->HasCrossSections()) {
        for (auto const & xs_list : interactions->GetCrossSectionsByTarget()) {
            for (auto const & xs : xs_list.second) {
                for (auto const & sig :
                     xs->GetPossibleSignaturesFromParents(
                         record.signature.primary_type,
                         xs_list.first)) {
                    if (sig == record.signature) {
                        ++matching_count;
                        direct_density = ConvertFinalStateDensity(
                            xs->FinalStateProbability(record),
                            xs->TopologyForSignature(sig),
                            xs->MeasureForSignature(sig),
                            convention, record);
                    }
                }
            }
        }
    }
    if (matching_count <= 1) return matching_count == 1 ? direct_density : 0.0;

    CandidateRates rates = AccumulateRates(detector_model, interactions, record);
    if (rates.selected <= 0.0 || !std::isfinite(rates.selected)) {
        throw siren::utilities::WeightCalculationError(
            "SelectedFinalStateProbability: cannot form a rate-conditional "
            "density from selected_rate=" + std::to_string(rates.selected));
    }
    double conditional_density = RateWeightedFinalStateProbability(
        rates.candidates, record, convention) / rates.selected;
    if (!std::isfinite(conditional_density)) {
        throw siren::utilities::WeightCalculationError(
            "SelectedFinalStateProbability: rate-conditional density is "
            "non-finite");
    }
    return conditional_density;
}

double ChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    CandidateRates rates = AccumulateRates(detector_model, interactions, record);
    if (rates.total == 0) return 0.0;
    return rates.selected / rates.total;
}

double FixedVertexChannelSelectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    CandidateRates rates = AccumulateRates(detector_model, interactions, record);
    double total_rate = rates.total;
    double selected_rate = rates.selected;
    std::size_t candidate_count = rates.candidates.size();

    if (candidate_count <= 1) {
        return 1.0;
    }

    // Multiple channels compete: fail loud rather than silently returning a
    // wrong factor.
    if (total_rate <= 0.0) {
        throw siren::utilities::WeightCalculationError(
            "FixedVertexChannelSelectionProbability: total interaction rate is "
            "non-positive (" + std::to_string(total_rate) + ") across "
            + std::to_string(candidate_count)
            + " competing channels; cannot form the channel-selection factor.");
    }
    double probability = selected_rate / total_rate;
    if (not std::isfinite(probability)) {
        throw siren::utilities::WeightCalculationError(
            "FixedVertexChannelSelectionProbability: non-finite channel-selection "
            "probability (selected_rate=" + std::to_string(selected_rate)
            + ", total_rate=" + std::to_string(total_rate) + ").");
    }
    return probability;
}

double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    return CrossSectionProbability(
        detector_model, interactions, record,
        SelectedFinalStateConvention(interactions, record));
}

double CrossSectionProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention)
{
    CandidateRates rates = AccumulateRates(detector_model, interactions, record);
    if (rates.total == 0) return 0.0;
    return RateWeightedFinalStateProbability(rates.candidates, record, convention)
        / rates.total;
}

double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    MultiChannelPhaseSpace const & phase_space)
{
    return CrossSectionProbabilityWithPhaseSpace(
        detector_model, interactions, record, phase_space,
        phase_space.CommonConvention());
}

double CrossSectionProbabilityWithPhaseSpace(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record,
    MultiChannelPhaseSpace const & phase_space,
    PhaseSpaceConvention const & convention)
{
    CandidateRates rates = AccumulateRates(detector_model, interactions, record);
    if (rates.total == 0 || rates.selected == 0) return 0.0;

    double mc_density = phase_space.DensityIn(
        detector_model, record, convention);

    return (rates.selected * mc_density) / rates.total;
}

PhaseSpaceConvention ProcessFinalStateConvention(
    PhysicalProcess const & process,
    siren::dataclasses::InteractionRecord const & record)
{
    auto phase_space = process.GetPhaseSpace(record.signature);
    if (phase_space) return phase_space->CommonConvention();
    if (!process.GetInteractions()) return PhaseSpaceConvention();
    return SelectedFinalStateConvention(process.GetInteractions(), record);
}

double ProcessFinalStateProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    PhysicalProcess const & process,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention)
{
    auto interactions = process.GetInteractions();
    auto phase_space = process.GetPhaseSpace(record.signature);
    if (process.GetWeightingMode().compute_interaction_probability) {
        return phase_space
            ? CrossSectionProbabilityWithPhaseSpace(
                detector_model, interactions, record, *phase_space, convention)
            : CrossSectionProbability(
                detector_model, interactions, record, convention);
    }
    double density = phase_space
        ? phase_space->DensityIn(detector_model, record, convention)
        : SelectedFinalStateProbability(
            detector_model, interactions, record, convention);
    return density * FixedVertexChannelSelectionProbability(
        detector_model, interactions, record);
}

} // namespace injection
} // namespace siren
