#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Errors.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <typeinfo>

namespace siren {
namespace injection {

using MType = PhaseSpaceMeasure::Type;

namespace {

int MeasurePriority(PhaseSpaceMeasure const & m) {
    switch (m.type) {
        case MType::SolidAngleRest:  return 0;
        case MType::Recursive2Body:  return 0;
        case MType::MandelstamQ2Phi:  return 1;
        case MType::FixedMassYPhi:    return 2;
        case MType::MandelstamQ2YPhi: return 3;
        case MType::BjorkenXYPhi:     return 4;
        case MType::MandelstamQ2:     return 5;
        case MType::FixedMassY:       return 6;
        case MType::MandelstamQ2Y:    return 7;
        case MType::BjorkenXY:        return 8;
        case MType::HelicityAngles:   return 9;
        case MType::DalitzPair:       return 10;
        case MType::SolidAngleLab:    return 11;
        case MType::CosThetaRest:    return 12;
        case MType::Unspecified:     return 13;
    }
    return 12;
}

// Checks that there are channels, one weight per channel, and that every
// weight is finite and nonnegative. Returns the weight sum.
double CheckedWeightSum(
    std::vector<std::shared_ptr<PhaseSpaceChannel>> const & channels,
    std::vector<double> const & weights,
    std::string const & where)
{
    if (channels.empty()) {
        throw siren::utilities::ConfigurationError(
            where + ": the mixture has no channels [siren-docs: errors#configuration]");
    }
    if (channels.size() != weights.size()) {
        throw siren::utilities::ConfigurationError(
            where + ": " + std::to_string(weights.size()) + " weights for "
            + std::to_string(channels.size()) + " channels [siren-docs: errors#configuration]");
    }
    for (size_t i = 0; i < weights.size(); ++i) {
        if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
            std::ostringstream oss;
            oss << where << ": weight " << i << " is " << weights[i]
                << " (must be finite and nonnegative) [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
    }
    return std::accumulate(weights.begin(), weights.end(), 0.0);
}

} // anonymous namespace

MultiChannelPhaseSpace::MultiChannelPhaseSpace(
    std::vector<std::shared_ptr<PhaseSpaceChannel>> channels,
    std::vector<double> weights)
    : channels_(std::move(channels))
{
    SetWeights(std::move(weights));
    RequireCompatibleConvention();
}

void MultiChannelPhaseSpace::SetWeights(std::vector<double> weights) {
    if (weights.empty() && !channels_.empty()) {
        weights.assign(channels_.size(), 1.0 / static_cast<double>(channels_.size()));
    } else {
        double sum = CheckedWeightSum(channels_, weights, "MultiChannelPhaseSpace::SetWeights");
        if (!std::isfinite(sum) || sum <= 0.0) {
            std::ostringstream oss;
            oss << "MultiChannelPhaseSpace::SetWeights: weight sum is " << sum
                << " (must be finite and positive) [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
        for (double & weight : weights) weight /= sum;
    }
    weights_ = std::move(weights);
}

void MultiChannelPhaseSpace::ValidateArchivedWeights() const {
    double sum = CheckedWeightSum(channels_, weights_, "MultiChannelPhaseSpace::load");
    if (!(std::abs(sum - 1.0) <= 1e-9)) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace::load: weights sum to " << sum
            << ", expected 1 [siren-docs: errors#configuration]";
        throw siren::utilities::ConfigurationError(oss.str());
    }
}

MultiChannelPhaseSpace::ConventionEvaluation
MultiChannelPhaseSpace::EvaluateConventions() const {
    if (channels_.empty()) {
        throw siren::utilities::ConfigurationError("MultiChannelPhaseSpace has no channels");
    }
    // Nested channels can be reassigned through shared ownership. Reject
    // recursive mixtures before invoking any sampler or density.
    thread_local std::vector<MultiChannelPhaseSpace const *> active;
    if (std::find(active.begin(), active.end(), this) != active.end()) {
        throw siren::utilities::ConfigurationError("Recursive phase-space mixture");
    }
    active.push_back(this);
    struct PopActive {
        std::vector<MultiChannelPhaseSpace const *> & stack;
        ~PopActive() { stack.pop_back(); }
    } pop{active};

    // Read the channel conventions once. Their backing models may change, and
    // sampling/density use the conventions observed for the current operation.
    ConventionEvaluation result;
    result.measures.reserve(channels_.size());
    for (std::size_t i = 0; i < channels_.size(); ++i) {
        auto const * channel = channels_[i].get();
        if (!channel) {
            throw siren::utilities::ConfigurationError("Null phase-space channel");
        }
        PhaseSpaceConvention convention;
        if (typeid(*channel) == typeid(NestedMixtureChannel)) {
            auto const & nested = static_cast<NestedMixtureChannel const &>(*channel);
            if (!nested.mixture) {
                throw siren::utilities::ConfigurationError(
                    "Nested phase-space channel has no mixture");
            }
            convention = nested.mixture->CommonConvention();
        } else {
            convention = {channel->Topology(), channel->Measure()};
        }
        result.measures.push_back(convention.measure);
        if (i == 0) {
            result.common = convention;
        } else if (convention.topology != result.common.topology) {
            std::ostringstream message;
            message << "Topology mismatch: channel 0 (" << channels_[0]->Name()
                    << ") is " << PhaseSpaceTopologyName(result.common.topology)
                    << " but channel " << i << " (" << channel->Name()
                    << ") is " << PhaseSpaceTopologyName(convention.topology)
                    << " [siren-docs: errors#measure-compat]";
            result.diagnostics.push_back(
                {ChannelDiagnostic::Severity::Fatal, message.str()});
        }
    }
    // A measure comparison is meaningful only when topologies agree.
    if (!result.diagnostics.empty()) return result;

    auto const & measures = result.measures;
    if (std::all_of(measures.begin(), measures.end(),
                    [&](auto const & measure) { return measure == measures.front(); })) {
        return result;
    }

    bool found_viable = false;
    int best_count = -1;
    int best_priority = 100;
    bool has_specified = std::any_of(
        measures.begin(), measures.end(),
        [](PhaseSpaceMeasure const & measure) {
            return measure.type != MType::Unspecified;
        });
    for (auto const & candidate : measures) {
        if (has_specified && candidate.type == MType::Unspecified) continue;
        bool viable = std::all_of(
            measures.begin(), measures.end(),
            [&](PhaseSpaceMeasure const & source) {
                return PhaseSpaceDensityConvertible(
                    result.common.topology, source, candidate);
            });
        if (found_viable && !viable) continue;
        if (!found_viable && viable) {
            found_viable = true;
            best_count = -1;
            best_priority = 100;
        }
        int count = static_cast<int>(std::count(
            measures.begin(), measures.end(), candidate));
        int priority = MeasurePriority(candidate);
        if (count > best_count ||
            (count == best_count && priority < best_priority)) {
            result.common.measure = candidate;
            best_count = count;
            best_priority = priority;
        }
    }

    for (std::size_t i = 0; i < channels_.size(); ++i) {
        if (measures[i] == result.common.measure) continue;
        std::ostringstream diagnostic;
        ChannelDiagnostic::Severity severity;
        if (!PhaseSpaceDensityConvertible(
                result.common.topology, measures[i], result.common.measure)) {
            diagnostic << "Measure incompatibility: channel " << i
                       << " (" << channels_[i]->Name() << ") uses "
                       << PhaseSpaceMeasureName(measures[i])
                       << " which is not convertible to "
                       << PhaseSpaceMeasureName(result.common.measure)
                       << " within "
                       << PhaseSpaceTopologyName(result.common.topology)
                       << " topology [siren-docs: errors#measure-compat]";
            severity = ChannelDiagnostic::Severity::Fatal;
        } else {
            diagnostic << "Channel " << i << " (" << channels_[i]->Name()
                       << ") uses " << PhaseSpaceMeasureName(measures[i])
                       << "; will auto-convert to "
                       << PhaseSpaceMeasureName(result.common.measure);
            severity = ChannelDiagnostic::Severity::Info;
        }
        result.diagnostics.push_back({severity, diagnostic.str()});
    }
    return result;
}

MultiChannelPhaseSpace::ConventionEvaluation
MultiChannelPhaseSpace::RequireCompatibleConvention() const {
    auto result = EvaluateConventions();
    std::string errors;
    for (auto const & diagnostic : result.diagnostics) {
        if (diagnostic.severity == ChannelDiagnostic::Severity::Fatal) {
            errors += "\n  - " + diagnostic.message;
        }
    }
    if (!errors.empty()) {
        throw siren::utilities::MeasureCompatibilityError(
            "MultiChannelPhaseSpace: incompatible channels:" + errors);
    }
    return result;
}

int MultiChannelPhaseSpace::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    RequireCompatibleConvention();

    double r = random->Uniform(0, 1);
    double cumulative = 0.0;
    int selected = static_cast<int>(channels_.size()) - 1;
    for (int i = 0; i < static_cast<int>(channels_.size()); ++i) {
        cumulative += weights_[i];
        if (r < cumulative) {
            selected = i;
            break;
        }
    }

    channels_[selected]->Sample(random, detector_model, record);
    return selected;
}

double MultiChannelPhaseSpace::ComputeDensity(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    ConventionEvaluation const & conventions,
    std::vector<double> * weighted,
    std::vector<double> * bare) const
{
    if (weighted) weighted->reserve(channels_.size());
    if (bare) {
        bare->clear();
        bare->reserve(channels_.size());
    }
    double density = 0.0;
    for (std::size_t i = 0; i < channels_.size(); ++i) {
        double d = channels_[i]->Density(detector_model, record);
        if (!std::isfinite(d) || d < 0.0) {
            throw siren::utilities::WeightCalculationError(
                "Phase-space channel density must be finite and nonnegative");
        }
        d = ConvertDensity(
            d, conventions.measures[i], conventions.common.measure,
            conventions.common.topology, record);
        if (bare) bare->push_back(d);
        double contribution = weights_[i] * d;
        if (weighted) weighted->push_back(contribution);
        density += contribution;
    }
    if (!std::isfinite(density) || density < 0.0) {
        throw siren::utilities::WeightCalculationError("Invalid mixture density");
    }
    return density;
}

double MultiChannelPhaseSpace::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    return ComputeDensity(detector_model, record, RequireCompatibleConvention());
}

double MultiChannelPhaseSpace::DensityIn(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceMeasure const & measure) const
{
    auto conventions = RequireCompatibleConvention();
    double density = ComputeDensity(detector_model, record, conventions);
    return ConvertDensity(
        density, conventions.common.measure, measure,
        conventions.common.topology, record);
}

double MultiChannelPhaseSpace::DensityIn(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention) const
{
    auto conventions = RequireCompatibleConvention();
    if (conventions.common.topology != convention.topology) {
        std::ostringstream oss;
        oss << "Phase-space topology "
            << PhaseSpaceTopologyName(conventions.common.topology)
            << " does not match the requested density convention "
            << PhaseSpaceTopologyName(convention.topology)
            << " [siren-docs: errors#measure-compat]";
        throw siren::utilities::MeasureCompatibilityError(oss.str());
    }
    double density = ComputeDensity(detector_model, record, conventions);
    return ConvertDensity(
        density, conventions.common.measure, convention.measure,
        convention.topology, record);
}

std::vector<double> MultiChannelPhaseSpace::DensityBreakdown(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    std::vector<double> weighted;
    ComputeDensity(detector_model, record, RequireCompatibleConvention(), &weighted);
    return weighted;
}

PhaseSpaceTopology MultiChannelPhaseSpace::CommonTopology() const {
    return CommonConvention().topology;
}

PhaseSpaceMeasure MultiChannelPhaseSpace::CommonMeasure() const {
    return CommonConvention().measure;
}

PhaseSpaceConvention MultiChannelPhaseSpace::CommonConvention() const {
    return RequireCompatibleConvention().common;
}

std::vector<MultiChannelPhaseSpace::ChannelDiagnostic>
MultiChannelPhaseSpace::ValidateChannelsDetailed() const {
    return EvaluateConventions().diagnostics;
}

void NestedMixtureChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    if (!mixture) {
        throw std::runtime_error("NestedMixtureChannel has no inner mixture");
    }
    mixture->Sample(random, detector_model, record);
}

double NestedMixtureChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!mixture) return 0.0;
    return mixture->Density(detector_model, record);
}

std::string NestedMixtureChannel::Name() const {
    return label;
}

PhaseSpaceTopology NestedMixtureChannel::Topology() const {
    if (!mixture) return PhaseSpaceTopology::Unspecified;
    return mixture->CommonTopology();
}

PhaseSpaceMeasure NestedMixtureChannel::Measure() const {
    if (!mixture) return PhaseSpaceMeasure::Unspecified();
    return mixture->CommonMeasure();
}

} // namespace injection
} // namespace siren
