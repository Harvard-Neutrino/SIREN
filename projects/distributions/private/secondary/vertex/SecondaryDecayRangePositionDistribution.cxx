#include "SIREN/distributions/secondary/vertex/SecondaryDecayRangePositionDistribution.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace distributions {

using detector::DetectorDirection;
using detector::DetectorPosition;

namespace {

struct Interval {
    double lower;
    double upper;
};

struct TransportParameters {
    std::vector<siren::dataclasses::ParticleType> targets;
    std::vector<double> total_cross_sections;
    double total_decay_length = std::numeric_limits<double>::infinity();
};

struct DepthSegment {
    double lower_depth;
    double upper_depth;
    double log_maximum_success;
    double log_envelope_weight;
    double log_integral;
    bool splittable = true;
};

double LogOneMinusExpNegative(double depth) {
    if(!(depth > 0.0)) {
        return -std::numeric_limits<double>::infinity();
    }
    if(depth < std::log(2.0)) {
        return std::log(-std::expm1(-depth));
    }
    return std::log1p(-std::exp(-depth));
}

double LogAddExp(double left, double right) {
    if(left == -std::numeric_limits<double>::infinity()) {
        return right;
    }
    if(right == -std::numeric_limits<double>::infinity()) {
        return left;
    }
    double high = std::max(left, right);
    double low = std::min(left, right);
    return high + std::log1p(std::exp(low - high));
}

template<typename Function>
double AdaptiveSimpsonRecursive(
    Function const & function,
    double lower,
    double upper,
    double lower_value,
    double midpoint_value,
    double upper_value,
    double estimate,
    double tolerance,
    unsigned int depth) {
    double midpoint = 0.5 * (lower + upper);
    double left_midpoint = 0.5 * (lower + midpoint);
    double right_midpoint = 0.5 * (midpoint + upper);
    double left_midpoint_value = function(left_midpoint);
    double right_midpoint_value = function(right_midpoint);
    double left_estimate = (midpoint - lower)
        * (lower_value + 4.0 * left_midpoint_value + midpoint_value) / 6.0;
    double right_estimate = (upper - midpoint)
        * (midpoint_value + 4.0 * right_midpoint_value + upper_value) / 6.0;
    double refined = left_estimate + right_estimate;
    double error = refined - estimate;
    if(depth == 0 || std::abs(error) <= 15.0 * tolerance) {
        return refined + error / 15.0;
    }
    return AdaptiveSimpsonRecursive(
               function, lower, midpoint,
               lower_value, left_midpoint_value, midpoint_value,
               left_estimate, 0.5 * tolerance, depth - 1)
        + AdaptiveSimpsonRecursive(
               function, midpoint, upper,
               midpoint_value, right_midpoint_value, upper_value,
               right_estimate, 0.5 * tolerance, depth - 1);
}

template<typename Function>
double AdaptiveIntegral(Function const & function) {
    constexpr double lower = 0.0;
    constexpr double upper = 1.0;
    double lower_value = function(lower);
    double midpoint_value = function(0.5);
    double upper_value = function(upper);
    double estimate = (lower_value + 4.0 * midpoint_value + upper_value) / 6.0;
    double scale = std::max({
        std::abs(lower_value), std::abs(midpoint_value),
        std::abs(upper_value)});
    if(scale == 0.0) {
        return 0.0;
    }
    return AdaptiveSimpsonRecursive(
        function, lower, upper,
        lower_value, midpoint_value, upper_value,
        estimate, 1e-9 * scale, 20);
}

std::vector<Interval> FiducialIntervals(
    std::shared_ptr<siren::geometry::Geometry> const & fiducial_volume,
    siren::math::Vector3D const & origin,
    siren::math::Vector3D const & direction,
    double max_length) {
    std::vector<Interval> intervals;
    if(!fiducial_volume || !(max_length > 0.0)) {
        return intervals;
    }

    auto intersections = fiducial_volume->Intersections(origin, direction);
    std::sort(intersections.begin(), intersections.end(),
        [](auto const & left, auto const & right) {
            return left.distance < right.distance;
        });

    bool active = fiducial_volume->IsInside(origin, direction);
    double lower = 0.0;
    constexpr double distance_tolerance = 1e-9;

    for(auto const & intersection : intersections) {
        double distance = intersection.distance;
        if(distance <= distance_tolerance) {
            continue;
        }
        if(distance >= max_length) {
            break;
        }

        if(intersection.entering) {
            if(!active) {
                lower = distance;
                active = true;
            }
        } else if(active) {
            if(distance > lower) {
                intervals.push_back({lower, distance});
            }
            active = false;
        }
    }

    if(active && std::isfinite(max_length) && max_length > lower) {
        intervals.push_back({lower, max_length});
    }
    return intervals;
}

TransportParameters MakeTransportParameters(
    std::shared_ptr<siren::detector::DetectorModel const> const & detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> const & interactions,
    siren::dataclasses::InteractionRecord const & record) {
    TransportParameters result;
    if(!interactions) {
        return result;
    }

    result.targets.assign(
        interactions->TargetTypes().begin(), interactions->TargetTypes().end());
    result.total_cross_sections.assign(result.targets.size(), 0.0);
    result.total_decay_length = interactions->TotalDecayLengthAllFinalStates(record);

    siren::dataclasses::InteractionRecord target_record = record;
    for(std::size_t i = 0; i < result.targets.size(); ++i) {
        auto target = result.targets[i];
        target_record.signature.target_type = target;
        target_record.target_mass = detector_model->GetTargetMass(target);
        for(auto const & cross_section : interactions->GetCrossSectionsForTarget(target)) {
            double value = cross_section->TotalCrossSectionAllFinalStates(target_record);
            if(std::isfinite(value) && value > 0.0) {
                result.total_cross_sections[i] += value;
            }
        }
    }
    return result;
}

siren::dataclasses::InteractionRecord MakeDaughterRecord(
    siren::dataclasses::InteractionRecord const & parent_record,
    std::shared_ptr<siren::interactions::InteractionCollection> const & daughter_interactions,
    double daughter_mass,
    double daughter_energy_fraction,
    siren::math::Vector3D const & direction) {
    siren::dataclasses::InteractionRecord daughter_record;
    daughter_record.signature.primary_type = daughter_interactions->GetPrimaryType();
    daughter_record.primary_initial_position = parent_record.primary_initial_position;
    daughter_record.interaction_vertex = parent_record.primary_initial_position;
    daughter_record.primary_mass = daughter_mass;
    daughter_record.primary_helicity = parent_record.primary_helicity;

    double parent_energy = parent_record.primary_momentum[0];
    double daughter_energy = daughter_energy_fraction * parent_energy;
    daughter_energy = std::max(daughter_energy, daughter_mass);
    double momentum_squared = std::max(
        0.0, daughter_energy * daughter_energy - daughter_mass * daughter_mass);
    double momentum = std::sqrt(momentum_squared);
    daughter_record.primary_momentum = {
        daughter_energy,
        momentum * direction.GetX(),
        momentum * direction.GetY(),
        momentum * direction.GetZ()};
    return daughter_record;
}

void AddEdge(std::vector<double> & edges, double value,
             double lower, double upper) {
    if(value > lower && value < upper && std::isfinite(value)) {
        edges.push_back(value);
    }
}

} // namespace

struct SecondaryDecayRangePositionDistribution::DensityContext {
    siren::math::Vector3D origin;
    siren::math::Vector3D direction;
    double support_lower = 0.0;
    double support_upper = 0.0;
    std::shared_ptr<siren::detector::Path> path;
    std::vector<Interval> fiducial_intervals;
    TransportParameters parent_transport;
    TransportParameters daughter_transport;
    std::vector<DepthSegment> segments;
    double log_normalization = -std::numeric_limits<double>::infinity();
    double log_envelope_normalization = -std::numeric_limits<double>::infinity();
};

SecondaryDecayRangePositionDistribution::
SecondaryDecayRangePositionDistribution() = default;

SecondaryDecayRangePositionDistribution::
SecondaryDecayRangePositionDistribution(
    std::shared_ptr<siren::geometry::Geometry> fiducial_volume,
    std::shared_ptr<siren::interactions::InteractionCollection> daughter_interactions,
    double daughter_mass,
    double daughter_energy_fraction,
    double max_length)
    : fiducial_volume_(std::move(fiducial_volume)),
      daughter_interactions_(std::move(daughter_interactions)),
      daughter_mass_(daughter_mass),
      daughter_energy_fraction_(daughter_energy_fraction),
      max_length_(max_length) {
    if(!fiducial_volume_) {
        throw std::invalid_argument(
            "SecondaryDecayRangePositionDistribution requires a fiducial volume");
    }
    if(!daughter_interactions_) {
        throw std::invalid_argument(
            "SecondaryDecayRangePositionDistribution requires daughter interactions");
    }
    if(!(daughter_mass_ > 0.0) || !std::isfinite(daughter_mass_)) {
        throw std::invalid_argument(
            "SecondaryDecayRangePositionDistribution requires a finite positive daughter mass");
    }
    if(!(daughter_energy_fraction_ > 0.0)
       || !(daughter_energy_fraction_ <= 1.0)) {
        throw std::invalid_argument(
            "SecondaryDecayRangePositionDistribution requires a daughter energy fraction in (0, 1]");
    }
    if(!(max_length_ > 0.0)) {
        throw std::invalid_argument(
            "SecondaryDecayRangePositionDistribution requires a positive maximum length");
    }
}

double SecondaryDecayRangePositionDistribution::LogDaughterSuccessProbability(
    DensityContext const & context, double distance) const {
    auto depth_at = [&context](double position) {
        return context.path->GetInteractionDepthFromStartInBounds(
            position - context.support_lower,
            context.daughter_transport.targets,
            context.daughter_transport.total_cross_sections,
            context.daughter_transport.total_decay_length);
    };

    double vertex_depth = depth_at(distance);
    double log_probability = -std::numeric_limits<double>::infinity();
    for(auto const & interval : context.fiducial_intervals) {
        double lower = std::max(distance, interval.lower);
        if(!(interval.upper > lower)) {
            continue;
        }
        double lower_depth = depth_at(lower);
        double upper_depth = depth_at(interval.upper);
        double survival_depth = std::max(0.0, lower_depth - vertex_depth);
        double process_depth = std::max(0.0, upper_depth - lower_depth);
        double log_term = -survival_depth
            + LogOneMinusExpNegative(process_depth);
        log_probability = LogAddExp(log_probability, log_term);
    }
    return std::min(0.0, log_probability);
}

std::shared_ptr<SecondaryDecayRangePositionDistribution::DensityContext>
SecondaryDecayRangePositionDistribution::BuildDensityContext(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    auto context = std::make_shared<DensityContext>();
    context->origin = siren::math::Vector3D(record.primary_initial_position);
    context->direction = siren::math::Vector3D(
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3]);
    if(!detector_model || !interactions || !fiducial_volume_
       || !daughter_interactions_ || context->direction.magnitude() == 0.0) {
        return context;
    }
    context->direction.normalize();

    auto fiducial_intervals = FiducialIntervals(
        fiducial_volume_, context->origin, context->direction, max_length_);
    if(fiducial_intervals.empty()) {
        return context;
    }

    double ray_end = fiducial_intervals.back().upper;
    if(!(ray_end > 0.0) || !std::isfinite(ray_end)) {
        return context;
    }

    siren::detector::Path initial_path(
        detector_model, DetectorPosition(context->origin),
        DetectorDirection(context->direction), ray_end);
    initial_path.ClipToOuterBounds();
    context->support_lower =
        (initial_path.GetFirstPoint().get() - context->origin) * context->direction;
    context->support_upper =
        (initial_path.GetLastPoint().get() - context->origin) * context->direction;
    context->support_lower = std::max(0.0, context->support_lower);
    context->support_upper = std::min(ray_end, context->support_upper);
    if(!(context->support_upper > context->support_lower)) {
        return context;
    }

    for(auto const & interval : fiducial_intervals) {
        double lower = std::max(interval.lower, context->support_lower);
        double upper = std::min(interval.upper, context->support_upper);
        if(upper > lower) {
            context->fiducial_intervals.push_back({lower, upper});
        }
    }
    if(context->fiducial_intervals.empty()) {
        return context;
    }
    context->support_upper = context->fiducial_intervals.back().upper;

    context->path = std::make_shared<siren::detector::Path>(
        detector_model,
        DetectorPosition(context->origin
            + context->support_lower * context->direction),
        DetectorPosition(context->origin
            + context->support_upper * context->direction));
    context->path->EnsureIntersections();

    context->parent_transport = MakeTransportParameters(
        detector_model, interactions, record);
    auto daughter_record = MakeDaughterRecord(
        record, daughter_interactions_, daughter_mass_,
        daughter_energy_fraction_, context->direction);
    context->daughter_transport = MakeTransportParameters(
        detector_model, daughter_interactions_, daughter_record);

    double total_parent_depth = context->path->GetInteractionDepthInBounds(
        context->parent_transport.targets,
        context->parent_transport.total_cross_sections,
        context->parent_transport.total_decay_length);
    if(!(total_parent_depth > 0.0) || !std::isfinite(total_parent_depth)) {
        return context;
    }

    std::vector<double> natural_edges{
        context->support_lower, context->support_upper};
    for(auto const & interval : context->fiducial_intervals) {
        AddEdge(natural_edges, interval.lower,
                context->support_lower, context->support_upper);
        AddEdge(natural_edges, interval.upper,
                context->support_lower, context->support_upper);
    }
    for(auto const & intersection : context->path->GetIntersections().intersections) {
        AddEdge(natural_edges,
                context->support_lower + intersection.distance,
                context->support_lower, context->support_upper);
    }
    std::sort(natural_edges.begin(), natural_edges.end());
    double edge_tolerance = 1e-12 * std::max(
        1.0, context->support_upper - context->support_lower);
    natural_edges.erase(
        std::unique(natural_edges.begin(), natural_edges.end(),
            [edge_tolerance](double left, double right) {
                return std::abs(left - right) <= edge_tolerance;
            }),
        natural_edges.end());

    auto parent_depth_at = [&context](double distance) {
        return context->path->GetInteractionDepthFromStartInBounds(
            distance - context->support_lower,
            context->parent_transport.targets,
            context->parent_transport.total_cross_sections,
            context->parent_transport.total_decay_length);
    };
    auto distance_at_parent_depth = [&context](double depth) {
        return context->support_lower
            + context->path->GetDistanceFromStartInBounds(
                depth,
                context->parent_transport.targets,
                context->parent_transport.total_cross_sections,
                context->parent_transport.total_decay_length);
    };
    auto daughter_depth_at = [&context](double distance) {
        return context->path->GetInteractionDepthFromStartInBounds(
            distance - context->support_lower,
            context->daughter_transport.targets,
            context->daughter_transport.total_cross_sections,
            context->daughter_transport.total_decay_length);
    };

    auto evaluate_segment = [&](double lower_depth, double upper_depth) {
        double depth_width = upper_depth - lower_depth;
        double lower_distance = distance_at_parent_depth(lower_depth);
        double upper_distance = distance_at_parent_depth(upper_depth);
        double middle_distance = distance_at_parent_depth(
            0.5 * (lower_depth + upper_depth));
        double log_maximum_success = std::max({
            LogDaughterSuccessProbability(*context, lower_distance),
            LogDaughterSuccessProbability(*context, middle_distance),
            LogDaughterSuccessProbability(*context, upper_distance)});
        if(!std::isfinite(log_maximum_success)) {
            return DepthSegment{
                lower_depth, upper_depth,
                -std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity(),
                -std::numeric_limits<double>::infinity(), false};
        }

        // Between fiducial and material boundaries the daughter's success
        // probability is monotone.  Endpoint values therefore bound the
        // segment; the midpoint protects against roundoff at a boundary.
        log_maximum_success = std::min(
            0.0, log_maximum_success
                + 4.0 * std::numeric_limits<double>::epsilon());
        double log_envelope_weight = -lower_depth
            + LogOneMinusExpNegative(depth_width)
            + log_maximum_success;

        auto scaled_integrand = [&, lower_depth, depth_width,
                                  log_maximum_success](
            double unit_coordinate) {
            double depth = lower_depth + depth_width * unit_coordinate;
            double distance = distance_at_parent_depth(depth);
            double log_success = LogDaughterSuccessProbability(
                *context, distance);
            return std::exp(
                -(depth - lower_depth)
                + log_success - log_maximum_success);
        };
        double scaled_integral = depth_width
            * AdaptiveIntegral(scaled_integrand);
        double envelope_integral = -std::expm1(-depth_width);
        scaled_integral = std::min(scaled_integral, envelope_integral);
        double log_integral = -std::numeric_limits<double>::infinity();
        if(scaled_integral > 0.0 && std::isfinite(scaled_integral)) {
            log_integral = -lower_depth
                + log_maximum_success + std::log(scaled_integral);
        }
        return DepthSegment{
            lower_depth, upper_depth, log_maximum_success,
            log_envelope_weight, log_integral, true};
    };

    std::vector<DepthSegment> working_segments;
    for(std::size_t i = 0; i + 1 < natural_edges.size(); ++i) {
        double lower_depth = parent_depth_at(natural_edges[i]);
        double upper_depth = parent_depth_at(natural_edges[i + 1]);
        if(upper_depth > lower_depth) {
            auto segment = evaluate_segment(lower_depth, upper_depth);
            if(std::isfinite(segment.log_envelope_weight)) {
                working_segments.push_back(segment);
            }
        }
    }

    // Refine only segments whose loose envelope materially hurts the global
    // rejection efficiency.  Splits are chosen in daughter interaction depth,
    // then mapped back through Path; no spatial integration grid is needed.
    constexpr std::size_t maximum_refinements = 4096;
    constexpr double target_log_acceptance = -1.3862943611198906; // log(1/4)
    double minimum_depth_width = 1e-12 * std::max(1.0, total_parent_depth);
    for(std::size_t refinement = 0;
        refinement < maximum_refinements && !working_segments.empty();
        ++refinement) {
        double log_integral = -std::numeric_limits<double>::infinity();
        double log_envelope = -std::numeric_limits<double>::infinity();
        for(auto const & segment : working_segments) {
            log_integral = LogAddExp(log_integral, segment.log_integral);
            log_envelope = LogAddExp(
                log_envelope, segment.log_envelope_weight);
        }
        if(std::isfinite(log_integral)
           && log_integral - log_envelope >= target_log_acceptance) {
            break;
        }

        std::size_t selected_index = working_segments.size();
        double selected_log_waste = -std::numeric_limits<double>::infinity();
        for(std::size_t i = 0; i < working_segments.size(); ++i) {
            auto const & segment = working_segments[i];
            if(!segment.splittable
               || segment.upper_depth - segment.lower_depth
                    <= minimum_depth_width) {
                continue;
            }
            double log_waste = segment.log_envelope_weight;
            if(std::isfinite(segment.log_integral)) {
                double log_ratio = std::max(
                    0.0, segment.log_envelope_weight - segment.log_integral);
                log_waste = segment.log_envelope_weight
                    + LogOneMinusExpNegative(log_ratio);
            }
            if(log_waste > selected_log_waste) {
                selected_log_waste = log_waste;
                selected_index = i;
            }
        }
        if(selected_index == working_segments.size()) {
            break;
        }

        auto & selected = working_segments[selected_index];
        double lower_distance = distance_at_parent_depth(selected.lower_depth);
        double upper_distance = distance_at_parent_depth(selected.upper_depth);
        double lower_daughter_depth = daughter_depth_at(lower_distance);
        double upper_daughter_depth = daughter_depth_at(upper_distance);
        double split_depth = 0.5 * (
            selected.lower_depth + selected.upper_depth);
        if(upper_daughter_depth > lower_daughter_depth) {
            double daughter_midpoint = 0.5 * (
                lower_daughter_depth + upper_daughter_depth);
            double daughter_distance = context->support_lower
                + context->path->GetDistanceFromStartInBounds(
                    daughter_midpoint,
                    context->daughter_transport.targets,
                    context->daughter_transport.total_cross_sections,
                    context->daughter_transport.total_decay_length);
            double candidate = parent_depth_at(daughter_distance);
            if(candidate > selected.lower_depth + minimum_depth_width
               && candidate < selected.upper_depth - minimum_depth_width) {
                split_depth = candidate;
            }
        }
        if(!(split_depth > selected.lower_depth + minimum_depth_width
             && split_depth < selected.upper_depth - minimum_depth_width)) {
            selected.splittable = false;
            continue;
        }

        auto left = evaluate_segment(selected.lower_depth, split_depth);
        auto right = evaluate_segment(split_depth, selected.upper_depth);
        working_segments[selected_index] = left;
        working_segments.insert(
            working_segments.begin() + selected_index + 1, right);
    }

    context->segments = std::move(working_segments);
    for(auto const & segment : context->segments) {
        context->log_normalization = LogAddExp(
            context->log_normalization, segment.log_integral);
        context->log_envelope_normalization = LogAddExp(
            context->log_envelope_normalization,
            segment.log_envelope_weight);
    }

    if(context->segments.empty()
       || !std::isfinite(context->log_normalization)
       || !std::isfinite(context->log_envelope_normalization)) {
        context->segments.clear();
    }
    return context;
}

void SecondaryDecayRangePositionDistribution::SampleVertex(
    std::shared_ptr<siren::utilities::SIREN_random> rand,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::SecondaryDistributionRecord & record) const {
    auto context = BuildDensityContext(detector_model, interactions, record.record);
    if(context->segments.empty()) {
        throw siren::utilities::InjectionFailure(
            siren::utilities::FailureReason::NoPathThroughVolume,
            "No parent vertex can produce a daughter process in the fiducial volume");
    }

    constexpr std::size_t maximum_attempts = 100000;
    for(std::size_t attempt = 0; attempt < maximum_attempts; ++attempt) {
        double draw = rand->Uniform();
        DepthSegment const * selected = &context->segments.back();
        double cumulative = 0.0;
        for(auto const & segment : context->segments) {
            cumulative += std::exp(
                segment.log_envelope_weight
                - context->log_envelope_normalization);
            if(draw <= cumulative) {
                selected = &segment;
                break;
            }
        }

        double depth_width = selected->upper_depth - selected->lower_depth;
        double mass = -std::expm1(-depth_width);
        double depth = selected->lower_depth
            - std::log1p(-rand->Uniform() * mass);
        double distance = context->support_lower
            + context->path->GetDistanceFromStartInBounds(
                depth,
                context->parent_transport.targets,
                context->parent_transport.total_cross_sections,
                context->parent_transport.total_decay_length);
        double log_success = LogDaughterSuccessProbability(*context, distance);
        double log_acceptance_draw = std::log(rand->Uniform())
            + selected->log_maximum_success;
        if(log_acceptance_draw <= log_success) {
            record.SetLength(distance);
            return;
        }
    }

    // Not a geometry failure: a path exists, but the refined proposal
    // envelope never reached a workable acceptance rate.
    throw siren::utilities::InjectionFailure(
        siren::utilities::FailureReason::Unspecified,
        "Conditional parent-vertex rejection sampling exhausted its attempt budget");
}

double SecondaryDecayRangePositionDistribution::GenerationProbability(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    auto context = BuildDensityContext(detector_model, interactions, record);
    if(context->segments.empty()) {
        return 0.0;
    }

    siren::math::Vector3D vertex(record.interaction_vertex);
    siren::math::Vector3D displacement = vertex - context->origin;
    double distance = displacement * context->direction;
    siren::math::Vector3D transverse =
        displacement - distance * context->direction;
    if(transverse.magnitude() > 1e-6
       || distance < context->support_lower
       || distance > context->support_upper) {
        return 0.0;
    }

    double path_distance = distance - context->support_lower;
    double parent_depth = context->path->GetInteractionDepthFromStartInBounds(
        path_distance,
        context->parent_transport.targets,
        context->parent_transport.total_cross_sections,
        context->parent_transport.total_decay_length);
    double log_success = LogDaughterSuccessProbability(*context, distance);
    double parent_density = detector_model->GetInteractionDensity(
        context->path->GetIntersections(), DetectorPosition(vertex),
        context->parent_transport.targets,
        context->parent_transport.total_cross_sections,
        context->parent_transport.total_decay_length);
    if(!(parent_density > 0.0) || !std::isfinite(log_success)
       || !std::isfinite(parent_density) || !std::isfinite(parent_depth)) {
        return 0.0;
    }

    double log_density = std::log(parent_density)
        - parent_depth + log_success - context->log_normalization;
    return std::exp(log_density);
}

std::string SecondaryDecayRangePositionDistribution::Name() const {
    return "SecondaryDecayRangePositionDistribution";
}

std::shared_ptr<SecondaryInjectionDistribution>
SecondaryDecayRangePositionDistribution::clone() const {
    return std::make_shared<SecondaryDecayRangePositionDistribution>(*this);
}

std::tuple<siren::math::Vector3D, siren::math::Vector3D>
SecondaryDecayRangePositionDistribution::InjectionBounds(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record) const {
    auto context = BuildDensityContext(detector_model, interactions, record);
    if(context->segments.empty()) {
        return {siren::math::Vector3D(), siren::math::Vector3D()};
    }
    return {
        context->origin + context->support_lower * context->direction,
        context->origin + context->support_upper * context->direction};
}

bool SecondaryDecayRangePositionDistribution::equal(
    WeightableDistribution const & other) const {
    auto const * distribution =
        dynamic_cast<SecondaryDecayRangePositionDistribution const *>(&other);
    if(!distribution) {
        return false;
    }
    bool same_fiducial = (!fiducial_volume_ && !distribution->fiducial_volume_)
        || (fiducial_volume_ && distribution->fiducial_volume_
            && *fiducial_volume_ == *distribution->fiducial_volume_);
    // InteractionCollection defines no ordering, so equal() and less() must
    // both compare the daughter collection by identity to stay mutually
    // consistent; weighting-side distribution cancellation is identity-based
    // as well.
    bool same_daughter =
        daughter_interactions_ == distribution->daughter_interactions_;
    return same_fiducial && same_daughter
        && std::tie(daughter_mass_, daughter_energy_fraction_, max_length_)
           == std::tie(distribution->daughter_mass_,
                       distribution->daughter_energy_fraction_,
                       distribution->max_length_);
}

bool SecondaryDecayRangePositionDistribution::less(
    WeightableDistribution const & other) const {
    auto const * distribution =
        dynamic_cast<SecondaryDecayRangePositionDistribution const *>(&other);
    if(!distribution) {
        return false;
    }
    auto values = std::tie(
        daughter_mass_, daughter_energy_fraction_, max_length_);
    auto other_values = std::tie(
        distribution->daughter_mass_,
        distribution->daughter_energy_fraction_,
        distribution->max_length_);
    if(values != other_values) {
        return values < other_values;
    }
    if(static_cast<bool>(fiducial_volume_)
       != static_cast<bool>(distribution->fiducial_volume_)) {
        return static_cast<bool>(fiducial_volume_)
            < static_cast<bool>(distribution->fiducial_volume_);
    }
    if(fiducial_volume_ && *fiducial_volume_ != *distribution->fiducial_volume_) {
        return *fiducial_volume_ < *distribution->fiducial_volume_;
    }
    return std::owner_less<
        std::shared_ptr<siren::interactions::InteractionCollection>>()(
            daughter_interactions_, distribution->daughter_interactions_);
}

} // namespace distributions
} // namespace siren
