#include <cmath>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/json.hpp>
#include <gtest/gtest.h>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/CartesianAxis1D.h"
#include "SIREN/detector/CartesianAxisPolynomialDensityDistribution.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/detector/MaterialModel.h"
#include "SIREN/detector/Path.h"
#include "SIREN/detector/PolynomialDistribution1D.h"
#include "SIREN/distributions/secondary/vertex/SecondaryDecayRangePositionDistribution.h"
#include "SIREN/distributions/serializable.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/HNLDipoleDecay.h"
#include "SIREN/interactions/DummyCrossSection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

namespace {

using siren::dataclasses::InteractionRecord;
using siren::dataclasses::InteractionSignature;
using siren::dataclasses::ParticleType;

class FixedLengthDecay : public siren::interactions::Decay {
public:
    explicit FixedLengthDecay(double length) : length_(length) {}

    bool equal(siren::interactions::Decay const & other) const override {
        auto const * decay = dynamic_cast<FixedLengthDecay const *>(&other);
        return decay && decay->length_ == length_;
    }
    double TotalDecayWidthAllFinalStates(InteractionRecord const &) const override {
        return 1.0;
    }
    double TotalDecayWidth(ParticleType) const override { return 1.0; }
    double TotalDecayWidth(InteractionRecord const &) const override { return 1.0; }
    double TotalDecayLengthAllFinalStates(InteractionRecord const &) const override {
        return length_;
    }
    double TotalDecayLength(InteractionRecord const &) const override {
        return length_;
    }
    double DifferentialDecayWidth(InteractionRecord const &) const override {
        return 1.0;
    }
    void SampleFinalState(
        siren::dataclasses::CrossSectionDistributionRecord &,
        std::shared_ptr<siren::utilities::SIREN_random>) const override {}
    std::vector<InteractionSignature> GetPossibleSignatures() const override {
        return {};
    }
    std::vector<InteractionSignature> GetPossibleSignaturesFromParent(
        ParticleType) const override {
        return {};
    }
    double FinalStateProbability(InteractionRecord const &) const override {
        return 1.0;
    }
    std::vector<std::string> DensityVariables() const override { return {}; }

private:
    double length_;
};

std::shared_ptr<siren::detector::DetectorModel> UniformDetector() {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = 0;
    world.level = 0;
    world.geo = siren::geometry::Sphere(100.0, 0.0).create();
    world.density = siren::detector::ConstantDensityDistribution(1.0).create();
    detector->AddSector(world);
    return detector;
}

std::shared_ptr<siren::detector::DetectorModel> LayeredDetector() {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    siren::detector::MaterialModel materials;
    materials.AddMaterial("hydrogen", {{1000010010, 1.0}});
    detector->SetMaterials(materials);

    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = materials.GetMaterialId("hydrogen");
    world.level = 0;
    world.geo = siren::geometry::Sphere(100.0, 0.0).create();
    world.density = siren::detector::ConstantDensityDistribution(1.0).create();
    detector->AddSector(world);

    siren::detector::DetectorSector dense_layer;
    dense_layer.name = "dense_layer";
    dense_layer.material_id = materials.GetMaterialId("hydrogen");
    dense_layer.level = 1;
    dense_layer.geo = siren::geometry::Box(
        siren::geometry::Placement(siren::math::Vector3D(3.0, 0.0, 0.0)),
        2.0, 10.0, 10.0).create();
    dense_layer.density =
        siren::detector::ConstantDensityDistribution(10.0).create();
    detector->AddSector(dense_layer);
    return detector;
}

std::shared_ptr<siren::detector::DetectorModel> PolynomialDensityDetector() {
    auto detector = std::make_shared<siren::detector::DetectorModel>();
    detector->ClearSectors();
    siren::detector::MaterialModel materials;
    materials.AddMaterial("hydrogen", {{1000010010, 1.0}});
    detector->SetMaterials(materials);

    siren::detector::DetectorSector world;
    world.name = "world";
    world.material_id = materials.GetMaterialId("hydrogen");
    world.level = 0;
    world.geo = siren::geometry::Sphere(100.0, 0.0).create();
    siren::detector::CartesianAxis1D axis(
        siren::math::Vector3D(1.0, 0.0, 0.0),
        siren::math::Vector3D(0.0, 0.0, 0.0));
    // rho(x) = 1 + 0.04 x^2 is positive throughout the sampled ray and makes
    // distance and interaction depth distinctly nonlinear.
    siren::detector::PolynomialDistribution1D density({1.0, 0.0, 0.04});
    world.density =
        siren::detector::CartesianAxisPolynomialDensityDistribution(
            axis, density).create();
    detector->AddSector(world);
    return detector;
}

std::shared_ptr<siren::geometry::Geometry> FiducialBox() {
    return siren::geometry::Box(
        siren::geometry::Placement(siren::math::Vector3D(6.0, 0.0, 0.0)),
        4.0, 10.0, 10.0).create();
}

// Union of two disjoint boxes spanning x in [4, 8] and [12, 16], so the ray
// from the origin along +x crosses two separated fiducial intervals.
std::shared_ptr<siren::geometry::Geometry> TwoBoxFiducial() {
    auto near_box = siren::geometry::Box(
        siren::geometry::Placement(siren::math::Vector3D(6.0, 0.0, 0.0)),
        4.0, 10.0, 10.0).create();
    auto far_box = siren::geometry::Box(
        siren::geometry::Placement(siren::math::Vector3D(14.0, 0.0, 0.0)),
        4.0, 10.0, 10.0).create();
    return siren::geometry::BooleanGeometry(
        siren::geometry::BooleanOperation::UNION, near_box, far_box).create();
}

// Probability that the daughter's first decay falls in either interval,
// given a start at `position` in a uniform unit-density medium.  Survival to
// the second interval correctly includes traversal of the first.
double TwoIntervalDaughterSuccess(double position, double daughter_length) {
    constexpr double first_lower = 4.0;
    constexpr double first_upper = 8.0;
    constexpr double second_lower = 12.0;
    constexpr double second_upper = 16.0;
    auto survive = [daughter_length](double from, double to) {
        return std::exp(-std::max(0.0, to - from) / daughter_length);
    };
    double result = 0.0;
    if(position < first_upper) {
        double start = std::max(position, first_lower);
        result += survive(position, start)
            * (1.0 - survive(start, first_upper));
    }
    if(position < second_upper) {
        double start = std::max(position, second_lower);
        result += survive(position, start)
            * (1.0 - survive(start, second_upper));
    }
    return result;
}

double TwoIntervalTarget(double position, double parent_length,
                         double daughter_length) {
    return std::exp(-position / parent_length) / parent_length
        * TwoIntervalDaughterSuccess(position, daughter_length);
}

double SimpsonIntegral(std::function<double(double)> const & function,
                       double lower, double upper, std::size_t intervals);

double TwoIntervalTargetIntegral(double lower, double upper,
                                 double parent_length,
                                 double daughter_length) {
    // Piecewise Simpson between the analytic kinks at the interval edges.
    constexpr double edges[] = {0.0, 4.0, 8.0, 12.0, 16.0};
    auto target = [&](double position) {
        return TwoIntervalTarget(position, parent_length, daughter_length);
    };
    double result = 0.0;
    for(std::size_t i = 0; i + 1 < 5; ++i) {
        double piece_lower = std::max(lower, edges[i]);
        double piece_upper = std::min(upper, edges[i + 1]);
        if(piece_upper > piece_lower) {
            result += SimpsonIntegral(target, piece_lower, piece_upper, 128);
        }
    }
    return result;
}

std::shared_ptr<siren::interactions::InteractionCollection> Decays(
    ParticleType type, double length) {
    std::vector<std::shared_ptr<siren::interactions::Decay>> decays{
        std::make_shared<FixedLengthDecay>(length)};
    return std::make_shared<siren::interactions::InteractionCollection>(
        type, decays);
}

InteractionRecord ParentRecord(double energy = 2.0) {
    InteractionRecord record;
    record.signature.primary_type = ParticleType::NuMu;
    record.primary_initial_position = {0.0, 0.0, 0.0};
    record.primary_mass = 1.0;
    record.primary_momentum = {
        energy, std::sqrt(energy * energy - 1.0), 0.0, 0.0};
    return record;
}

double ConditionalTarget(double position, double parent_length,
                         double daughter_length) {
    constexpr double lower = 4.0;
    constexpr double upper = 8.0;
    double daughter_probability;
    if(position < lower) {
        daughter_probability = std::exp(-(lower - position) / daughter_length)
            - std::exp(-(upper - position) / daughter_length);
    } else {
        daughter_probability =
            1.0 - std::exp(-(upper - position) / daughter_length);
    }
    return std::exp(-position / parent_length)
        * daughter_probability / parent_length;
}

double ExpIntegral(double rate, double lower, double upper) {
    if(!(upper > lower)) {
        return 0.0;
    }
    if(std::abs(rate) < 1e-12) {
        return upper - lower;
    }
    return std::exp(rate * lower)
        * std::expm1(rate * (upper - lower)) / rate;
}

double ConditionalTargetIntegral(double lower, double upper,
                                 double parent_length,
                                 double daughter_length) {
    constexpr double fiducial_lower = 4.0;
    constexpr double fiducial_upper = 8.0;
    lower = std::max(0.0, lower);
    upper = std::min(fiducial_upper, upper);
    if(!(upper > lower)) {
        return 0.0;
    }

    double result = 0.0;
    double mixed_rate = -1.0 / parent_length + 1.0 / daughter_length;
    double outside_upper = std::min(upper, fiducial_lower);
    if(outside_upper > lower) {
        double coefficient = (
            std::exp(-fiducial_lower / daughter_length)
            - std::exp(-fiducial_upper / daughter_length)) / parent_length;
        result += coefficient
            * ExpIntegral(mixed_rate, lower, outside_upper);
    }

    double inside_lower = std::max(lower, fiducial_lower);
    if(upper > inside_lower) {
        result += ExpIntegral(
            -1.0 / parent_length, inside_lower, upper) / parent_length;
        result -= std::exp(-fiducial_upper / daughter_length)
            * ExpIntegral(mixed_rate, inside_lower, upper) / parent_length;
    }
    return result;
}

double SimpsonIntegral(std::function<double(double)> const & function,
                       double lower, double upper, std::size_t intervals) {
    if(intervals % 2 != 0) {
        ++intervals;
    }
    double step = (upper - lower) / static_cast<double>(intervals);
    double sum = function(lower) + function(upper);
    for(std::size_t i = 1; i < intervals; ++i) {
        sum += (i % 2 == 0 ? 2.0 : 4.0)
            * function(lower + static_cast<double>(i) * step);
    }
    return sum * step / 3.0;
}

double IndependentPhysicalTarget(
    std::shared_ptr<siren::detector::DetectorModel> const & detector,
    std::shared_ptr<siren::interactions::InteractionCollection> const & parent,
    InteractionRecord const & record,
    double position,
    double daughter_length) {
    siren::math::Vector3D origin(record.primary_initial_position);
    siren::math::Vector3D direction(
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3]);
    direction.normalize();
    siren::detector::Path path(
        detector,
        siren::detector::DetectorPosition(origin),
        siren::detector::DetectorDirection(direction),
        8.0);
    path.EnsureIntersections();

    std::vector<ParticleType> targets(
        parent->TargetTypes().begin(), parent->TargetTypes().end());
    std::vector<double> cross_sections(targets.size(), 0.0);
    InteractionRecord target_record = record;
    for(std::size_t i = 0; i < targets.size(); ++i) {
        target_record.signature.target_type = targets[i];
        target_record.target_mass = detector->GetTargetMass(targets[i]);
        for(auto const & cross_section : parent->GetCrossSectionsForTarget(targets[i])) {
            cross_sections[i] +=
                cross_section->TotalCrossSectionAllFinalStates(target_record);
        }
    }
    double decay_length = parent->TotalDecayLengthAllFinalStates(record);
    double depth = path.GetInteractionDepthFromStartInBounds(
        position, targets, cross_sections, decay_length);
    auto point = siren::detector::DetectorPosition(
        origin + position * direction);
    double hazard = detector->GetInteractionDensity(
        path.GetIntersections(), point, targets, cross_sections, decay_length);

    double daughter_success;
    if(position < 4.0) {
        daughter_success = std::exp(-(4.0 - position) / daughter_length)
            - std::exp(-(8.0 - position) / daughter_length);
    } else {
        daughter_success = 1.0
            - std::exp(-(8.0 - position) / daughter_length);
    }
    return hazard * std::exp(-depth) * daughter_success;
}

} // namespace

TEST(SecondaryDecayRangePositionDistribution,
     UniformLimitMatchesContinuousConditionalDensity) {
    constexpr double parent_length = 12.0;
    constexpr double daughter_length = 2.5;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    InteractionRecord record = ParentRecord();
    double normalization = ConditionalTargetIntegral(
        0.0, 8.0, parent_length, daughter_length);
    // Irrational-looking positions intentionally avoid any regular spatial
    // table and require GenerationProbability to evaluate the continuous PDF.
    for(double position : {0.371, 2.137, 3.913, 4.071, 5.319, 7.777}) {
        record.interaction_vertex = {position, 0.0, 0.0};
        double expected = ConditionalTarget(
            position, parent_length, daughter_length) / normalization;
        EXPECT_NEAR(
            distribution.GenerationProbability(detector, parent, record),
            expected,
            2e-6 * std::max(1.0, expected))
            << "position = " << position;
    }
}

TEST(SecondaryDecayRangePositionDistribution, ContinuousDensityIsNormalized) {
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, 12.0);
    auto daughter = Decays(ParticleType::N4, 2.5);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    InteractionRecord record = ParentRecord();
    auto density = [&](double position) {
        record.interaction_vertex = {position, 0.0, 0.0};
        return distribution.GenerationProbability(detector, parent, record);
    };
    // Split at the fiducial entrance, where the analytic branch changes.
    double integral = SimpsonIntegral(density, 0.0, 4.0, 24)
        + SimpsonIntegral(density, 4.0, 8.0, 24);
    EXPECT_NEAR(integral, 1.0, 2e-5);

    record.interaction_vertex = {-1e-6, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        distribution.GenerationProbability(detector, parent, record), 0.0);
    record.interaction_vertex = {8.0 + 1e-6, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        distribution.GenerationProbability(detector, parent, record), 0.0);
}

TEST(SecondaryDecayRangePositionDistribution,
     NegligibleParentInteractionDepthMatchesUniformLimit) {
    constexpr double parent_length = 1e12;
    constexpr double daughter_length = 2.0;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    double normalization = ConditionalTargetIntegral(
        0.0, 8.0, parent_length, daughter_length);
    InteractionRecord record = ParentRecord();
    for(double position : {0.417, 3.271, 4.613, 7.319}) {
        record.interaction_vertex = {position, 0.0, 0.0};
        double expected = ConditionalTarget(
            position, parent_length, daughter_length) / normalization;
        EXPECT_NEAR(
            distribution.GenerationProbability(detector, parent, record),
            expected,
            3e-6 * std::max(1.0, expected));
    }
}

TEST(SecondaryDecayRangePositionDistribution, FollowsLayeredMaterialDensity) {
    auto detector = LayeredDetector();
    std::vector<std::shared_ptr<siren::interactions::CrossSection>> cross_sections{
        std::make_shared<siren::interactions::DummyCrossSection>()};
    auto parent = std::make_shared<siren::interactions::InteractionCollection>(
        ParticleType::NuMu, cross_sections);
    auto daughter = Decays(ParticleType::N4, 2.0);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    InteractionRecord record = ParentRecord();
    record.interaction_vertex = {1.05, 0.0, 0.0};
    double sparse_density = distribution.GenerationProbability(
        detector, parent, record);
    record.interaction_vertex = {3.05, 0.0, 0.0};
    double dense_density = distribution.GenerationProbability(
        detector, parent, record);

    EXPECT_GT(sparse_density, 0.0);
    EXPECT_GT(dense_density, 5.0 * sparse_density);
}

TEST(SecondaryDecayRangePositionDistribution,
     NonConstantDensityMatchesIndependentInteractionDepthCalculation) {
    constexpr double daughter_length = 2.0;
    auto detector = PolynomialDensityDetector();
    std::vector<std::shared_ptr<siren::interactions::CrossSection>> cross_sections{
        std::make_shared<siren::interactions::DummyCrossSection>()};
    auto parent = std::make_shared<siren::interactions::InteractionCollection>(
        ParticleType::NuMu, cross_sections);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    // At this energy DummyCrossSection gives an O(1 m) interaction length at
    // unit density, so both the nonlinear depth and its attenuation matter.
    InteractionRecord record = ParentRecord(1e11);
    constexpr double first_position = 0.731;
    constexpr double second_position = 2.617;
    double first_target = IndependentPhysicalTarget(
        detector, parent, record, first_position, daughter_length);
    double second_target = IndependentPhysicalTarget(
        detector, parent, record, second_position, daughter_length);

    record.interaction_vertex = {first_position, 0.0, 0.0};
    double first_density = distribution.GenerationProbability(
        detector, parent, record);
    record.interaction_vertex = {second_position, 0.0, 0.0};
    double second_density = distribution.GenerationProbability(
        detector, parent, record);

    ASSERT_GT(first_target, 0.0);
    ASSERT_GT(second_target, 0.0);
    ASSERT_GT(first_density, 0.0);
    ASSERT_GT(second_density, 0.0);
    EXPECT_NEAR(
        first_density / second_density,
        first_target / second_target,
        2e-5 * std::max(1.0, first_target / second_target));
}

TEST(SecondaryDecayRangePositionDistribution, SamplingMatchesAnalyticCdf) {
    constexpr double parent_length = 12.0;
    constexpr double daughter_length = 2.5;
    constexpr double threshold = 4.0;
    constexpr std::size_t samples = 4000;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    double total = ConditionalTargetIntegral(
        0.0, 8.0, parent_length, daughter_length);
    double expected_below = ConditionalTargetIntegral(
        0.0, threshold, parent_length, daughter_length) / total;

    auto random = std::make_shared<siren::utilities::SIREN_random>(12345);
    std::size_t observed_below = 0;
    for(std::size_t i = 0; i < samples; ++i) {
        InteractionRecord sample_record = ParentRecord();
        siren::dataclasses::SecondaryDistributionRecord secondary(sample_record);
        distribution.SampleVertex(random, detector, parent, secondary);
        if(secondary.GetLength() < threshold) {
            ++observed_below;
        }
    }
    double observed_fraction = static_cast<double>(observed_below)
        / static_cast<double>(samples);
    double standard_error = std::sqrt(
        expected_below * (1.0 - expected_below) / static_cast<double>(samples));
    EXPECT_NEAR(observed_fraction, expected_below, 6.0 * standard_error);
}

TEST(SecondaryDecayRangePositionDistribution,
     MultiIntervalFiducialMatchesAnalyticDensity) {
    constexpr double parent_length = 12.0;
    constexpr double daughter_length = 2.5;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        TwoBoxFiducial(), daughter, 1.0, 1.0, 20.0);

    double normalization = TwoIntervalTargetIntegral(
        0.0, 16.0, parent_length, daughter_length);

    InteractionRecord record = ParentRecord();
    // One probe in each analytic regime: before the first box, inside it, in
    // the gap where only the second box remains reachable, and inside the
    // second box.
    for(double position : {2.317, 5.129, 9.531, 13.717}) {
        record.interaction_vertex = {position, 0.0, 0.0};
        double expected = TwoIntervalTarget(
            position, parent_length, daughter_length) / normalization;
        ASSERT_GT(expected, 0.0);
        EXPECT_NEAR(
            distribution.GenerationProbability(detector, parent, record),
            expected,
            2e-6 * std::max(1.0, expected))
            << "position = " << position;
    }

    auto density = [&](double position) {
        record.interaction_vertex = {position, 0.0, 0.0};
        return distribution.GenerationProbability(detector, parent, record);
    };
    double integral = SimpsonIntegral(density, 0.0, 4.0, 24)
        + SimpsonIntegral(density, 4.0, 8.0, 24)
        + SimpsonIntegral(density, 8.0, 12.0, 24)
        + SimpsonIntegral(density, 12.0, 16.0, 24);
    EXPECT_NEAR(integral, 1.0, 2e-5);

    record.interaction_vertex = {16.0 + 1e-6, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        distribution.GenerationProbability(detector, parent, record), 0.0);
}

TEST(SecondaryDecayRangePositionDistribution,
     MultiIntervalSamplingMatchesAnalyticCdf) {
    constexpr double parent_length = 12.0;
    constexpr double daughter_length = 2.5;
    constexpr double threshold = 10.0;
    constexpr std::size_t samples = 4000;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        TwoBoxFiducial(), daughter, 1.0, 1.0, 20.0);

    double total = TwoIntervalTargetIntegral(
        0.0, 16.0, parent_length, daughter_length);
    double expected_below = TwoIntervalTargetIntegral(
        0.0, threshold, parent_length, daughter_length) / total;

    auto random = std::make_shared<siren::utilities::SIREN_random>(24680);
    std::size_t observed_below = 0;
    for(std::size_t i = 0; i < samples; ++i) {
        InteractionRecord sample_record = ParentRecord();
        siren::dataclasses::SecondaryDistributionRecord secondary(sample_record);
        distribution.SampleVertex(random, detector, parent, secondary);
        if(secondary.GetLength() < threshold) {
            ++observed_below;
        }
    }
    double observed_fraction = static_cast<double>(observed_below)
        / static_cast<double>(samples);
    double standard_error = std::sqrt(
        expected_below * (1.0 - expected_below)
        / static_cast<double>(samples));
    EXPECT_NEAR(observed_fraction, expected_below, 6.0 * standard_error);
}

TEST(SecondaryDecayRangePositionDistribution,
     OffAxisVertexHasZeroDensityAndBoundsSpanFiducialIntervals) {
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, 12.0);
    auto daughter = Decays(ParticleType::N4, 2.5);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        TwoBoxFiducial(), daughter, 1.0, 1.0, 20.0);

    InteractionRecord record = ParentRecord();
    record.interaction_vertex = {5.0, 0.5, 0.0};
    EXPECT_DOUBLE_EQ(
        distribution.GenerationProbability(detector, parent, record), 0.0);

    auto bounds = distribution.InjectionBounds(detector, parent, record);
    siren::math::Vector3D lower = std::get<0>(bounds);
    siren::math::Vector3D upper = std::get<1>(bounds);
    EXPECT_NEAR(lower.GetX(), 0.0, 1e-9);
    EXPECT_NEAR(lower.GetY(), 0.0, 1e-9);
    EXPECT_NEAR(lower.GetZ(), 0.0, 1e-9);
    EXPECT_NEAR(upper.GetX(), 16.0, 1e-9);
    EXPECT_NEAR(upper.GetY(), 0.0, 1e-9);
    EXPECT_NEAR(upper.GetZ(), 0.0, 1e-9);
}

TEST(SecondaryDecayRangePositionDistribution,
     ComparisonOperatorsAreMutuallyConsistent) {
    auto shared_daughter = Decays(ParticleType::N4, 2.5);
    siren::distributions::SecondaryDecayRangePositionDistribution first(
        FiducialBox(), shared_daughter, 1.0, 1.0, 20.0);
    siren::distributions::SecondaryDecayRangePositionDistribution second(
        FiducialBox(), shared_daughter, 1.0, 1.0, 20.0);
    EXPECT_TRUE(first == second);
    EXPECT_FALSE(first < second);
    EXPECT_FALSE(second < first);

    // Deep-equal but distinct daughter collections compare by identity:
    // unequal, with a definite relative order.
    siren::distributions::SecondaryDecayRangePositionDistribution third(
        FiducialBox(), Decays(ParticleType::N4, 2.5), 1.0, 1.0, 20.0);
    EXPECT_FALSE(first == third);
    EXPECT_NE(first < third, third < first);
}

TEST(SecondaryDecayRangePositionDistribution,
     SamplesWhenBothInteractionLengthsAreVeryShort) {
    constexpr double parent_length = 0.005;
    constexpr double daughter_length = 0.01;
    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, parent_length);
    auto daughter = Decays(ParticleType::N4, daughter_length);
    siren::distributions::SecondaryDecayRangePositionDistribution distribution(
        FiducialBox(), daughter, 1.0, 1.0, 20.0);

    // Before the fiducial entrance the conditional density falls as
    // exp[-x(1/lambda_parent - 1/lambda_daughter)].  This regime used to
    // underflow the daughter probability and produce a practically useless
    // rejection envelope spanning hundreds of interaction depths.
    InteractionRecord record = ParentRecord();
    record.interaction_vertex = {0.001, 0.0, 0.0};
    double first = distribution.GenerationProbability(detector, parent, record);
    record.interaction_vertex = {0.011, 0.0, 0.0};
    double second = distribution.GenerationProbability(detector, parent, record);
    ASSERT_GT(first, 0.0);
    ASSERT_GT(second, 0.0);
    EXPECT_NEAR(first / second, std::exp(1.0), 2e-4);

    auto random = std::make_shared<siren::utilities::SIREN_random>(8675309);
    InteractionRecord sample_record = ParentRecord();
    siren::dataclasses::SecondaryDistributionRecord secondary(sample_record);
    EXPECT_NO_THROW(distribution.SampleVertex(
        random, detector, parent, secondary));
    EXPECT_LT(secondary.GetLength(), 0.2);
}

TEST(SecondaryDecayRangePositionDistribution, SerializationRoundTrip) {
    std::vector<std::shared_ptr<siren::interactions::Decay>> decays{
        std::make_shared<siren::interactions::HNLDipoleDecay>(
            0.1, 2.5e-6,
            siren::interactions::HNLDipoleDecay::ChiralNature::Dirac)};
    auto daughter = std::make_shared<siren::interactions::InteractionCollection>(
        ParticleType::N4, decays);
    std::shared_ptr<siren::distributions::SecondaryInjectionDistribution> original =
        std::make_shared<
            siren::distributions::SecondaryDecayRangePositionDistribution>(
                FiducialBox(), daughter, 0.1, 0.9, 20.0);

    std::stringstream buffer;
    {
        cereal::JSONOutputArchive archive(buffer);
        archive(cereal::make_nvp("Distribution", original));
    }
    EXPECT_EQ(buffer.str().find("IntegrationSteps"), std::string::npos);

    std::shared_ptr<siren::distributions::SecondaryInjectionDistribution> restored;
    {
        cereal::JSONInputArchive archive(buffer);
        archive(cereal::make_nvp("Distribution", restored));
    }
    ASSERT_NE(restored, nullptr);
    EXPECT_EQ(restored->Name(), "SecondaryDecayRangePositionDistribution");
    auto restored_distribution = std::dynamic_pointer_cast<
        siren::distributions::SecondaryDecayRangePositionDistribution>(restored);
    ASSERT_NE(restored_distribution, nullptr);

    auto detector = UniformDetector();
    auto parent = Decays(ParticleType::NuMu, 12.0);
    InteractionRecord record = ParentRecord();
    record.interaction_vertex = {5.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(
        original->GenerationProbability(detector, parent, record),
        restored_distribution->GenerationProbability(detector, parent, record));
}

int main(int argc, char ** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
