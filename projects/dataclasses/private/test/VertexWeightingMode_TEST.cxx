#include <sstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>

#include "SIREN/dataclasses/VertexWeightingMode.h"

using siren::dataclasses::VertexWeightingMode;
using siren::dataclasses::VertexWeightingModeName;

namespace {

VertexWeightingMode RoundTrip(VertexWeightingMode const & mode) {
    std::stringstream stream;
    {
        cereal::JSONOutputArchive output(stream);
        output(cereal::make_nvp("Mode", mode));
    }
    VertexWeightingMode loaded = VertexWeightingMode::Fixed();
    {
        cereal::JSONInputArchive input(stream);
        input(cereal::make_nvp("Mode", loaded));
    }
    return loaded;
}

VertexWeightingMode FromJSON(std::string const & json) {
    std::stringstream stream(json);
    cereal::JSONInputArchive input(stream);
    VertexWeightingMode loaded;
    input(cereal::make_nvp("Mode", loaded));
    return loaded;
}

} // namespace

TEST(VertexWeightingMode, PresetsRoundTripThroughVersion1) {
    for(auto mode : {VertexWeightingMode::Propagated(), VertexWeightingMode::Fixed(),
                     VertexWeightingMode::ExternalBounds(),
                     VertexWeightingMode::PropagatedFromCreation()}) {
        EXPECT_EQ(RoundTrip(mode), mode) << VertexWeightingModeName(mode);
    }
    EXPECT_TRUE(RoundTrip(VertexWeightingMode::PropagatedFromCreation()).survival_from_creation);
}

TEST(VertexWeightingMode, PropagatedFromCreationDiffersOnlyInSurvival) {
    VertexWeightingMode mode = VertexWeightingMode::PropagatedFromCreation();
    EXPECT_NE(mode, VertexWeightingMode::Propagated());
    mode.survival_from_creation = false;
    EXPECT_EQ(mode, VertexWeightingMode::Propagated());
    EXPECT_EQ(VertexWeightingModeName(VertexWeightingMode::PropagatedFromCreation()), "PropagatedFromCreation");
    EXPECT_EQ(VertexWeightingModeName(VertexWeightingMode::Propagated()), "Propagated");
}

TEST(VertexWeightingMode, Version0ArchiveLoadsWithoutSurvival) {
    VertexWeightingMode loaded = FromJSON(R"JSON({"Mode": {
        "cereal_class_version": 0,
        "ComputeInteractionProbability": true,
        "ComputePositionProbability": true,
        "BoundSource": 0}})JSON");
    EXPECT_EQ(loaded, VertexWeightingMode::Propagated());
    EXPECT_FALSE(loaded.survival_from_creation);
}

TEST(VertexWeightingMode, InconsistentArchiveIsRejected) {
    EXPECT_THROW(FromJSON(R"JSON({"Mode": {
        "cereal_class_version": 1,
        "ComputeInteractionProbability": false,
        "ComputePositionProbability": false,
        "BoundSource": 2,
        "SurvivalFromCreation": true}})JSON"), std::runtime_error);
}

TEST(VertexWeightingMode, NewerArchiveIsRejected) {
    EXPECT_THROW(FromJSON(R"JSON({"Mode": {
        "cereal_class_version": 2,
        "ComputeInteractionProbability": true,
        "ComputePositionProbability": true,
        "BoundSource": 0,
        "SurvivalFromCreation": false}})JSON"), std::runtime_error);
}

TEST(VertexWeightingMode, SurvivalRequiresBothFactorsAndGeometryBounds) {
    EXPECT_NO_THROW(VertexWeightingMode::PropagatedFromCreation().Validate());
    for(auto mode : {VertexWeightingMode::Fixed(), VertexWeightingMode::ExternalBounds()}) {
        EXPECT_NO_THROW(mode.Validate());
        mode.survival_from_creation = true;
        EXPECT_THROW(mode.Validate(), std::runtime_error) << VertexWeightingModeName(mode);
    }
    VertexWeightingMode partial = VertexWeightingMode::PropagatedFromCreation();
    partial.compute_position_probability = false;
    EXPECT_THROW(partial.Validate(), std::runtime_error);
}
