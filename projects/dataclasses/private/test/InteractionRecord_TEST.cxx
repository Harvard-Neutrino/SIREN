
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <gtest/gtest.h>

#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"

using namespace LI::dataclasses;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Comparison, Comparison_equal)
{
    InteractionRecord A;
    InteractionRecord B;
    EXPECT_TRUE(A == B);

    A.signature.primary_type = Particle::ParticleType::EMinus;
    EXPECT_FALSE(A == B);
    B.signature.primary_type = Particle::ParticleType::EMinus;
    EXPECT_TRUE(A == B);

    A.signature.target_type = Particle::ParticleType::PPlus;
    EXPECT_FALSE(A == B);
    B.signature.target_type = Particle::ParticleType::PPlus;
    EXPECT_TRUE(A == B);

    A.signature.secondary_types.push_back(Particle::ParticleType::EMinus);
    EXPECT_FALSE(A == B);
    B.signature.secondary_types.push_back(Particle::ParticleType::EMinus);
    EXPECT_TRUE(A == B);

    A.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A == B);
    B.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_TRUE(A == B);

    A.primary_mass = 1.0;
    EXPECT_FALSE(A == B);
    B.primary_mass = 1.0;
    EXPECT_TRUE(A == B);

    A.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A == B);
    B.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_TRUE(A == B);

    A.primary_helicity = 1.0;
    EXPECT_FALSE(A == B);
    B.primary_helicity = 1.0;
    EXPECT_TRUE(A == B);

    A.target_mass = 1.0;
    EXPECT_FALSE(A == B);
    B.target_mass = 1.0;
    EXPECT_TRUE(A == B);

    A.target_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A == B);
    B.target_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_TRUE(A == B);

    A.target_helicity = 1.0;
    EXPECT_FALSE(A == B);
    B.target_helicity = 1.0;
    EXPECT_TRUE(A == B);

    A.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A == B);
    B.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_TRUE(A == B);

    A.secondary_masses.push_back(1.0);
    EXPECT_FALSE(A == B);
    B.secondary_masses.push_back(1.0);
    EXPECT_TRUE(A == B);

    A.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_FALSE(A == B);
    B.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_TRUE(A == B);

    A.secondary_helicity.push_back(1.0);
    EXPECT_FALSE(A == B);
    B.secondary_helicity.push_back(1.0);
    EXPECT_TRUE(A == B);

    A.interaction_parameters["test"] = 1.0;
    EXPECT_FALSE(A == B);
    B.interaction_parameters["test"] = 1.0;
    EXPECT_TRUE(A == B);

    A.interaction_parameters["test"] = 2.0;
    EXPECT_FALSE(A == B);
    B.interaction_parameters["test"] = 2.0;
    EXPECT_TRUE(A == B);

}

TEST(Comparison, LessThan)
{
    InteractionRecord A;
    InteractionRecord B;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.signature.primary_type = Particle::ParticleType::EMinus;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.signature.primary_type = Particle::ParticleType::EMinus;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.signature.target_type = Particle::ParticleType::PPlus;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.signature.target_type = Particle::ParticleType::PPlus;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.signature.secondary_types.push_back(Particle::ParticleType::EMinus);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.signature.secondary_types.push_back(Particle::ParticleType::EMinus);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.primary_initial_position = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.primary_mass = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.primary_mass = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.primary_helicity = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.primary_helicity = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.target_mass = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.target_mass = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.target_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.target_momentum = {1.0, 2.0, 3.0, 4.0};
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.target_helicity = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.target_helicity = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.interaction_vertex = {1.0, 2.0, 3.0};
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.secondary_masses.push_back(1.0);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.secondary_masses.push_back(1.0);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);
    
    B.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.secondary_helicity.push_back(1.0);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.secondary_helicity.push_back(1.0);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.interaction_parameters["test"] = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.interaction_parameters["test"] = 1.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.interaction_parameters["test"] = 2.0;
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.interaction_parameters["test"] = 2.0;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);
}

TEST(Serialization, Save)
{
    InteractionRecord record;
    record.signature.primary_type = Particle::ParticleType::EMinus;
    record.signature.target_type = Particle::ParticleType::PPlus;
    record.signature.secondary_types.push_back(Particle::ParticleType::EMinus);
    record.primary_mass = 1.0;
    record.primary_momentum = {1.0, 2.0, 3.0, 4.0};
    record.primary_helicity = 1.0;
    record.target_mass = 1.0;
    record.target_momentum = {1.0, 2.0, 3.0, 4.0};
    record.target_helicity = 1.0;

    record.interaction_vertex = {1.0, 2.0, 3.0};
    record.secondary_masses.push_back(1.0);
    record.secondary_momenta.push_back({1.0, 2.0, 3.0, 4.0});
    record.secondary_helicity.push_back(1.0);
    record.interaction_parameters["test"] = 1.0;

    std::stringstream ss;
    {
        cereal::JSONOutputArchive oarchive(ss);
        oarchive(record);
    }

    std::string expected = "{\n    \"cereal_class_version\": 0,\n    \"InteractionSignature\": {\n        \"cereal_class_version\": 0,\n        \"PrimaryType\": -11,\n        \"TargetType\": 2212,\n        \"SecondaryTypes\": [\n            -11\n        ]\n    },\n    \"PrimaryInitialPosition\": [\n        1,\n        2,\n        3\n    ],\n    \"PrimaryMass\": 1,\n    \"PrimaryMomentum\": [\n        1,\n        2,\n        3,\n        4\n    ],\n    \"PrimaryHelicity\": 1,\n    \"TargetMass\": 1,\n    \"TargetMomentum\": [\n        1,\n        2,\n        3,\n        4\n    ],\n    \"TargetHelicity\": 1,\n    \"InteractionVertex\": [\n        1,\n        2,\n        3\n    ],\n    \"SecondaryMasses\": [\n        1\n    ],\n    \"SecondaryMomenta\": [\n        [\n            1,\n            2,\n            3,\n            4\n        ]\n    ],\n    \"SecondaryHelicity\": [\n        1\n    ],\n    \"InteractionParameters\": {\n        \"test\": 1\n    }\n}";
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

