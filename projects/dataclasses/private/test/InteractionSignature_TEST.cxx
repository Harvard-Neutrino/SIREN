
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

#include "SIREN/dataclasses/InteractionSignature.h"

using namespace siren::dataclasses;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Comparison, Comparison_equal)
{
    InteractionSignature A;
    InteractionSignature B;
    EXPECT_TRUE(A == B);

    A.primary_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_FALSE(A == B);
    B.primary_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_TRUE(A == B);

    A.target_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_FALSE(A == B);
    B.target_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_TRUE(A == B);

    A.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    EXPECT_FALSE(A == B);
    B.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    EXPECT_TRUE(A == B);

    A.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);
    EXPECT_FALSE(A == B);
    B.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);
    EXPECT_TRUE(A == B);

    std::swap(A.secondary_types[0], A.secondary_types[1]);
    EXPECT_FALSE(A == B);
    std::swap(A.secondary_types[0], A.secondary_types[1]);
    EXPECT_TRUE(A == B);
}

TEST(Comparison, LessThan)
{
    InteractionSignature A;
    InteractionSignature B;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.primary_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_TRUE(A < B);
    EXPECT_FALSE(B < A);

    B.primary_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.target_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_TRUE(A < B);
    EXPECT_FALSE(B < A);

    B.target_type = siren::dataclasses::ParticleType::EPlus;
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    //
    InteractionSignature C;
    InteractionSignature D;
    EXPECT_FALSE(C < D);
    EXPECT_FALSE(D < C);

    C.primary_type = siren::dataclasses::ParticleType::EMinus;
    EXPECT_FALSE(C < D);
    EXPECT_TRUE(D < C);

    D.primary_type = siren::dataclasses::ParticleType::EMinus;
    EXPECT_FALSE(C < D);
    EXPECT_FALSE(D < C);

    C.target_type = siren::dataclasses::ParticleType::EMinus;
    EXPECT_FALSE(C < D);
    EXPECT_TRUE(D < C);

    D.target_type = siren::dataclasses::ParticleType::EMinus;
    EXPECT_FALSE(C < D);
    EXPECT_FALSE(D < C);
    //

    A.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    A.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    B.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);

    std::swap(A.secondary_types[0], A.secondary_types[1]);
    EXPECT_FALSE(A < B);
    EXPECT_TRUE(B < A);

    std::swap(A.secondary_types[0], A.secondary_types[1]);
    EXPECT_FALSE(A < B);
    EXPECT_FALSE(B < A);
}

TEST(Serialization, Save)
{
    InteractionSignature A;
    A.primary_type = siren::dataclasses::ParticleType::EPlus;
    A.target_type = siren::dataclasses::ParticleType::EPlus;
    A.secondary_types.push_back(siren::dataclasses::ParticleType::EPlus);
    A.secondary_types.push_back(siren::dataclasses::ParticleType::EMinus);

    std::stringstream ss;
    {
        cereal::JSONOutputArchive oarchive(ss);
        oarchive(cereal::make_nvp("InteractionSignature", A));
    }

    std::string expected = "{\n\
    \"InteractionSignature\": {\n\
        \"cereal_class_version\": 0,\n\
        \"PrimaryType\": -11,\n\
        \"TargetType\": -11,\n\
        \"SecondaryTypes\": [\n\
            -11,\n\
            11\n\
        ]\n\
    }\n\
}";
    EXPECT_EQ(ss.str(), expected);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

