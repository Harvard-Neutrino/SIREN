#include <cmath>
#include <tuple>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "LeptonInjector/distributions/primary/energy/Monoenergetic.h"

using namespace LI::distributions;
using namespace LI::dataclasses;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}


TEST(Monoenergetic, Constructor) {
    double energy = 0;
    ASSERT_NO_THROW(Monoenergetic A(energy););
}

TEST(Monoenergetic, Sample) {
    size_t N = 1000;
    size_t M = 10000;
    for(size_t i=0; i<N; ++i) {
        double energy = RandomDouble();
        Monoenergetic A(energy);
        for(size_t j=0; j<M; ++j) {
            InteractionRecord record;
            A.Sample(nullptr, nullptr, nullptr, record);
            double test_energy = record.primary_momentum[0];
            EXPECT_NEAR(energy, test_energy, energy * 1e-6);
        }
    }
}

TEST(Monoenergetic, GenerationProbability) {
    size_t N = 1000;
    size_t M = 10000;
    for(size_t i=0; i<N; ++i) {
        double energy = RandomDouble();
        Monoenergetic A(energy);
        for(size_t j=0; j<M; ++j) {
            InteractionRecord record;
            double test_energy = RandomDouble();
            record.primary_momentum[0] = test_energy;
            double density = A.GenerationProbability(nullptr, nullptr, record);
            if(std::abs(energy - test_energy) / energy < 1e-6) {
                EXPECT_TRUE(density == 1);
            } else {
                EXPECT_TRUE(density == 0);
            }
            test_energy = energy;
            record.primary_momentum[0] = test_energy;
            density = A.GenerationProbability(nullptr, nullptr, record);
            EXPECT_TRUE(density == 1);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

