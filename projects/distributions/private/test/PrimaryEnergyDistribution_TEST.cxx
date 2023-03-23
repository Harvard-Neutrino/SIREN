#include <cmath>
#include <tuple>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Integration.h"

#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "LeptonInjector/distributions/primary/energy/Monoenergetic.h"
#include "LeptonInjector/distributions/primary/energy/PowerLaw.h"

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

TEST(PowerLaw, ConstructorMono) {
    double gamma = 0;
    double energyMin = 0;
    double energyMax = 0;
    ASSERT_NO_THROW(PowerLaw A(gamma, energyMin, energyMax););
}

TEST(PowerLaw, ConstructorLog) {
    double gamma = 1;
    double energyMin = 1;
    double energyMax = 2;
    ASSERT_NO_THROW(PowerLaw A(gamma, energyMin, energyMax););
}

TEST(PowerLaw, ConstructorGeneral) {
    double gamma = 2;
    double energyMin = 1;
    double energyMax = 2;
    ASSERT_NO_THROW(PowerLaw A(gamma, energyMin, energyMax););
}

TEST(PowerLaw, SampleBounds) {
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    size_t N = 1000;
    size_t M = 10000;
    for(size_t i=0; i<N; ++i) {
        double gamma = (RandomDouble() - 0.5) + 1;
        double energyMin = RandomDouble();
        double energyRange = RandomDouble();
        double energyMax = energyMin + energyRange;
        PowerLaw A(gamma, energyMin, energyMax);
        for(size_t j=0; j<M; ++j) {
            InteractionRecord record;
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.primary_momentum[0];
            EXPECT_TRUE(test_energy >= energyMin);
            EXPECT_TRUE(test_energy <= energyMax);
        }
    }
}

TEST(PowerLaw, SampleDistribution) {
    size_t N = 1000;
    size_t M = 10000;
    size_t n_one_sigma = 0;
    size_t n_two_sigma = 0;
    size_t n_three_sigma = 0;
    size_t n_four_sigma = 0;
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        double gamma = (RandomDouble() - 0.5) + 1;
        double energyMin = RandomDouble() * 100 + 10;
        double energyRange = RandomDouble() * 100;
        double energyMax = energyMin + energyRange;

        PowerLaw A(gamma, energyMin, energyMax);

        size_t n_bins = (RandomDouble() * M / 500) + 1;
        double bin_min = log(energyMin);
        double bin_max = log(energyMax);
        DistributionTest test(bin_min, bin_max, n_bins);
        std::vector<double> expectation(n_bins);
        std::function<double(double)> f = [&](double y)->double {
            return std::pow(exp(y), -gamma) * exp(y);
        };
        double tol = std::abs(log(std::min(std::pow(energyMin, -gamma), std::pow(energyMax, -gamma)))) * 1e-8;
        double total_integral = 0;
        for(size_t j=0; j<n_bins; ++j) {
            double min = test.bin_edges[j];
            double max = test.bin_edges[j+1];
            double integral = LI::utilities::rombergIntegrate(f, min, max, 1e-8);
            expectation[j] = integral;
            total_integral += integral;
        }
        for(size_t j=0; j<n_bins; ++j) {
            expectation[j] /= total_integral;
            expectation[j] *= M;
        }

        for(size_t j=0; j<M; ++j) {
            LI::dataclasses::InteractionRecord record;
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.primary_momentum[0];
            test.AddValue(log(test_energy));
        }

        if(not test.TestContents(4, expectation)) {
            n_four_sigma += 1;
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(3, expectation)) {
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(2, expectation)) {
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(1, expectation)) {
            n_one_sigma += 1;
        }
    }
    EXPECT_TRUE(double(n_one_sigma) / N <= (1.0 - 0.682689492137086) * (1.0 + sqrt(n_one_sigma)/N));
    EXPECT_TRUE(double(n_two_sigma) / N <= (1.0 - 0.954499736103642) * (1.0 + sqrt(n_two_sigma)/N));
    EXPECT_TRUE(double(n_three_sigma) / N <= (1.0 - 0.997300203936740) * (1.0 + sqrt(n_three_sigma)/N));
    EXPECT_TRUE(double(n_four_sigma) / N <= (1.0 - 0.999936657516334) * (1.0 + sqrt(n_four_sigma)/N));
}

TEST(PowerLaw, GenerationProbability) {
    size_t N = 1000;
    size_t M = 10000;
    std::shared_ptr<LI::utilities::LI_random> rand = std::make_shared<LI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        double gamma = (RandomDouble() - 0.5) + 1;
        double energyMin = RandomDouble() * 100 + 10;
        double energyRange = RandomDouble() * 100;
        double energyMax = energyMin + energyRange;

        PowerLaw A(gamma, energyMin, energyMax);

        std::function<double(double)> f = [&](double y)->double {
            return std::pow(exp(y), -gamma) * exp(y);
        };

        double total_integral = LI::utilities::rombergIntegrate(f, log(energyMin), log(energyMax), 1e-8);

        for(size_t j=0; j<M; ++j) {
            LI::dataclasses::InteractionRecord record;
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.primary_momentum[0];
            double density = A.GenerationProbability(nullptr, nullptr, record);
            double expected_density = std::pow(test_energy, -gamma) / total_integral;
            EXPECT_NEAR(density, expected_density, expected_density * 1e-8);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

