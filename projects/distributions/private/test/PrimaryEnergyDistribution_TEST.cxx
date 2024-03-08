#include <cmath>
#include <tuple>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Integration.h"

#include "SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "SIREN/distributions/primary/energy/Monoenergetic.h"
#include "SIREN/distributions/primary/energy/PowerLaw.h"
#include "SIREN/distributions/primary/energy/ModifiedMoyalPlusExponentialEnergyDistribution.h"

using namespace SI::distributions;
using namespace SI::dataclasses;

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
            SI::dataclasses::PrimaryDistributionRecord record(SI::dataclasses::ParticleType::NuMu);
            A.Sample(nullptr, nullptr, nullptr, record);
            double test_energy = record.GetEnergy();
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
    std::shared_ptr<SI::utilities::LI_random> rand = std::make_shared<SI::utilities::LI_random>();
    size_t N = 1000;
    size_t M = 10000;
    for(size_t i=0; i<N; ++i) {
        double gamma = (RandomDouble() - 0.5) + 1;
        double energyMin = RandomDouble();
        double energyRange = RandomDouble();
        double energyMax = energyMin + energyRange;
        PowerLaw A(gamma, energyMin, energyMax);
        for(size_t j=0; j<M; ++j) {
            SI::dataclasses::PrimaryDistributionRecord record(SI::dataclasses::ParticleType::NuMu);
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.GetEnergy();
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
    std::shared_ptr<SI::utilities::LI_random> rand = std::make_shared<SI::utilities::LI_random>();
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
        double total_integral = 0;
        for(size_t j=0; j<n_bins; ++j) {
            double min = test.bin_edges[j];
            double max = test.bin_edges[j+1];
            double integral = SI::utilities::rombergIntegrate(f, min, max, 1e-8);
            expectation[j] = integral;
            total_integral += integral;
        }
        for(size_t j=0; j<n_bins; ++j) {
            expectation[j] /= total_integral;
            expectation[j] *= M;
        }

        for(size_t j=0; j<M; ++j) {
            SI::dataclasses::PrimaryDistributionRecord record(SI::dataclasses::ParticleType::NuMu);
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.GetEnergy();
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
    std::shared_ptr<SI::utilities::LI_random> rand = std::make_shared<SI::utilities::LI_random>();
    for(size_t i=0; i<N; ++i) {
        double gamma = (RandomDouble() - 0.5) + 1;
        double energyMin = RandomDouble() * 100 + 10;
        double energyRange = RandomDouble() * 100;
        double energyMax = energyMin + energyRange;

        PowerLaw A(gamma, energyMin, energyMax);

        std::function<double(double)> f = [&](double y)->double {
            return std::pow(exp(y), -gamma) * exp(y);
        };

        double total_integral = SI::utilities::rombergIntegrate(f, log(energyMin), log(energyMax), 1e-8);

        for(size_t j=0; j<M; ++j) {
            SI::dataclasses::PrimaryDistributionRecord record(SI::dataclasses::ParticleType::NuMu);
            A.Sample(rand, nullptr, nullptr, record);
            double test_energy = record.GetEnergy();
            SI::dataclasses::InteractionRecord interaction_record;
            record.FinalizeAvailable(interaction_record);
            double density = A.GenerationProbability(nullptr, nullptr, interaction_record);
            double expected_density = std::pow(test_energy, -gamma) / total_integral;
            EXPECT_NEAR(density, expected_density, expected_density * 1e-8);
        }
    }
}

TEST(ModifiedMoyalPlusExponentialEnergyDistribution, ConstructorNonPhysical) {
    size_t N = 1000;
    double energyMin;
    double energyMax;
    double mu;
    double sigma;
    double A;
    double l;
    double B;

    for(size_t i=0; i<N; ++i) {
        energyMin = RandomDouble() * 10;
        energyMax = energyMin + RandomDouble() * 10;
        mu = energyMin + RandomDouble() * (energyMax - energyMin);
        sigma = RandomDouble();
        A = RandomDouble() * 0.1;
        l = RandomDouble() * 5;
        B = RandomDouble() * 0.1;

        ModifiedMoyalPlusExponentialEnergyDistribution dist(energyMin, energyMax, mu, sigma, A, l, B);
        bool has_norm = dist.IsNormalizationSet();
        EXPECT_TRUE(not has_norm);
    }
}

TEST(ModifiedMoyalPlusExponentialEnergyDistribution, ConstructorPhysical) {
    size_t N = 1000;
    double energyMin;
    double energyMax;
    double mu;
    double sigma;
    double A;
    double l;
    double B;
    bool has_physical_normalization = true;

    for(size_t i=0; i<N; ++i) {
        energyMin = RandomDouble() * 10;
        energyMax = energyMin + RandomDouble() * 10;
        mu = energyMin + RandomDouble() * (energyMax - energyMin);
        sigma = RandomDouble();
        A = RandomDouble() * 0.1;
        l = RandomDouble() * 5;
        B = RandomDouble() * 0.1;

        ModifiedMoyalPlusExponentialEnergyDistribution dist(energyMin, energyMax, mu, sigma, A, l, B, has_physical_normalization);
        bool has_norm = dist.IsNormalizationSet();
        EXPECT_TRUE(has_norm);
    }
}

TEST(ModifiedMoyalPlusExponentialEnergyDistribution, Normalization) {
    size_t N = 1000;
    double energyMin;
    double energyMax;
    double mu;
    double sigma;
    double A;
    double l;
    double B;
    bool has_physical_normalization = true;

    std::function<double(double)> unnormed_pdf = [&](double x)->double {
        double z = (x - mu) / sigma;
        double exponential = B * exp(-x/l) / l;
        double moyal = A * exp(0.5 * -(exp(-z) + z)) / (2.0 * M_PI * sigma);
        return exponential + moyal;
    };

    std::function<double()> normalization = [&]()->double {
        double exponential = B * (exp(-energyMin/l) - exp(-energyMax/l));
        double moyal = A * (std::erf(exp((mu - energyMin)/(2.0 * sigma)) / sqrt(2.0)) - std::erf(exp((mu - energyMax)/(2.0 * sigma)) / sqrt(2.0)));
        return exponential + moyal;
    };


    for(size_t i=0; i<N; ++i) {
        energyMin = RandomDouble() * 10 + 1e-3;
        energyMax = energyMin + RandomDouble() * 10;
        mu = energyMin + RandomDouble() * (energyMax - energyMin);
        sigma = (RandomDouble()*10 + 1e-2);
        A = RandomDouble() * 0.1 + 1e-3;
        l = RandomDouble() * 5 + 1e-3;
        B = RandomDouble() * 0.1 + 1e-3;

        ModifiedMoyalPlusExponentialEnergyDistribution dist(energyMin, energyMax, mu, sigma, A, l, B, has_physical_normalization);

        std::function<double(double)> pdf = [&](double y)->double {
            InteractionRecord record;
            record.primary_momentum[0] = exp(y);
            return dist.GenerationProbability(nullptr, nullptr, record) * exp(y);
        };

        double norm = SI::utilities::rombergIntegrate(pdf, log(energyMin), log(energyMax), 1e-6);
        EXPECT_NEAR(norm, 1.0, 2e-3);

        double expected_norm = dist.GetNormalization();
        EXPECT_NEAR(expected_norm, normalization(), normalization() * 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

