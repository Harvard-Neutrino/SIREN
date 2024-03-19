
#include <cmath>
#include <tuple>
#include <math.h>
#include <random>
#include <iostream>
#include <algorithm>

#include <gtest/gtest.h>

#include "DistributionTest.h"

#include "SIREN/utilities/Random.h"

#include "SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "SIREN/distributions/primary/direction/Cone.h"
#include "SIREN/distributions/primary/direction/FixedDirection.h"
#include "SIREN/distributions/primary/direction/IsotropicDirection.h"

using namespace siren::math;
using namespace siren::distributions;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Cone, Constructor) {
    Vector3D direction(0,0,1);
    double opening_angle = M_PI;
    Cone A(direction, opening_angle);
}

TEST(Cone, SampleBounds) {
    size_t N = 100;
    size_t M = 10000;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        direction.normalize();

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        for(size_t j=0; j<M; ++j) {
            siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
            record.SetEnergy(1);
            record.SetMass(0);
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D vec(record.GetDirection());
            double angle = acos(scalar_product(vec, direction));
            EXPECT_TRUE(angle <= opening_angle);
            EXPECT_TRUE(angle >= 0);
        }
    }
}

TEST(Cone, SampleDistributionTheta) {
    size_t N = 1000;
    size_t M = 10000;
    size_t n_one_sigma = 0;
    size_t n_two_sigma = 0;
    size_t n_three_sigma = 0;
    size_t n_four_sigma = 0;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        direction.normalize();

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        size_t n_bins = (RandomDouble() * M / 500) + 1;
        double bin_max = 1.0;
        double bin_min = cos(opening_angle);
        DistributionTest test(bin_min, bin_max, n_bins);
        for(size_t j=0; j<M; ++j) {
            siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
            record.SetEnergy(1);
            record.SetMass(0);
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D vec(record.GetDirection());
            double c = scalar_product(vec, direction);
            test.AddValue(c);
            double angle = acos(scalar_product(vec, direction));
            EXPECT_TRUE(angle <= opening_angle);
            EXPECT_TRUE(angle >= 0);
        }
        double expected_contents = double(M) / double(n_bins);
        std::vector<double> expect(n_bins, expected_contents);
        if(not test.TestContents(4, expect)) {
            n_four_sigma += 1;
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(3, expect)) {
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(2, expect)) {
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(1, expect)) {
            n_one_sigma += 1;
        }
    }
    EXPECT_TRUE(double(n_one_sigma) / N <= (1.0 - 0.682689492137086) * (1.0 + sqrt(n_one_sigma)/N));
    EXPECT_TRUE(double(n_two_sigma) / N <= (1.0 - 0.954499736103642) * (1.0 + sqrt(n_two_sigma)/N));
    EXPECT_TRUE(double(n_three_sigma) / N <= (1.0 - 0.997300203936740) * (1.0 + sqrt(n_three_sigma)/N));
    EXPECT_TRUE(double(n_four_sigma) / N <= (1.0 - 0.999936657516334) * (1.0 + sqrt(n_four_sigma)/N));
}

TEST(Cone, SampleDistributionPhi) {
    size_t N = 1000;
    size_t M = 10000;
    size_t n_one_sigma = 0;
    size_t n_two_sigma = 0;
    size_t n_three_sigma = 0;
    size_t n_four_sigma = 0;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        direction.normalize();

        Vector3D ortho_1(RandomDouble(), RandomDouble(), RandomDouble());
        while(true) {
            ortho_1.normalize();
            double magnitude = ortho_1.magnitude();
            if(std::isnan(ortho_1.magnitude()) or magnitude == 0)
                ortho_1 = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
            else
                break;
        }
        ortho_1 = cross_product(direction, ortho_1).normalized();
        Vector3D ortho_2 = cross_product(direction, ortho_1);

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        size_t n_bins = (RandomDouble() * M / 500) + 1;
        double bin_max = M_PI;
        double bin_min = -M_PI;
        DistributionTest test(bin_min, bin_max, n_bins);
        for(size_t j=0; j<M; ++j) {
            siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
            record.SetEnergy(1);
            record.SetMass(0);
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D sample_vec(record.GetDirection());
            double phi = std::atan2(scalar_product(sample_vec, ortho_2), scalar_product(sample_vec, ortho_1));
            test.AddValue(phi);
        }
        double expected_contents = double(M) / double(n_bins);
        std::vector<double> expect(n_bins, expected_contents);
        if(not test.TestContents(4, expect)) {
            n_four_sigma += 1;
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(3, expect)) {
            n_three_sigma += 1;
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(2, expect)) {
            n_two_sigma += 1;
            n_one_sigma += 1;
        } else if(not test.TestContents(1, expect)) {
            n_one_sigma += 1;
        }
    }
    EXPECT_TRUE(double(n_one_sigma) / N <= (1.0 - 0.682689492137086) * (1.0 + sqrt(n_one_sigma)/N));
    EXPECT_TRUE(double(n_two_sigma) / N <= (1.0 - 0.954499736103642) * (1.0 + sqrt(n_two_sigma)/N));
    EXPECT_TRUE(double(n_three_sigma) / N <= (1.0 - 0.997300203936740) * (1.0 + sqrt(n_three_sigma)/N));
    EXPECT_TRUE(double(n_four_sigma) / N <= (1.0 - 0.999936657516334) * (1.0 + sqrt(n_four_sigma)/N));
}

TEST(Cone, GenerationProbability) {
    size_t N = 1000;
    size_t M = 10000;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        while(true) {
            direction.normalize();
            double magnitude = direction.magnitude();
            if(std::isnan(direction.magnitude()) or magnitude == 0)
                direction = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
            else
                break;
        }

        double opening_angle = M_PI * RandomDouble();
        Cone A(direction, opening_angle);
        double expected_density = 1.0 / ((2.0 * M_PI) * (1.0 - cos(opening_angle)));
        for(size_t j=0; j<M; ++j) {
            double input_angle = RandomDouble() * opening_angle;
            Vector3D vec(RandomDouble(), RandomDouble(), RandomDouble());
            while(true) {
                vec.normalize();
                double magnitude = vec.magnitude();
                if(std::isnan(vec.magnitude()) or magnitude == 0)
                    vec = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
                else
                    break;
            }
            vec = cross_product(direction, vec).normalized();
            Quaternion q(vec * sin(input_angle));
            q.SetW(1.0 + cos(input_angle));
            q.normalize();
            vec = q.rotate(direction, false);
            vec.normalize();

            double pre_c = scalar_product(vec, direction);
            if(pre_c > 1)
                pre_c = 1;
            EXPECT_TRUE(acos(pre_c) <= opening_angle);

            siren::dataclasses::InteractionRecord record;
            record.primary_momentum[1] = vec.GetX();
            record.primary_momentum[2] = vec.GetY();
            record.primary_momentum[3] = vec.GetZ();

            double c = scalar_product(direction.normalized(), vec.normalized());
            double output_angle;
            if(c > 1)
                output_angle = 0;
            else
                output_angle = acos(c);

            double resolution;
            if(input_angle == 0)
                resolution = 0;
            else
                resolution = 1.0 / sqrt(input_angle) * 1e-8;

            EXPECT_NEAR(input_angle, output_angle, resolution);

            double density = A.GenerationProbability(nullptr, nullptr, record);

            if(output_angle < opening_angle)
                EXPECT_NEAR(density, expected_density, expected_density * 1e-8);
            else
                EXPECT_TRUE(density == 0);
        }
    }
}

TEST(FixedDirection, Constructor) {
    Vector3D direction(0,0,1);
    FixedDirection A(direction);
}

TEST(FixedDirection, Sample) {
    size_t N = 10000;
    size_t M = 10;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    for(size_t i=0; i<N; ++i) {
        Vector3D direction(RandomDouble(), RandomDouble(), RandomDouble());
        while(true) {
            direction.normalize();
            double magnitude = direction.magnitude();
            if(std::isnan(direction.magnitude()) or magnitude == 0)
                direction = Vector3D(RandomDouble(), RandomDouble(), RandomDouble());
            else
                break;
        }
        FixedDirection A(direction);
        for(size_t j=0; j<M; ++j) {
            siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
            record.SetEnergy(1);
            record.SetMass(0);
            A.Sample(rand, nullptr, nullptr, record);
            Vector3D vec(record.GetDirection());
            double dot = scalar_product(direction, vec);
            EXPECT_TRUE(std::abs(dot - 1.0) < 1e-9);
        }
    }
}

TEST(IsotropicDirection, Constructor) {
    IsotropicDirection A;
}

TEST(IsotropicDirection, SampleDistributionTheta) {
    size_t N = 100000;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    IsotropicDirection A;
    double bin_max = 1.0;
    double bin_min = -1.0;
    Vector3D direction(0,0,1);
    size_t n_bins = (N / 5000) + 1;
    DistributionTest test(bin_min, bin_max, n_bins);
    for(size_t i=0; i<N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        A.Sample(rand, nullptr, nullptr, record);
        Vector3D vec(record.GetDirection());
        double c = scalar_product(vec, direction);
        test.AddValue(c);
        double angle = acos(c);
        EXPECT_TRUE(angle <= M_PI);
        EXPECT_TRUE(angle >= 0);
    }
    double expected_contents = double(N) / double(n_bins);
    std::vector<double> expect(n_bins, expected_contents);
    EXPECT_TRUE(test.TestContents(3, expect));
}

TEST(IsotropicDirection, SampleDistributionPhi) {
    size_t N = 100000;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    IsotropicDirection A;
    double bin_max = M_PI;
    double bin_min = -M_PI;
    Vector3D direction(0,0,1);
    Vector3D ortho_1(1,0,0);
    Vector3D ortho_2(0,1,0);
    size_t n_bins = (N / 5000) + 1;
    DistributionTest test(bin_min, bin_max, n_bins);
    for(size_t i=0; i<N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        A.Sample(rand, nullptr, nullptr, record);
        Vector3D vec(record.GetDirection());
        double phi = atan2(scalar_product(ortho_2, vec), scalar_product(ortho_1, vec));
        test.AddValue(phi);
        EXPECT_TRUE(phi <= M_PI);
        EXPECT_TRUE(phi >= -M_PI);
    }
    double expected_contents = double(N) / double(n_bins);
    std::vector<double> expect(n_bins, expected_contents);
    EXPECT_TRUE(test.TestContents(3, expect));
}

TEST(IsotropicDirection, GenerationProbability) {
    size_t N = 100000;
    std::shared_ptr<siren::utilities::SIREN_random> rand = std::make_shared<siren::utilities::SIREN_random>();
    IsotropicDirection A;
    double expected_density = 1.0 / (4.0 * M_PI);
    for(size_t i=0; i<N; ++i) {
        siren::dataclasses::PrimaryDistributionRecord record(siren::dataclasses::ParticleType::NuMu);
        record.SetEnergy(1);
        record.SetMass(0);
        A.Sample(rand, nullptr, nullptr, record);
        Vector3D vec(record.GetDirection());
        siren::dataclasses::InteractionRecord interaction_record;
        record.FinalizeAvailable(interaction_record);
        double density = A.GenerationProbability(nullptr, nullptr, interaction_record);
        EXPECT_NEAR(density, expected_density, expected_density * 1e-8);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

