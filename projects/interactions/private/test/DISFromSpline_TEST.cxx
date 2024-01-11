#include <cmath>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>

#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/DISFromSpline.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"

using namespace LI::interactions;
using namespace LI::dataclasses;
using namespace LI::utilities;

TEST(DISFromSpline, Constructor)
{
    std::string differential_xs = "/home/austin/programs/LIDUNE/sources/DUNEAtmo/cross_sections/csms_differential_v1.0/dsdxdy_nu_CC_iso.fits";
    std::string total_xs = "/home/austin/programs/LIDUNE/sources/DUNEAtmo/cross_sections/csms_differential_v1.0/sigma_nu_CC_iso.fits";
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = {Particle::ParticleType::PPlus, Particle::ParticleType::Neutron, Particle::ParticleType::Nucleon};
    std::shared_ptr<DISFromSpline> dis_xs = std::make_shared<DISFromSpline>(differential_xs, total_xs, primary_types, target_types);
    std::shared_ptr<CrossSection> xs = dis_xs;

    InteractionSignature signature;
    signature.primary_type = Particle::ParticleType::NuE;
    signature.target_type = Particle::ParticleType::Nucleon;
    signature.secondary_types = {Particle::ParticleType::EMinus, Particle::ParticleType::Hadrons};
    InteractionRecord event;
    event.signature = signature;
    double energy = 1e4; // 10TeV
    event.primary_momentum[0] = energy; // 10TeV
    double target_mass = dis_xs->GetTargetMass();
    event.target_mass = target_mass;
    event.target_momentum[0] = target_mass;

    std::shared_ptr<LI_random> rand = std::make_shared<LI_random>();

    double z = rand->Uniform(-1.0, 1.0);
    double rho = sqrt(1.0 - z*z);
    double theta = rand->Uniform(0.0, 2.0*M_PI);
    double x = rho * cos(theta);
    double y = rho * sin(theta);

    event.primary_momentum[1] = x * energy;
    event.primary_momentum[2] = y * energy;
    event.primary_momentum[3] = z * energy;

    xs->SampleFinalState(event, rand);
    cereal::JSONOutputArchive output(std::cout);
    //output(xs);
    output(event);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

