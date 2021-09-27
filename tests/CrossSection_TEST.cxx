
#include <cmath>
#include <math.h>
#include <memory>
#include <iostream>

#include <gtest/gtest.h>
#include <cereal/archives/binary.hpp>

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"

using namespace LeptonInjector;

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
    event.target_momentum[0] = target_mass;
    event.target_momentum[1] = 0.1;

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
}

TEST(DipoleFromTable, Constructor)
{
    double hnl_mass = 0.001;
    std::string differential_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/dxsec_Z_82_A_208_mHNL_0.001_hf.dat";
    std::string total_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/xsec_Z_82_A_208_mHNL_0.001_hf.dat";
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = {Particle::ParticleType::PPlus, Particle::ParticleType::Neutron, Particle::ParticleType::Nucleon};
    std::shared_ptr<DipoleFromTable> dipole_xs = std::make_shared<DipoleFromTable>(hnl_mass);
    dipole_xs->AddDifferentialCrossSectionFile(differential_xs, Particle::ParticleType::Pb208Nucleus);
    dipole_xs->AddTotalCrossSectionFile(total_xs, Particle::ParticleType::Pb208Nucleus);
    std::shared_ptr<CrossSection> xs = dipole_xs;

    InteractionSignature signature;
    signature.primary_type = Particle::ParticleType::NuE;
    signature.target_type = Particle::ParticleType::Pb208Nucleus;
    signature.secondary_types = {Particle::ParticleType::NuF4, Particle::ParticleType::Pb208Nucleus};
    InteractionRecord event;
    event.signature = signature;
    double energy = 10; // 10GeV
    event.primary_momentum[0] = energy; // 10GeV
    double target_mass = 976652.005 * 0.9314941 * 1e-6;
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

    for(unsigned int i=0; i<10; ++i) {
        xs->SampleFinalState(event, rand);
        std::cerr << event << std::endl;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

