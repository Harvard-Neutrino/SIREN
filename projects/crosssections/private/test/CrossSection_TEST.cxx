#include <cmath>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>

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

TEST(DipoleFromTable, Constructor)
{
    double hnl_mass = 0.001;
    std::string differential_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/dxsec_Z_6_A_12_mHNL_0.001_hf.dat";
    std::string total_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/xsec_Z_6_A_12_mHNL_0.001_hf.dat";
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = {Particle::ParticleType::PPlus, Particle::ParticleType::Neutron, Particle::ParticleType::Nucleon};
    std::shared_ptr<DipoleFromTable> dipole_xs = std::make_shared<DipoleFromTable>(hnl_mass, 1e-7, DipoleFromTable::HelicityChannel::Flipping);
    dipole_xs->AddDifferentialCrossSectionFile(differential_xs, Particle::ParticleType::C12Nucleus);
    dipole_xs->AddTotalCrossSectionFile(total_xs, Particle::ParticleType::C12Nucleus);
    std::shared_ptr<CrossSection> xs = dipole_xs;

    std::cerr << "Test cross section" << std::endl << "y    XS" << std::endl;
    for(unsigned int i=0; i<100; ++i) {
        double y = 1e-8 * std::pow(10, i/100.0);
        double test_cross_section = dipole_xs->DifferentialCrossSection(Particle::ParticleType::NuE, 10.0, Particle::ParticleType::C12Nucleus, 12., y);
        std::cerr << y << " " << test_cross_section << std::endl;
    }
    // return;

    InteractionSignature signature;
    signature.primary_type = Particle::ParticleType::NuE;
    signature.target_type = Particle::ParticleType::C12Nucleus;
    signature.secondary_types = {Particle::ParticleType::NuF4, Particle::ParticleType::C12Nucleus};
    InteractionRecord event;
    event.signature = signature;
    double energy = 10; // 10GeV
    event.primary_momentum[0] = energy; // 10GeV
    double carbon_mass = 12000000.0;
    double target_mass = carbon_mass * 0.9314941 * 1e-6;
    event.target_mass = target_mass;
    event.target_momentum[0] = target_mass;
    event.target_mass = target_mass;

    std::shared_ptr<LI_random> rand = std::make_shared<LI_random>();

    double z = rand->Uniform(-1.0, 1.0);
    double rho = sqrt(1.0 - z*z);
    double theta = rand->Uniform(0.0, 2.0*M_PI);
    double x = rho * cos(theta);
    double y = rho * sin(theta);

    event.primary_momentum[1] = x * energy;
    event.primary_momentum[2] = y * energy;
    event.primary_momentum[3] = z * energy;

    cereal::JSONOutputArchive err(std::cerr);
    err(*dipole_xs);

    std::ofstream out;
    out.open("carbon_test.json");

    cereal::JSONOutputArchive output(out);
    unsigned int total_events = 10000;
    output(cereal::make_size_tag(static_cast<size_t>(total_events)));
    for(unsigned int i=0; i<total_events; ++i) {
        xs->SampleFinalState(event, rand);
        //std::cerr << event << std::endl;
        output(event);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

