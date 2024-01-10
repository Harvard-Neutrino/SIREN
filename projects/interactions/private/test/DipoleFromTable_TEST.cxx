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
#include "LeptonInjector/interactions/DipoleFromTable.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/dataclasses/InteractionSignature.h"

using namespace LI::interactions;
using namespace LI::dataclasses;
using namespace LI::utilities;

TEST(DipoleFromTable, Constructor)
{
    double hnl_mass = 0.001;
    std::string differential_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/dxsec_Z_6_A_12_mHNL_0.001_hf.dat";
    std::string total_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/xsec_Z_6_A_12_mHNL_0.001_hf.dat";
    std::string proton_differential_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/dxsec_Z_1_A_1_mHNL_0.001_hf.dat";
    std::string proton_total_xs = "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/xsec_Z_1_A_1_mHNL_0.001_hf.dat";
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = {Particle::ParticleType::PPlus, Particle::ParticleType::Neutron, Particle::ParticleType::Nucleon};
    std::shared_ptr<DipoleFromTable> dipole_xs = std::make_shared<DipoleFromTable>(hnl_mass, 1e-7, DipoleFromTable::HelicityChannel::Flipping);
    dipole_xs->AddDifferentialCrossSectionFile(differential_xs, Particle::ParticleType::C12Nucleus);
    dipole_xs->AddTotalCrossSectionFile(total_xs, Particle::ParticleType::C12Nucleus);
    dipole_xs->AddDifferentialCrossSectionFile(proton_differential_xs, Particle::ParticleType::HNucleus);
    dipole_xs->AddTotalCrossSectionFile(proton_total_xs, Particle::ParticleType::HNucleus);
    std::shared_ptr<CrossSection> xs = dipole_xs;

    //std::cerr << "Test cross section" << std::endl << "y    XS" << std::endl;
    for(unsigned int i=0; i<100; ++i) {
        double y = 1e-8 * std::pow(10, i/100.0);
        double test_cross_section = dipole_xs->DifferentialCrossSection(Particle::ParticleType::NuE, 10.0, Particle::ParticleType::C12Nucleus, 12., y);
        //std::cerr << y << " " << test_cross_section << std::endl;
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

