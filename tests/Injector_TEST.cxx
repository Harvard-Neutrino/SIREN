
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Constants.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/Controller.h"

using namespace LeptonInjector;

std::string diff_xs(int Z, int A, std::string mHNL) {
  std::stringstream ss;
  // ss << "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
  ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::string tot_xs(int Z, int A, std::string mHNL) {
  std::stringstream ss;
  // ss << "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/";
  ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/";
    ss << "xsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::vector<std::array<int, 2>> gen_ZA() {
    return std::vector<std::array<int, 2>>{
        {1, 1},
        {6, 12},
        {8, 16},
        {13, 27},
        {14, 28},
        {20, 40},
        {26, 56},
        {29, 63},
        {29, 65},
        {82, 208},
    };
}

std::vector<LeptonInjector::Particle::ParticleType> gen_TargetPIDs() {
    using ParticleType = LeptonInjector::Particle::ParticleType;
    return std::vector<ParticleType>{
        ParticleType::HNucleus,
        ParticleType::C12Nucleus,
        ParticleType::O16Nucleus,
        ParticleType::Al27Nucleus,
        ParticleType::Si28Nucleus,
        ParticleType::Ca40Nucleus,
        ParticleType::Fe56Nucleus,
        ParticleType::Cu63Nucleus,
        ParticleType::Cu65Nucleus,
        ParticleType::Pb208Nucleus
    };
}

std::vector<std::string> gen_diff_xs_hf(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hf(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_diff_xs_hc(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hc(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

TEST(Injector, Generation)
{
    using ParticleType = LeptonInjector::Particle::ParticleType;

    // std::string material_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    // std::string earth_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    std::string material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    std::string earth_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    double powerLawIndex = 2;
    double energyMin = 1; // in GeV
    double energyMax = 1; // in GeV

    double hnl_mass = 0.4; // in GeV; The HNL mass we are injecting
    double d = 1e-7; // in GeV^-1; the effective dipole coupling strength
    std::string mHNL = "0.4";

    // Decay parameters used to set the max range when injecting an HNL, decay width is likely wrong, set accordingly...
    double HNL_decay_width = std::pow(d,2)*std::pow(hnl_mass,3)/(4*Constants::pi); // in GeV; decay_width = d^2 m^3 / (4 * pi)
    double n_decay_lengths = 3.0;

    // This should encompass Minerva, should probably be smaller? Depends on how long Minerva is...
    double disk_radius = 1; // in meters
    double endcap_length = 5; // in meters


    // Events to inject
    unsigned int events_to_inject = 1e3;
    Particle::ParticleType primary_type = ParticleType::NuE;

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = gen_TargetPIDs();
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass, DipoleFromTable::HelicityChannel::Flipping);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass, DipoleFromTable::HelicityChannel::Conserving);
    std::vector<std::string> hf_diff_fnames = gen_diff_xs_hf(mHNL);
    std::vector<std::string> hc_diff_fnames = gen_diff_xs_hc(mHNL);
    std::vector<std::string> hf_tot_fnames = gen_tot_xs_hf(mHNL);
    std::vector<std::string> hc_tot_fnames = gen_tot_xs_hc(mHNL);
    for(unsigned int i=0; i < target_types.size(); ++i) {
        std::cerr << hf_diff_fnames[i] << std::endl;
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        std::cerr << hf_tot_fnames[i] << std::endl;
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        std::cerr << hc_diff_fnames[i] << std::endl;
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        std::cerr << hc_tot_fnames[i] << std::endl;
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    // Load the earth model
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->LoadMaterialModel(material_file);
    earth_model->LoadEarthModel(earth_file);

    // Setup power law
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();
    std::shared_ptr<LeptonInjector::PowerLaw> power_law = std::make_shared<LeptonInjector::PowerLaw>();
    power_law->powerLawIndex = powerLawIndex;
    power_law->energyMin = energyMin;
    power_law->energyMax = energyMax;
    std::shared_ptr<PrimaryEnergyDistribution> edist = power_law;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::FixedDirection>(earthmodel::Vector3D{0.0, 0.0, 1.0});

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to the decay distribution
    std::shared_ptr<RangeFunction> range_func = std::make_shared<LeptonInjector::DecayRangeFunction>(hnl_mass, HNL_decay_width, n_decay_lengths);

    // Spin distribution
    std::shared_ptr<PrimaryNeutrinoSpinDistribution> spin_distribution = std::make_shared<LeptonInjector::PrimaryNeutrinoSpinDistribution>();

    // Put it all together!
    //RangedLeptonInjector injector(events_to_inject, primary_type, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length);
    std::shared_ptr<InjectorBase> injector = std::make_shared<RangedLeptonInjector>(events_to_inject, primary_type, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length, spin_distribution);

    /*
    Controller cont(injector);
    cont.NameOutfile("injector_test_events.h5");
    cont.NameLicFile("injector_test_events.lic");

    // Run the program.
    cont.Execute();
    */

    std::ofstream myFile("injector_test_events.csv");
    myFile << std::fixed << std::setprecision(4);
    myFile << "intX intY intZ ";
    myFile << "decX decY decZ ";
    myFile << "ppX ppY ppZ ";
    myFile << "p4nu_0 p4nu_1 p4nu_2 p4nu_3 ";
    myFile << "snuX snuY snuZ ";
    myFile << "p4itgt_0 p4itgt_1 p4itgt_2 p4itgt_3 ";
    myFile << "sitgtX sitgtY sitgtZ ";
    myFile << "p4hnl_0 p4hnl_1 p4hnl_2 p4hnl_3 ";
    myFile << "shnlX shnlY shnlZ ";
    myFile << "p4ftgt_0 p4ftgt_1 p4ftgt_2 p4ftgt_3 ";
    myFile << "sftgtX sftgtY sftgtZ ";
    myFile << "p4gamma_0 p4gamma_1 p4gamma_2 p4gamma_3 ";
    myFile << "sgammaX sgammaY sgammaZ ";
    myFile << "decay_length prob_nopairprod target\n";
    int i = 0;
    while(*injector) {
        LeptonInjector::InteractionRecord event = injector->GenerateEvent();
        if(event.secondary_momenta.size() > 0)
        {
          myFile << event.interaction_vertex[0] << " ";
          myFile << event.interaction_vertex[1] << " ";
          myFile << event.interaction_vertex[2] << " ";

          myFile << event.decay_vertex[0] << " ";
          myFile << event.decay_vertex[1] << " ";
          myFile << event.decay_vertex[2] << " ";

          myFile << event.pairprod_vertex[0] << " ";
          myFile << event.pairprod_vertex[1] << " ";
          myFile << event.pairprod_vertex[2] << " ";

          myFile << event.primary_momentum[0] << " ";
          myFile << event.primary_momentum[1] << " ";
          myFile << event.primary_momentum[2] << " ";
          myFile << event.primary_momentum[3] << " ";

          myFile << event.primary_spin[0] << " ";
          myFile << event.primary_spin[1] << " ";
          myFile << event.primary_spin[2] << " ";

          myFile << event.target_momentum[0] << " ";
          myFile << event.target_momentum[1] << " ";
          myFile << event.target_momentum[2] << " ";
          myFile << event.target_momentum[3] << " ";

          myFile << event.target_spin[0] << " ";
          myFile << event.target_spin[1] << " ";
          myFile << event.target_spin[2] << " ";

          myFile << event.secondary_momenta[0][0] << " ";
          myFile << event.secondary_momenta[0][1] << " ";
          myFile << event.secondary_momenta[0][2] << " ";
          myFile << event.secondary_momenta[0][3] << " ";

          myFile << event.secondary_spin[0][0] << " ";
          myFile << event.secondary_spin[0][1] << " ";
          myFile << event.secondary_spin[0][2] << " ";

          myFile << event.secondary_momenta[1][0] << " ";
          myFile << event.secondary_momenta[1][1] << " ";
          myFile << event.secondary_momenta[1][2] << " ";
          myFile << event.secondary_momenta[1][3] << " ";

          myFile << event.secondary_spin[1][0] << " ";
          myFile << event.secondary_spin[1][1] << " ";
          myFile << event.secondary_spin[1][2] << " ";

          myFile << event.secondary_momenta[2][0] << " ";
          myFile << event.secondary_momenta[2][1] << " ";
          myFile << event.secondary_momenta[2][2] << " ";
          myFile << event.secondary_momenta[2][3] << " ";

          myFile << event.secondary_spin[2][0] << " ";
          myFile << event.secondary_spin[2][1] << " ";
          myFile << event.secondary_spin[2][2] << " ";

          myFile << event.decay_length << " ";
          myFile << event.prob_nopairprod << " ";
          myFile << event.signature.target_type << "\n";
        }
        std::cout << ++i << std::endl;
    }
    myFile.close();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

