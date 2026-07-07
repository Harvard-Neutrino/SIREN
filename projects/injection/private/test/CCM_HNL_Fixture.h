#pragma once
#ifndef SIREN_CCM_HNL_FIXTURE_H
#define SIREN_CCM_HNL_FIXTURE_H

// Shared assembly for the CCM dipole-HNL upper-injector chain: detector model,
// primary injection/physical processes, secondary N4 decay, and the resulting
// Injector + Weighter pair.  MakeCCMHNLFixture() is the single construction
// point for a generation-capable injector/weighter pair in CCM_HNL_TEST.cxx.
//
// Data-dependent: MakeCCMHNLFixture() loads material/detector files and
// cross-section tables from hardcoded /home/nwkamp paths and cannot run
// outside that environment.

#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <fstream>

#include "SIREN/utilities/Random.h"
#include "SIREN/dataclasses/Particle.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/injection/Injector.h"
#include "SIREN/injection/Process.h"
#include "SIREN/injection/Weighter.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"

#include "SIREN/distributions/primary/energy/Monoenergetic.h"
#include "SIREN/distributions/primary/direction/Cone.h"
#include "SIREN/distributions/primary/direction/IsotropicDirection.h"
#include "SIREN/distributions/primary/vertex/PointSourcePositionDistribution.h"
#include "SIREN/distributions/primary/helicity/PrimaryNeutrinoHelicityDistribution.h"
#include "SIREN/distributions/secondary/vertex/SecondaryPhysicalVertexDistribution.h"

#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/HNLDipoleFromTable.h"
#include "SIREN/interactions/HNLDipoleDecay.h"
#include "SIREN/interactions/Decay.h"

namespace siren { namespace injection { namespace test {

// generation settings
static const bool ccm_hnl_z_samp = false;
static const bool ccm_hnl_in_invGeV = true; // are the xsec tables in GeV^-1?
static const bool ccm_hnl_inelastic = true;

static const std::string ccm_hnl_tot_xsec_table_path = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/tot_xsec_Enu/";
static const std::string ccm_hnl_diff_xsec_table_path = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/";

static const std::string ccm_hnl_material_file = "/home/nwkamp/Research/CCM/DipoleAnalysis/sources/SIRENDevPrivate/resources/Detectors/materials/CCM.dat";
static const std::string ccm_hnl_detector_file = "/home/nwkamp/Research/CCM/DipoleAnalysis/sources/SIRENDevPrivate/resources/Detectors/densities/PREM_ccm.dat";

static const double ccm_hnl_mass = 0.01375; // in GeV; The HNL mass we are injecting
static const double ccm_hnl_dipole_coupling = 1.0e-6; // in GeV^-1; the effective dipole coupling strength
static const std::string ccm_hnl_mHNL = "0.01375";

// Events to inject
static const unsigned int ccm_hnl_events_to_inject = 100u;
static const siren::dataclasses::ParticleType ccm_hnl_primary_type = siren::dataclasses::ParticleType::NuMu;

inline std::vector<double> CCMHNLDipoleCouplingVec() {
    return std::vector<double>{
        (ccm_hnl_primary_type == siren::dataclasses::ParticleType::NuE) ? ccm_hnl_dipole_coupling : 0,
        (ccm_hnl_primary_type == siren::dataclasses::ParticleType::NuMu) ? ccm_hnl_dipole_coupling : 0,
        (ccm_hnl_primary_type == siren::dataclasses::ParticleType::NuTau) ? ccm_hnl_dipole_coupling : 0};
}

// Functions for loading xsec tables

inline std::string CCMHNLDiffXsPath(int Z, int A, std::string mHNL) {
    std::stringstream ss;
    ss << ccm_hnl_diff_xsec_table_path;
    if(ccm_hnl_z_samp) ss << "diff_xsec_z_Enu/";
    else ss << "diff_xsec_y_Enu/";
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

inline std::string CCMHNLTotXsPath(int Z, int A, std::string mHNL) {
    std::stringstream ss;
    ss << ccm_hnl_tot_xsec_table_path;
    ss << "xsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

inline std::vector<std::array<int, 2>> CCMHNLTargetZA() {
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

inline std::vector<siren::dataclasses::ParticleType> CCMHNLTargetPIDs() {
    using ParticleType = siren::dataclasses::ParticleType;
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

inline std::vector<std::string> CCMHNLDiffXsHfFiles(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : CCMHNLTargetZA()) {
        res.push_back(CCMHNLDiffXsPath(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

inline std::vector<std::string> CCMHNLTotXsHfFiles(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : CCMHNLTargetZA()) {
        res.push_back(CCMHNLTotXsPath(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

inline std::vector<std::string> CCMHNLDiffXsHcFiles(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : CCMHNLTargetZA()) {
        res.push_back(CCMHNLDiffXsPath(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

inline std::vector<std::string> CCMHNLTotXsHcFiles(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : CCMHNLTargetZA()) {
        res.push_back(CCMHNLTotXsPath(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

// True only when every hardcoded input MakeCCMHNLFixture() loads is present:
// the material/detector files and a representative differential/total
// cross-section table.  A test gates on this to GTEST_SKIP() rather than abort
// on environments without the private data.
inline bool CCMHNLDataPresent() {
    std::vector<std::string> const diff = CCMHNLDiffXsHfFiles(ccm_hnl_mHNL);
    std::vector<std::string> const tot = CCMHNLTotXsHfFiles(ccm_hnl_mHNL);
    return std::ifstream(ccm_hnl_material_file).good()
        && std::ifstream(ccm_hnl_detector_file).good()
        && !diff.empty() && std::ifstream(diff.front()).good()
        && !tot.empty() && std::ifstream(tot.front()).good();
}

// Bundles the pieces a BreakdownInvariant-style test needs: the upper
// injector (targets the upper tungsten converter), its matching physical
// weighter, and the shared detector model.
struct CCMHNLFixture {
    std::shared_ptr<siren::injection::Injector> injector;
    std::shared_ptr<siren::injection::Weighter> weighter;
    std::shared_ptr<siren::detector::DetectorModel> detector_model;
};

// Builds the upper-target injector/weighter pair using the same materials,
// cross sections, distributions, and stopping condition as the reference
// generation test.  Requires the hardcoded /home/nwkamp data files to exist.
inline CCMHNLFixture MakeCCMHNLFixture() {
    using ParticleType = siren::dataclasses::ParticleType;
    using namespace siren::distributions;
    using namespace siren::interactions;
    using namespace siren::injection;
    using namespace siren::detector;
    using namespace siren::utilities;

    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    detector_model->LoadMaterialModel(ccm_hnl_material_file);
    detector_model->LoadDetectorModel(ccm_hnl_detector_file);

    std::shared_ptr<SIREN_random> random = std::make_shared<SIREN_random>();

    std::shared_ptr<PrimaryInjectionProcess> primary_injection_process =
        std::make_shared<PrimaryInjectionProcess>();
    std::shared_ptr<PhysicalProcess> primary_physical_process =
        std::make_shared<PhysicalProcess>();
    std::vector<std::shared_ptr<PhysicalProcess>> secondary_physical_processes;
    primary_injection_process->SetPrimaryType(ccm_hnl_primary_type);
    primary_physical_process->SetPrimaryType(ccm_hnl_primary_type);

    // Cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<ParticleType> target_types = CCMHNLTargetPIDs();
    std::shared_ptr<HNLDipoleFromTable> hf_xs = std::make_shared<HNLDipoleFromTable>(
        ccm_hnl_mass, ccm_hnl_dipole_coupling, HNLDipoleFromTable::HelicityChannel::Flipping,
        ccm_hnl_z_samp, ccm_hnl_in_invGeV, ccm_hnl_inelastic);
    std::shared_ptr<HNLDipoleFromTable> hc_xs = std::make_shared<HNLDipoleFromTable>(
        ccm_hnl_mass, ccm_hnl_dipole_coupling, HNLDipoleFromTable::HelicityChannel::Conserving,
        ccm_hnl_z_samp, ccm_hnl_in_invGeV, ccm_hnl_inelastic);
    std::vector<std::string> hf_diff_fnames = CCMHNLDiffXsHfFiles(ccm_hnl_mHNL);
    std::vector<std::string> hc_diff_fnames = CCMHNLDiffXsHcFiles(ccm_hnl_mHNL);
    std::vector<std::string> hf_tot_fnames = CCMHNLTotXsHfFiles(ccm_hnl_mHNL);
    std::vector<std::string> hc_tot_fnames = CCMHNLTotXsHcFiles(ccm_hnl_mHNL);
    for(unsigned int i = 0; i < target_types.size(); ++i) {
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    std::shared_ptr<InteractionCollection> primary_interactions =
        std::make_shared<InteractionCollection>(ccm_hnl_primary_type, cross_sections);
    primary_injection_process->SetInteractions(primary_interactions);
    primary_physical_process->SetInteractions(primary_interactions);

    // Primary energy distribution: pion decay-at-rest
    double nu_energy = 0.02965;
    std::shared_ptr<PrimaryEnergyDistribution> edist = std::make_shared<Monoenergetic>(nu_energy);
    primary_injection_process->AddPrimaryInjectionDistribution(edist);
    primary_physical_process->AddPhysicalDistribution(edist);

    // Flux normalization: 2105.14020, 4.74e9 nu/m^2/s / (6.2e14 POT/s) * 4*pi*20m^2 -> nu/POT
    std::shared_ptr<WeightableDistribution> flux_units = std::make_shared<NormalizationConstant>(3.76e-2);
    primary_physical_process->AddPhysicalDistribution(flux_units);

    // Primary direction: cone toward the detector from the upper target
    double opening_angle = std::atan(2./23.);
    siren::math::Vector3D upper_target_origin(0, 0, 0.1375);
    siren::math::Vector3D detector_origin(23, 0, -0.65);
    siren::math::Vector3D upper_dir = detector_origin - upper_target_origin;
    upper_dir.normalize();
    std::shared_ptr<PrimaryDirectionDistribution> upper_inj_ddist = std::make_shared<Cone>(upper_dir, opening_angle);
    std::shared_ptr<PrimaryDirectionDistribution> phys_ddist = std::make_shared<IsotropicDirection>();
    primary_injection_process->AddPrimaryInjectionDistribution(upper_inj_ddist);
    primary_physical_process->AddPhysicalDistribution(phys_ddist);

    // Helicity distribution: this is a neutrino
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution =
        std::make_shared<PrimaryNeutrinoHelicityDistribution>();
    primary_injection_process->AddPrimaryInjectionDistribution(helicity_distribution);
    primary_physical_process->AddPhysicalDistribution(helicity_distribution);

    // Primary position distribution: treat the upper target as a point source
    double max_dist = 25; // m
    std::shared_ptr<VertexPositionDistribution> upper_pos_dist =
        std::make_shared<PointSourcePositionDistribution>(upper_target_origin, max_dist);
    primary_injection_process->AddPrimaryInjectionDistribution(upper_pos_dist);

    // Secondary process: dipole portal HNL (N4) decay
    std::shared_ptr<SecondaryInjectionProcess> secondary_decay_inj_process =
        std::make_shared<SecondaryInjectionProcess>();
    std::shared_ptr<PhysicalProcess> secondary_decay_phys_process =
        std::make_shared<PhysicalProcess>();
    secondary_decay_inj_process->SetPrimaryType(ParticleType::N4);
    secondary_decay_phys_process->SetPrimaryType(ParticleType::N4);

    std::shared_ptr<HNLDipoleDecay> sec_decay = std::make_shared<HNLDipoleDecay>(
        ccm_hnl_mass, CCMHNLDipoleCouplingVec(), HNLDipoleDecay::ChiralNature::Majorana);
    std::vector<std::shared_ptr<Decay>> sec_decays = {sec_decay};
    std::shared_ptr<InteractionCollection> secondary_interactions =
        std::make_shared<InteractionCollection>(ParticleType::N4, sec_decays);
    secondary_decay_inj_process->SetInteractions(secondary_interactions);
    secondary_decay_phys_process->SetInteractions(secondary_interactions);

    std::shared_ptr<SecondaryVertexPositionDistribution> secondary_pos_dist =
        std::make_shared<SecondaryPhysicalVertexDistribution>();
    secondary_decay_inj_process->AddSecondaryInjectionDistribution(secondary_pos_dist);

    std::vector<std::shared_ptr<SecondaryInjectionProcess>> secondary_injection_processes;
    secondary_injection_processes.push_back(secondary_decay_inj_process);
    secondary_physical_processes.push_back(secondary_decay_phys_process);

    std::shared_ptr<Injector> upper_injector = std::make_shared<Injector>(
        ccm_hnl_events_to_inject, detector_model, primary_injection_process,
        secondary_injection_processes, random);

    std::function<bool(siren::dataclasses::InteractionTree const &, std::shared_ptr<siren::dataclasses::InteractionTreeDatum>, size_t)> stopping_condition =
        [&](siren::dataclasses::InteractionTree const & tree, std::shared_ptr<siren::dataclasses::InteractionTreeDatum> datum, size_t i) {
            if(datum->depth(tree) >= 1) return true;
            return false;
        };
    upper_injector->SetStoppingCondition(stopping_condition);

    std::shared_ptr<Weighter> upper_weighter = std::make_shared<Weighter>(
        std::vector<std::shared_ptr<Injector>>{upper_injector}, detector_model,
        primary_physical_process, secondary_physical_processes);

    CCMHNLFixture fixture;
    fixture.injector = upper_injector;
    fixture.weighter = upper_weighter;
    fixture.detector_model = detector_model;
    return fixture;
}

} // namespace test
} // namespace injection
} // namespace siren

#endif // SIREN_CCM_HNL_FIXTURE_H
