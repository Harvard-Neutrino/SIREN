
#include <ios>
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "LeptonInjector/interactions/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Constants.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/injection/Injector.h"
#include "LeptonInjector/injection/RangedLeptonInjector.h"
#include "LeptonInjector/injection/Weighter.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/ExtrPoly.h"
#include "LeptonInjector/geometry/Sphere.h"
#include "LeptonInjector/math/EulerQuaternionConversions.h"
#include "LeptonInjector/geometry/Placement.h"

#include "LeptonInjector/distributions/primary/type/PrimaryInjector.h"
#include "LeptonInjector/distributions/primary/energy/ModifiedMoyalPlusExponentialEnergyDistribution.h"
#include "LeptonInjector/distributions/primary/energy/TabulatedFluxDistribution.h"
#include "LeptonInjector/distributions/primary/direction/FixedDirection.h"
#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/DipoleFromTable.h"

#define AUSTIN

using namespace LI::math;
using namespace LI::geometry;
using namespace LI::detector;
using namespace LI::injection;
using namespace LI::dataclasses;
using namespace LI::interactions;
using namespace LI::utilities;
using namespace LI::distributions;

bool z_samp = false;
bool in_invGeV = true;
bool inelastic = false;
bool miniboone = false;


std::string diff_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
#ifdef AUSTIN
    ss << "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
#else
    if(z_samp) ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/diff_xsec_z_Enu/";
    else ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
#endif
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::string tot_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
#ifdef AUSTIN
    ss << "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/";
#else
    ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/tot_xsec_Enu/";
#endif
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

std::vector<ParticleType> gen_TargetPIDs() {
    using ParticleType = ParticleType;
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

bool inFiducial(std::array<double,3> & int_vtx, ExtrPoly & fidVol) {
    Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

bool inFiducial(std::array<double,3> & int_vtx, Sphere & fidVol) {
    Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

double ComputeInteractionLengths(std::shared_ptr<DetectorModel const> detector_model, std::shared_ptr<InteractionCollection const> interactions, std::tuple<Vector3D, Vector3D> const & bounds, InteractionRecord const & record) {
    Vector3D interaction_vertex = record.interaction_vertex;
    Vector3D direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    direction.normalize();

    Geometry::IntersectionList intersections = detector_model->GetIntersections(detector_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), direction);
	std::map<ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = interactions->GetCrossSectionsByTarget();
    std::vector<double> total_cross_sections;
    std::vector<ParticleType> targets;
	InteractionRecord fake_record = record;
	for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
		fake_record.target_mass = detector_model->GetTargetMass(target_xs.first);
		std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
		double total_xs = 0.0;
		for(auto const & xs : xs_list) {
			std::vector<InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
			for(auto const & signature : signatures) {
				fake_record.signature = signature;
				// Add total cross section
				total_xs += xs->TotalCrossSection(fake_record);
			}
		}
		total_cross_sections.push_back(total_xs);
	}
    std::vector<double> particle_depths = detector_model->GetParticleColumnDepth(intersections, std::get<0>(bounds), std::get<1>(bounds), targets);

    double interaction_depth = 0.0;
    for(unsigned int i=0; i<targets.size(); ++i) {
        interaction_depth += particle_depths[i] * total_cross_sections[i];
    }
    return interaction_depth;
}

std::vector<double> p_LE_FHC_nue = {1.94e+00, 9.57e-01, 3.86e-01, 1.38e+01, 1.41e-01};
std::vector<double> p_LE_FHC_numu = {2.06e+00, 6.52e-01, 3.36e-01, 7.50e+00, 1.19e-01};
std::vector<double> p_LE_FHC_nuebar = {1.80e+00, 2.95e+00, 3.80e-01, 2.19e+01, 3.12e-01};
std::vector<double> p_LE_FHC_numubar = {2.75e+00, 2.46e+00, 4.90e-01, 4.44e+00, 5.15e-02};

std::vector<double> p_LE_RHC_nue = {2.25e+00, 4.38e+00, 5.82e-01, 2.33e+02, 1.00e+00};
std::vector<double> p_LE_RHC_numu = {3.75e+00, 3.04e+00, 5.53e-01, 1.50e+02, 3.11e-14};
std::vector<double> p_LE_RHC_nuebar = {1.89e+00, 9.06e-01, 3.95e-01, 8.79e+00, 1.02e-01};
std::vector<double> p_LE_RHC_numubar = {1.95e+00, 6.09e-01, 3.49e-01, 5.74e+00, 8.92e-02};

std::vector<double> p_ME_FHC_numu = {4.65e+00, 1.35e+00, 7.24e-02, 3.07e+00, 4.45e-03};


TEST(Injector, Generation)
{
    using ParticleType = ParticleType;

#ifdef AUSTIN
    std::string material_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/Detectors/materials/Minerva.dat";
    std::string detector_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/Detectors/densities/PREM_minerva.dat";
    std::string flux_file = "/home/austin/nu-dipole/fluxes/LE_FHC_numu.txt";
    z_samp = false;
    in_invGeV = false;
#else
    std::string material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/Detectors/materials/Minerva.dat";
    std::string detector_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/Detectors/densities/PREM_minerva.dat";
    std::string flux_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/NUMI_Flux_Tables/ME_FHC_numu.txt";
    if(miniboone) {
			material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/Detectors/materials/MiniBooNE.dat";
			detector_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/Detectors/densities/PREM_miniboone.dat";
			flux_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/BNB_Flux_Tables/BNB_numu_flux.txt";
			inelastic = true;
    }
#endif

    double hnl_mass = 0.4; // in GeV; The HNL mass we are injecting
    double dipole_coupling = 1.0e-6; // in GeV^-1; the effective dipole coupling strength
    std::string mHNL = "0.4";

    // Decay parameters used to set the max range when injecting an HNL
    double HNL_decay_width = std::pow(dipole_coupling,2)*std::pow(hnl_mass,3)/(4*Constants::pi); // in GeV; decay_width = d^2 m^3 / (4 * pi)
    double n_decay_lengths = 3.0; // Number of decay lengths to consider
    double max_distance = 240; // Maximum distance, set by distance from MiniBooNE to the decay pipe
    // This should encompass Minerva, should probably be smaller? Depends on how long Minerva is...
    double disk_radius = 1.4; // in meters
    double endcap_length = 5; // in meters
    if(miniboone) {
			max_distance = 541;
			disk_radius = 6.2;
			endcap_length = 6.2;
    }


    // Events to inject
    unsigned int events_to_inject = 5e5;
    ParticleType primary_type = ParticleType::NuMu;

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<ParticleType> primary_types = {ParticleType::NuE, ParticleType::NuMu, ParticleType::NuTau};
    std::vector<ParticleType> target_types = gen_TargetPIDs();
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass, dipole_coupling, DipoleFromTable::HelicityChannel::Flipping, z_samp, in_invGeV, inelastic);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass, dipole_coupling, DipoleFromTable::HelicityChannel::Conserving, z_samp, in_invGeV, inelastic);
    std::vector<std::string> hf_diff_fnames = gen_diff_xs_hf(mHNL);
    std::vector<std::string> hc_diff_fnames = gen_diff_xs_hc(mHNL);
    std::vector<std::string> hf_tot_fnames = gen_tot_xs_hf(mHNL);
    std::vector<std::string> hc_tot_fnames = gen_tot_xs_hc(mHNL);
    for(unsigned int i=0; i < target_types.size(); ++i) {
        //std::cerr << hf_diff_fnames[i] << std::endl;
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        ////std::cerr << hf_tot_fnames[i] << std::endl;
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        //std::cerr << hc_diff_fnames[i] << std::endl;
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        //std::cerr << hc_tot_fnames[i] << std::endl;
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    // Load the detector model
    std::shared_ptr<DetectorModel> detector_model = std::make_shared<DetectorModel>();
    detector_model->LoadMaterialModel(material_file);
    detector_model->LoadDetectorModel(detector_file);

    // Setup the primary type and mass
    //std::shared_ptr<PrimaryInjector> primary_injector = std::make_shared<PrimaryInjector>(primary_type, hnl_mass);
    std::shared_ptr<PrimaryInjector> primary_injector = std::make_shared<PrimaryInjector>(primary_type, 0);

    // Setup power law
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();

    std::vector<double> moyal_exp_params = p_ME_FHC_numu;

    // Setup NUMI flux
    std::shared_ptr<ModifiedMoyalPlusExponentialEnergyDistribution> pdf = std::make_shared<ModifiedMoyalPlusExponentialEnergyDistribution>(hnl_mass, 20, moyal_exp_params[0], moyal_exp_params[1], moyal_exp_params[2], moyal_exp_params[3], moyal_exp_params[4]);

    // Setup tabulated flux
    std::shared_ptr<TabulatedFluxDistribution> tab_pdf = std::make_shared<TabulatedFluxDistribution>(flux_file, true);
    std::shared_ptr<TabulatedFluxDistribution> tab_pdf_gen = std::make_shared<TabulatedFluxDistribution>(hnl_mass, 10, flux_file);

    // Change the flux units from cm^-2 to m^-2
    std::shared_ptr<WeightableDistribution> flux_units = std::make_shared<NormalizationConstant>(1e4);

    // Pick energy distribution
    std::shared_ptr<PrimaryEnergyDistribution> edist = tab_pdf_gen;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<FixedDirection>(Vector3D{0.0, 0.0, 1.0});

    // Let us inject according to the decay distribution
    std::shared_ptr<RangeFunction> range_func = std::make_shared<DecayRangeFunction>(hnl_mass, HNL_decay_width, n_decay_lengths, max_distance);

    // Helicity distribution
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<PrimaryNeutrinoHelicityDistribution>();

    // Put it all together!
    std::shared_ptr<Injector> injector = std::make_shared<RangedLeptonInjector>(events_to_inject, primary_injector, cross_sections, detector_model, random, edist, ddist, range_func, disk_radius, endcap_length, helicity_distribution);

    std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions = {
        std::shared_ptr<WeightableDistribution>(tab_pdf),
        std::shared_ptr<WeightableDistribution>(flux_units),
        std::shared_ptr<WeightableDistribution>(ddist),
        std::shared_ptr<WeightableDistribution>(helicity_distribution)
    };

    LeptonWeighter weighter(std::vector<std::shared_ptr<Injector>>{injector}, detector_model, injector->GetInteractions(), physical_distributions);

    std::vector<std::vector<double>> poly;
    // MINERvA Fiducial Volume
    // 88.125 cm apothem
    //poly.push_back({0.0, 1.01758});
    //poly.push_back({0.88125, 0.50879});
    //poly.push_back({0.88125, -0.50879});
    //poly.push_back({0.0, -1.01758});
    //poly.push_back({-0.88125, -0.50879});
    //poly.push_back({-0.88125, 0.50879});

    // 81.125 cm apothem
    poly.push_back({0.0, 0.93675});
    poly.push_back({0.81125, 0.46838});
    poly.push_back({0.81125, -0.46838});
    poly.push_back({0.0, -0.93675});
    poly.push_back({-0.81125, -0.46838});
    poly.push_back({-0.81125, 0.46838});

    double offset[2];
    offset[0] = 0;
    offset[1] = 0;
    std::vector<ExtrPoly::ZSection> zsecs;
    zsecs.push_back(ExtrPoly::ZSection(1.45,offset,1));
    zsecs.push_back(ExtrPoly::ZSection(4.0,offset,1));
    Placement placement(Vector3D(0,0,0), QFromZXZr(0,0,0));
    ExtrPoly MINERvA_fiducial = ExtrPoly(placement, poly, zsecs);

    // MiniBooNE Fiducial Volume
    Placement placementMB(Vector3D(0,0,0), QFromZXZr(0,0,0));
    Sphere MiniBooNE_fiducial = Sphere(placement, 5.0, 0.0);

    std::ofstream myFile("injector_test_events.csv");
    // myFile << std::fixed << std::setprecision(6);
    myFile << std::scientific << std::setprecision(16);
    myFile << "intX intY intZ ";
    myFile << "decX decY decZ ";
    myFile << "ppX ppY ppZ ";
    myFile << "p4nu_0 p4nu_1 p4nu_2 p4nu_3 ";
    myFile << "helnu ";
    myFile << "p4itgt_0 p4itgt_1 p4itgt_2 p4itgt_3 ";
    myFile << "helitgt ";
    myFile << "p4hnl_0 p4hnl_1 p4hnl_2 p4hnl_3 ";
    myFile << "helhnl ";
    myFile << "p4ftgt_0 p4ftgt_1 p4ftgt_2 p4ftgt_3 ";
    myFile << "helftgt ";
    myFile << "p4gamma_0 p4gamma_1 p4gamma_2 p4gamma_3 ";
    myFile << "p4gamma_hnlRest_0 p4gamma_hnlRest_1 p4gamma_hnlRest_2 p4gamma_hnlRest_3 ";
    myFile << "helgamma ";
    myFile << "decay_length decay_fid_weight decay_ang_weight prob_nopairprod basic_weight simplified_weight interaction_lengths interaction_prob y target fid\n";
    myFile << std::endl;
    int i = 0;
    while(*injector) {
        InteractionRecord event = injector->GenerateEvent();
        DecayRecord decay;
        InteractionRecord pair_prod;
        double basic_weight, simplified_weight, interaction_lengths, interaction_prob = 0;
        if(event.signature.target_type != ParticleType::unknown) {
            if(miniboone) injector->SampleSecondaryDecay(event, decay, HNL_decay_width, 1, 0, &MiniBooNE_fiducial, 0.1);
            else injector->SampleSecondaryDecay(event, decay, HNL_decay_width, 1, 0, &MINERvA_fiducial, 1.0);
            //injector->SamplePairProduction(decay, pair_prod);
            //basic_weight = weighter.EventWeight(event);
            simplified_weight = weighter.SimplifiedEventWeight(event);
            interaction_lengths = ComputeInteractionLengths(detector_model, injector->GetInteractions(), injector->InjectionBounds(event), event);
            interaction_prob = weighter.InteractionProbability(injector->InjectionBounds(event), event);
        }
        if(event.secondary_momenta.size() > 0) {
            myFile << event.interaction_vertex[0] << " ";
            myFile << event.interaction_vertex[1] << " ";
            myFile << event.interaction_vertex[2] << " ";

            myFile << decay.decay_vertex[0] << " ";
            myFile << decay.decay_vertex[1] << " ";
            myFile << decay.decay_vertex[2] << " ";

            myFile << pair_prod.interaction_vertex[0] << " ";
            myFile << pair_prod.interaction_vertex[1] << " ";
            myFile << pair_prod.interaction_vertex[2] << " ";

            myFile << event.primary_momentum[0] << " ";
            myFile << event.primary_momentum[1] << " ";
            myFile << event.primary_momentum[2] << " ";
            myFile << event.primary_momentum[3] << " ";

            myFile << event.primary_helicity << " ";

            myFile << event.target_helicity << " ";

            myFile << event.secondary_momenta[0][0] << " ";
            myFile << event.secondary_momenta[0][1] << " ";
            myFile << event.secondary_momenta[0][2] << " ";
            myFile << event.secondary_momenta[0][3] << " ";

            myFile << event.secondary_helicities[0] << " ";

            myFile << event.secondary_momenta[1][0] << " ";
            myFile << event.secondary_momenta[1][1] << " ";
            myFile << event.secondary_momenta[1][2] << " ";
            myFile << event.secondary_momenta[1][3] << " ";

            myFile << event.secondary_helicities[1] << " ";

            myFile << decay.secondary_momenta[0][0] << " ";
            myFile << decay.secondary_momenta[0][1] << " ";
            myFile << decay.secondary_momenta[0][2] << " ";
            myFile << decay.secondary_momenta[0][3] << " ";

            myFile << decay.secondary_momenta[1][0] << " ";
            myFile << decay.secondary_momenta[1][1] << " ";
            myFile << decay.secondary_momenta[1][2] << " ";
            myFile << decay.secondary_momenta[1][3] << " ";

            myFile << decay.secondary_helicities[0] << " ";

            myFile << decay.decay_parameters[3] << " "; // decay enter
            myFile << decay.decay_parameters[4] << " "; // decay exit
            myFile << decay.decay_parameters[0] << " "; // decay length
            myFile << decay.decay_parameters[1] << " "; // decay fid weight
            myFile << decay.decay_parameters[2] << " "; // decay anglular weight
            myFile << pair_prod.interaction_parameters[0] << " "; // probability of no pair production
            myFile << basic_weight << " ";
            myFile << simplified_weight << " ";
            myFile << interaction_lengths << " ";
            myFile << interaction_prob << " ";
            myFile << event.interaction_parameters[1] << " "; // sampled y
            myFile << event.signature.target_type << " "; // target type
            if(miniboone) myFile << int(inFiducial(pair_prod.interaction_vertex, MiniBooNE_fiducial)) << "\n"; // fid vol
            else myFile << int(inFiducial(pair_prod.interaction_vertex, MINERvA_fiducial)) << "\n"; // fid vol
            myFile << "\n";
        }
        if((++i) % (events_to_inject/10)==0)
            std::cerr << (int)(i*100 / (events_to_inject)) << "%" << std::endl;
    }
    myFile.close();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

