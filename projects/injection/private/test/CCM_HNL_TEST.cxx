
#include <ios>
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "LeptonInjector/crosssections/CrossSection.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Constants.h"
#include "LeptonInjector/dataclasses/Particle.h"
#include "LeptonInjector/injection/InjectorBase.h"
#include "LeptonInjector/injection/RangedLeptonInjector.h"
#include "LeptonInjector/injection/TreeWeighter.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/ExtrPoly.h"
#include "LeptonInjector/geometry/Sphere.h"
#include "LeptonInjector/math/EulerQuaternionConversions.h"
#include "LeptonInjector/geometry/Placement.h"

#include "LeptonInjector/distributions/primary/energy/Monoenergetic.h"
#include "LeptonInjector/distributions/primary/direction/Cone.h"
#include "LeptonInjector/distributions/primary/direction/IsotropicDirection.h"
#include "LeptonInjector/distributions/primary/vertex/PointSourcePositionDistribution.h"
#include "LeptonInjector/distributions/primary/vertex/SecondaryPositionDistribution.h"

#include "LeptonInjector/crosssections/CrossSectionCollection.h"
#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/DipoleFromTable.h"
#include "LeptonInjector/crosssections/Decay.h"
#include "LeptonInjector/crosssections/NeutrissimoDecay.h"

using namespace LI::math;
using namespace LI::geometry;
using namespace LI::detector;
using namespace LI::injection;
using namespace LI::dataclasses;
using namespace LI::crosssections;
using namespace LI::utilities;
using namespace LI::distributions;

using ParticleType = Particle::ParticleType;

// Variables which the user will have to specify

// generation settings
bool z_samp = false;
bool in_invGeV = true; // are the xsec tables in GeV^-1?
bool inelastic = true;

std::string tot_xsec_table_path = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/tot_xsec_Enu/";
std::string diff_xsec_table_path = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables/";

std::string material_file = "/home/nwkamp/Research/CCM/DipoleAnalysis/sources/LeptonInjectorDevPrivate/resources/earthparams/materials/CCM.dat";
std::string earth_file = "/home/nwkamp/Research/CCM/DipoleAnalysis/sources/LeptonInjectorDevPrivate/resources/earthparams/densities/PREM_ccm.dat";
    
double hnl_mass = 0.01375; // in GeV; The HNL mass we are injecting
double dipole_coupling = 1.0e-6; // in GeV^-1; the effective dipole coupling strength
std::string mHNL = "0.01375";
    
// Events to inject
unsigned int events_to_inject = 5e5;
Particle::ParticleType primary_type = ParticleType::NuMu;
std::vector<double> dipole_coupling_vec = {(primary_type==ParticleType::NuE) ? dipole_coupling : 0,
                                           (primary_type==ParticleType::NuMu) ? dipole_coupling : 0,
                                           (primary_type==ParticleType::NuTau) ? dipole_coupling : 0};

// Functions for loading xsec tables
std::string diff_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
    ss << diff_xsec_table_path;
    if(z_samp) ss << "diff_xsec_z_Enu/";
    else ss << "diff_xsec_y_Enu/";
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::string tot_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
    ss << tot_xsec_table_path;
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

std::vector<Particle::ParticleType> gen_TargetPIDs() {
    using ParticleType = Particle::ParticleType;
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

double ComputeInteractionLengths(std::shared_ptr<EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::pair<Vector3D, Vector3D> const & bounds, InteractionRecord const & record) {
    Vector3D interaction_vertex = record.interaction_vertex;
    Vector3D direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    direction.normalize();

    Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), direction);
	std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<double> total_cross_sections;
    std::vector<Particle::ParticleType> targets;
	InteractionRecord fake_record = record;
	for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
		fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
		fake_record.target_momentum = {fake_record.target_mass,0,0,0};
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
    std::vector<double> particle_depths = earth_model->GetParticleColumnDepth(intersections, bounds.first, bounds.second, targets);

    double interaction_depth = 0.0;
    for(unsigned int i=0; i<targets.size(); ++i) {
        interaction_depth += particle_depths[i] * total_cross_sections[i];
    }
    return interaction_depth;
}


TEST(Injector, Generation)
{
    using ParticleType = Particle::ParticleType;


    // Load the earth model
    std::shared_ptr<EarthModel> earth_model = std::make_shared<EarthModel>();
    std::cout << "LoadMaterialModel...\n";
    earth_model->LoadMaterialModel(material_file);
    std::cout << "LoadEarthModel...\n";
    earth_model->LoadEarthModel(earth_file);
    std::cout << "Loaded EarthModel!\n";

    // random class instance
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();

    // let's make the process instances
    // Injection processes
    std::shared_ptr<InjectionProcess> primary_injection_process_upper_injector = std::make_shared<InjectionProcess>(); // will inject in upper tungsten target
    std::shared_ptr<InjectionProcess> primary_injection_process_lower_injector = std::make_shared<InjectionProcess>(); // will inject in lower tungsten target
    std::vector<std::shared_ptr<InjectionProcess>> secondary_injection_processes; // common to both injectors
    // Physical processes
    std::shared_ptr<PhysicalProcess> primary_physical_process_upper_injector = std::make_shared<PhysicalProcess>(); // will inject in upper tungsten target
    std::shared_ptr<PhysicalProcess> primary_physical_process_lower_injector = std::make_shared<PhysicalProcess>(); // will inject in lower tungsten target
    std::vector<std::shared_ptr<PhysicalProcess>> secondary_physical_processes; // common to both injectors
    primary_injection_process_upper_injector->primary_type = primary_type;
    primary_injection_process_lower_injector->primary_type = primary_type;
    primary_physical_process_upper_injector->primary_type = primary_type;
    primary_physical_process_lower_injector->primary_type = primary_type;
    
    std::cout << "LoadingCrossSections...\n";
    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<Particle::ParticleType> target_types = gen_TargetPIDs();
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
    
    std::cout << "GotCrossSections!\n";

    std::shared_ptr<CrossSectionCollection> primary_cross_sections = std::make_shared<CrossSectionCollection>(primary_type, cross_sections);
    primary_injection_process_upper_injector->cross_sections = primary_cross_sections;
    primary_injection_process_lower_injector->cross_sections = primary_cross_sections;
    primary_physical_process_upper_injector->cross_sections = primary_cross_sections;
    primary_physical_process_lower_injector->cross_sections = primary_cross_sections;

    std::cout << "PrimaryEnergyDistribution...\n";
    // Primary energy distribution: pion decay-at-rest
    double nu_energy = 0.02965;
    std::shared_ptr<PrimaryEnergyDistribution> edist = std::make_shared<Monoenergetic>(nu_energy); // this creates a monoenergetic numu distribution
    primary_injection_process_upper_injector->injection_distributions.push_back(edist);
    primary_injection_process_lower_injector->injection_distributions.push_back(edist);
    primary_physical_process_upper_injector->physical_distributions.push_back(edist);
    primary_physical_process_lower_injector->physical_distributions.push_back(edist);

    // Flux normalization: using the number quoted in 2105.14020, 4.74e5 nu/cm^2/s / (6.2e14 POT/s)
    std::shared_ptr<WeightableDistribution> flux_units = std::make_shared<NormalizationConstant>(7.48e-10);
    primary_physical_process_upper_injector->physical_distributions.push_back(flux_units);
    primary_physical_process_lower_injector->physical_distributions.push_back(flux_units);

    std::cout << "PrimaryDirectionDistribution...\n";
    // Primary direction: cone
    double opening_angle = std::cos(std::atan(5./23.)); // slightly larger than CCM xsec
    std::shared_ptr<PrimaryDirectionDistribution> inj_ddist = std::make_shared<Cone>(Vector3D{1.0, 0.0, 0.0},opening_angle);
    std::shared_ptr<PrimaryDirectionDistribution> phys_ddist = std::make_shared<IsotropicDirection>(); // truly we are isotropic
    primary_injection_process_upper_injector->injection_distributions.push_back(inj_ddist);
    primary_injection_process_lower_injector->injection_distributions.push_back(inj_ddist);
    primary_physical_process_upper_injector->physical_distributions.push_back(phys_ddist);
    primary_physical_process_lower_injector->physical_distributions.push_back(phys_ddist);

    // Target momentum distribution: assume stationary for simplicity
    std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution = std::make_shared<TargetAtRest>();
    primary_injection_process_upper_injector->injection_distributions.push_back(target_momentum_distribution);
    primary_injection_process_lower_injector->injection_distributions.push_back(target_momentum_distribution);
    primary_physical_process_upper_injector->physical_distributions.push_back(target_momentum_distribution);
    primary_physical_process_lower_injector->physical_distributions.push_back(target_momentum_distribution);

    // Helicity distribution: this is a neutrino
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<PrimaryNeutrinoHelicityDistribution>();
    primary_injection_process_upper_injector->injection_distributions.push_back(helicity_distribution);
    primary_injection_process_lower_injector->injection_distributions.push_back(helicity_distribution);
    primary_physical_process_upper_injector->physical_distributions.push_back(helicity_distribution);
    primary_physical_process_lower_injector->physical_distributions.push_back(helicity_distribution);

    // Primary position distribution: treat targets as point sources, generate from center
    LI::math::Vector3D upper_target_origin(0, 0, 0.1375);
    LI::math::Vector3D lower_target_origin(0, 0, -0.241);
    double max_dist = 25; // m
    std::shared_ptr<VertexPositionDistribution> upper_pos_dist = std::make_shared<PointSourcePositionDistribution>(upper_target_origin, max_dist, primary_cross_sections->TargetTypes()); 
    std::shared_ptr<VertexPositionDistribution> lower_pos_dist = std::make_shared<PointSourcePositionDistribution>(lower_target_origin, max_dist, primary_cross_sections->TargetTypes());
    primary_injection_process_upper_injector->injection_distributions.push_back(upper_pos_dist);
    primary_injection_process_lower_injector->injection_distributions.push_back(lower_pos_dist);
    primary_physical_process_upper_injector->physical_distributions.push_back(upper_pos_dist);
    primary_physical_process_lower_injector->physical_distributions.push_back(lower_pos_dist);

    // Secondary process
    std::shared_ptr<InjectionProcess> secondary_decay_inj_process = std::make_shared<InjectionProcess>();
    std::shared_ptr<PhysicalProcess> secondary_decay_phys_process = std::make_shared<PhysicalProcess>();
    secondary_decay_inj_process->primary_type = ParticleType::NuF4;
    secondary_decay_phys_process->primary_type = ParticleType::NuF4;
    
    // Secondary cross sections: neutrissimo decay
    // Assume dirac HNL for now
    std::shared_ptr<NeutrissimoDecay> sec_decay = std::make_shared<NeutrissimoDecay>(hnl_mass, dipole_coupling_vec, NeutrissimoDecay::ChiralNature::Dirac);  
    std::vector<std::shared_ptr<Decay>> sec_decays = {sec_decay};
    std::shared_ptr<CrossSectionCollection> secondary_cross_sections = std::make_shared<CrossSectionCollection>(ParticleType::NuF4, sec_decays);
    secondary_decay_inj_process->cross_sections = secondary_cross_sections;
    secondary_decay_phys_process->cross_sections = secondary_cross_sections;

    // Secondary physical distribution
    std::shared_ptr<VertexPositionDistribution> secondary_pos_dist = std::make_shared<SecondaryPositionDistribution>();
    secondary_decay_inj_process->injection_distributions.push_back(secondary_pos_dist);

    secondary_injection_processes.push_back(secondary_decay_inj_process);
    secondary_physical_processes.push_back(secondary_decay_phys_process);

    // Put it all together!
    std::shared_ptr<InjectorBase> upper_injector = std::make_shared<InjectorBase>(events_to_inject, earth_model, primary_injection_process_upper_injector, secondary_injection_processes, random);
    std::shared_ptr<InjectorBase> lower_injector = std::make_shared<InjectorBase>(events_to_inject, earth_model, primary_injection_process_lower_injector, secondary_injection_processes, random);

    // Set stopping condition
    std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> stopping_condition = 
      [&] (std::shared_ptr<LI::dataclasses::InteractionTreeDatum> datum) {
        if(datum->depth() >=1) return true;
        return false;
    };
    upper_injector->SetStoppingCondition(stopping_condition);
    lower_injector->SetStoppingCondition(stopping_condition);

    std::shared_ptr<LeptonTreeWeighter> upper_weighter = std::make_shared<LeptonTreeWeighter>(std::vector<std::shared_ptr<InjectorBase>>{upper_injector}, earth_model, primary_physical_process_upper_injector, secondary_physical_processes);
    std::shared_ptr<LeptonTreeWeighter> lower_weighter = std::make_shared<LeptonTreeWeighter>(std::vector<std::shared_ptr<InjectorBase>>{lower_injector}, earth_model, primary_physical_process_lower_injector, secondary_physical_processes);


    int i = 0;
    while(*upper_injector) {
        std::cout << "\nEvent " << i << std::endl;
        InteractionTree tree = upper_injector->GenerateEvent();
        for(auto datum : tree.tree) {
          std::cout << "\nDatum for " << datum->record.signature.primary_type << std::endl;
          std::cout << "Vertex: ";
          std::cout << datum->record.interaction_vertex[0] << " ";
          std::cout << datum->record.interaction_vertex[1] << " ";
          std::cout << datum->record.interaction_vertex[2] << std::endl;
          std::cout << "Momentum: ";
          std::cout << datum->record.primary_momentum[0] << " ";
          std::cout << datum->record.primary_momentum[1] << " ";
          std::cout << datum->record.primary_momentum[2] << " ";
          std::cout << datum->record.primary_momentum[3] << std::endl;
        }
        double weight = upper_weighter->EventWeight(tree);
        ++i;
    }
    while(*lower_injector) {
        InteractionTree tree = lower_injector->GenerateEvent();
        double weight = lower_weighter->EventWeight(tree);
        ++i;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

