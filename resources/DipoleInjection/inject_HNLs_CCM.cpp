#include <ios>
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <argagg.hpp>
#include <string>
#include <iomanip>
#include <memory>
#include <chrono>
#include <ctime>


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

std::string diff_xs(int Z, int A, std::string mHNL, std::string diff_path) {
    std::stringstream ss;
    ss << diff_path;
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::string tot_xs(int Z, int A, std::string mHNL, std::string tot_path) {
    std::stringstream ss;
    ss << tot_path;
    ss << "xsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::vector<std::array<int, 2>> gen_ZA() {
    return std::vector<std::array<int, 2>>{
        {1, 1},
        {4, 9},
        {6, 12},
        {7, 14},
        {8, 16},
        {11, 23},
        {13, 27},
        {14, 28},
        {18, 40},
        {20, 40},
        {25, 55},
        {26, 56},
        {29, 63},
        {29, 65},
        {74, 183},
        {82, 208},
    };
}

std::vector<Particle::ParticleType> gen_TargetPIDs() {
    using ParticleType = Particle::ParticleType;
    return std::vector<ParticleType>{
        ParticleType::HNucleus,
        ParticleType::Be9Nucleus,
        ParticleType::C12Nucleus,
        ParticleType::N14Nucleus,
        ParticleType::O16Nucleus,
        ParticleType::Na23Nucleus,
        ParticleType::Al27Nucleus,
        ParticleType::Si28Nucleus,
        ParticleType::Ar40Nucleus,
        ParticleType::Ca40Nucleus,
        ParticleType::Mn55Nucleus,
        ParticleType::Fe56Nucleus,
        ParticleType::Cu63Nucleus,
        ParticleType::Cu65Nucleus,
        ParticleType::W183Nucleus,
        ParticleType::Pb208Nucleus
    };
}

std::vector<std::string> gen_diff_xs_hf(std::string mHNL, std::string diff_path) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL, diff_path) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hf(std::string mHNL, std::string tot_path) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL, tot_path) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_diff_xs_hc(std::string mHNL, std::string diff_path) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL, diff_path) + "_hc.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hc(std::string mHNL, std::string tot_path) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL, tot_path) + "_hc.dat");
    }
    return res;
}

bool inFiducial(std::array<double,3> & int_vtx, Geometry & fidVol) {
    Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

int main(int argc, char ** argv) {
    argagg::parser argparser = {{
        {
            "help", {"-h", "--help"},
            "Print help and exit", 0,
        },
        {
            "output", {"--output"},
            "Path and prefix for the output. Used for h5 and lic.", 1,
        },
        {
            "lepton_injector_path", {"--li", "--li-path", "--lepton-injector-path"},
            "Path to the lepton injector source directory. Used for earth model resources.", 1,
        },
        {
            "earth_model", {"--earth-model", "--earth-density", "--earth-density-model"},
            "Name of the earth density model.", 1
        },
        {
            "materials_model", {"--materials-model", "--materials-density", "--materials-density-model"},
            "Name of the materials density model.", 1,
        },
        {
            "tot_xsec_path", {"--tot-xsec-path"},
            "Path to tot xsec tables", 1,
        },
        {
            "dif_xsec_path", {"--diff-xsec-path"},
            "Path to diff xsec tables", 1,
        },
        {
            "z_samp", {"--z-samp","-z"},
            "flag to sample the differential cross section in z", 1,
        },
        {
            "sparse", {"--sparse","-s"},
            "flag for sparse output", 1,
        },
        {
            "n_inj", {"--n-inj","-n"},
            "number of events to inject", 1,
        },
        {
            "hnl_mass", {"--hnl-mass","-m"},
            "HNL mass in GeV", 1,
        },
        {
            "d_dipole", {"--d-dipole","-d"},
            "HNL effective dipole coupling in GeV^-1", 1,
        },
        {
            "seed", {"--seed"},
            "Seed for the injector", 1,
        },
    }};

    std::ostringstream usage;
    usage << argv[0] << std::endl;

    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        argagg::fmt_ostream fmt(std::cerr);
        fmt << usage.str() << argparser << std::endl
            << "Encountered exception while parsing arguments: " << e.what()
            << std::endl;
        return EXIT_FAILURE;
    }

    if(args["help"]) {
        std::cerr << argparser;
        return EXIT_SUCCESS;
    }

    double seed = args["seed"].as<int>(-1);

    if(seed < 0) {
        std::cerr << "--seed requires positive integer!" << std::endl << argparser;
        return EXIT_FAILURE;
    }

    if(not args["output"]) {
        std::cerr << "--output required!" << std::endl << argparser;
        return EXIT_FAILURE;
    }
    std::string output = args["output"].as<std::string>("./injected/output_DUNE");

    std::string path;
    if(args["lepton_injector_path"]) {
        path = args["lepton_injector_path"].as<std::string>();
    }
    else {
        if (const char* env_p = getenv("GOLEMSOURCEPATH")){
            path = std::string( env_p ) + "/LeptonInjectorDUNE/";
        }
        else {
            std::cerr << "WARNING no lepton injector path specified and GOLEMSOURCEPATH not set! Assuming earth model information is in ./resources/earthparams/" << std::endl;
            path = "./";
        }
    }
    if(path[-1] != '/')
        path += "/";
    path += "resources/earthparams";

    std::string earth_file = args["earth_model"].as<std::string>("PREM_ccm");
    std::string materials_file = args["materials_model"].as<std::string>("CCM");
    
    // Load the earth model
    std::shared_ptr<EarthModel> earth_model = std::make_shared<EarthModel>();
    earth_model->SetPath(path);
    earth_model->LoadMaterialModel(materials_file);
    earth_model->LoadEarthModel(earth_file);
    
    // random class instance
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();
    
    // HNL parameters
    double hnl_mass = args["hnl_mass"].as<double>(0.001);
    std::string mHNL = args["hnl_mass"].as<std::string>("0.001");
    double d_dipole = args["d_dipole"].as<double>(1e-6);

    // Events to inject
    int events_to_inject = int(args["n_inj"].as<float>(float(1e5)));

    // Choose injection type (assumes only one is specified)
    Particle::ParticleType primary_type=Particle::ParticleType::NuMu;

    // let's make the process instances
    // Injection processes
    std::shared_ptr<InjectionProcess> primary_injection_process_upper_target = std::make_shared<InjectionProcess>(); // will inject in upper tungsten target
    std::shared_ptr<InjectionProcess> primary_injection_process_lower_target = std::make_shared<InjectionProcess>(); // will inject in lower tungsten target
    std::vector<std::shared_ptr<InjectionProcess>> secondary_injection_processes; // common to both injectors
    // Physical processes
    std::shared_ptr<PhysicalProcess> primary_physical_process_upper_target = std::make_shared<PhysicalProcess>(); // will inject in upper tungsten target
    std::shared_ptr<PhysicalProcess> primary_physical_process_lower_target = std::make_shared<PhysicalProcess>(); // will inject in lower tungsten target
    std::vector<std::shared_ptr<PhysicalProcess>> secondary_physical_processes; // common to both injectors
    primary_injection_process_upper_target->primary_type = primary_type;
    primary_injection_process_lower_target->primary_type = primary_type;
    primary_physical_process_upper_target->primary_type = primary_type;
    primary_physical_process_lower_target->primary_type = primary_type;

    std::string tot_xs_base;
    std::string dif_xs_base;
    if(args["tot_xsec_path"]) {
        tot_xs_base = args["tot_xsec_path"].as<std::string>();
    }
		else {
				std::cerr << "WARNING no tot xsec path specified" << std::endl;
				exit(0);
		}
    if(args["dif_xsec_path"]) {
        dif_xs_base = args["dif_xsec_path"].as<std::string>();
    }
		else {
				std::cerr << "WARNING no tot xsec path specified" << std::endl;
				exit(0);
		} 
    if(tot_xs_base[-1] != '/') tot_xs_base += "/";
    if(dif_xs_base[-1] != '/') dif_xs_base += "/";

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<Particle::ParticleType> target_types = gen_TargetPIDs();
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass, d_dipole, DipoleFromTable::HelicityChannel::Flipping, z_samp, in_invGeV, inelastic);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass, d_dipole, DipoleFromTable::HelicityChannel::Conserving, z_samp, in_invGeV, inelastic);
    std::vector<std::string> hf_diff_fnames = gen_diff_xs_hf(mHNL,dif_xs_base);
    std::vector<std::string> hc_diff_fnames = gen_diff_xs_hc(mHNL,dif_xs_base);
    std::vector<std::string> hf_tot_fnames = gen_tot_xs_hf(mHNL,tot_xs_base);
    std::vector<std::string> hc_tot_fnames = gen_tot_xs_hc(mHNL,tot_xs_base);
    for(unsigned int i=0; i < target_types.size(); ++i) {
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    std::shared_ptr<CrossSectionCollection> primary_cross_sections = std::make_shared<CrossSectionCollection>(primary_type, cross_sections);
    primary_injection_process_upper_target->cross_sections = primary_cross_sections;
    primary_injection_process_lower_target->cross_sections = primary_cross_sections;
    primary_physical_process_upper_target->cross_sections = primary_cross_sections;
    primary_physical_process_lower_target->cross_sections = primary_cross_sections;

    // Primary energy distribution: pion decay-at-rest
    double nu_energy = 0.02965;
    std::shared_ptr<PrimaryEnergyDistribution> edist = std::make_shared<Monoenergetic>(nu_energy); // this creates a monoenergetic numu distribution
    primary_injection_process_upper_target->injection_distributions.push_back(edist);
    primary_injection_process_lower_target->injection_distributions.push_back(edist);
    primary_physical_process_upper_target->physical_distributions.push_back(edist);
    primary_physical_process_lower_target->physical_distributions.push_back(edist);
    
    // Flux normalization: using the number quoted in 2105.14020, 4.74e9 nu/m^2/s / (6.2e14 POT/s) * 4*pi*20m^2 to get nu/POT
    std::shared_ptr<WeightableDistribution> flux_units = std::make_shared<NormalizationConstant>(3.76e-2);
    primary_physical_process_upper_target->physical_distributions.push_back(flux_units);
    primary_physical_process_lower_target->physical_distributions.push_back(flux_units);

    // Primary direction: cone
    double opening_angle = std::atan(2./23.); // slightly larger than CCM xsec
    LI::math::Vector3D upper_target_origin(0, 0, 0.1375);
    LI::math::Vector3D lower_target_origin(0, 0, -0.241);
    LI::math::Vector3D detector_origin(23, 0, -0.65);
    LI::math::Vector3D upper_dir = detector_origin - upper_target_origin;
    upper_dir.normalize();
    LI::math::Vector3D lower_dir = detector_origin - lower_target_origin;
    lower_dir.normalize();
    std::shared_ptr<PrimaryDirectionDistribution> upper_inj_ddist = std::make_shared<Cone>(upper_dir,opening_angle);
    std::shared_ptr<PrimaryDirectionDistribution> lower_inj_ddist = std::make_shared<Cone>(lower_dir,opening_angle);
    std::shared_ptr<PrimaryDirectionDistribution> phys_ddist = std::make_shared<IsotropicDirection>(); // truly we are isotropic
    primary_injection_process_upper_target->injection_distributions.push_back(upper_inj_ddist);
    primary_injection_process_lower_target->injection_distributions.push_back(lower_inj_ddist);
    primary_physical_process_upper_target->physical_distributions.push_back(phys_ddist);
    primary_physical_process_lower_target->physical_distributions.push_back(phys_ddist);

    // Target momentum distribution: assume stationary for simplicity
    std::shared_ptr<TargetMomentumDistribution> target_momentum_distribution = std::make_shared<TargetAtRest>();
    primary_injection_process_upper_target->injection_distributions.push_back(target_momentum_distribution);
    primary_injection_process_lower_target->injection_distributions.push_back(target_momentum_distribution);
    primary_physical_process_upper_target->physical_distributions.push_back(target_momentum_distribution);
    primary_physical_process_lower_target->physical_distributions.push_back(target_momentum_distribution);

    // Helicity distribution: this is a neutrino
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<PrimaryNeutrinoHelicityDistribution>();
    primary_injection_process_upper_target->injection_distributions.push_back(helicity_distribution);
    primary_injection_process_lower_target->injection_distributions.push_back(helicity_distribution);
    primary_physical_process_upper_target->physical_distributions.push_back(helicity_distribution);
    primary_physical_process_lower_target->physical_distributions.push_back(helicity_distribution);

    // Primary position distribution: treat targets as point sources, generate from center
    double max_dist = 25; // m
    std::shared_ptr<VertexPositionDistribution> upper_pos_dist = std::make_shared<PointSourcePositionDistribution>(upper_target_origin, max_dist, primary_cross_sections->TargetTypes());
    std::shared_ptr<VertexPositionDistribution> lower_pos_dist = std::make_shared<PointSourcePositionDistribution>(lower_target_origin, max_dist, primary_cross_sections->TargetTypes());
    primary_injection_process_upper_target->injection_distributions.push_back(upper_pos_dist);
    primary_injection_process_lower_target->injection_distributions.push_back(lower_pos_dist);

    // Secondary process
    std::shared_ptr<InjectionProcess> secondary_decay_inj_process = std::make_shared<InjectionProcess>();
    std::shared_ptr<PhysicalProcess> secondary_decay_phys_process = std::make_shared<PhysicalProcess>();
    secondary_decay_inj_process->primary_type = ParticleType::NuF4;
    secondary_decay_phys_process->primary_type = ParticleType::NuF4;

    // Secondary cross sections: neutrissimo decay
    // Assume dirac HNL for now
    std::shared_ptr<NeutrissimoDecay> sec_decay = std::make_shared<NeutrissimoDecay>(hnl_mass, d_dipole, NeutrissimoDecay::ChiralNature::Majorana);
    std::vector<std::shared_ptr<Decay>> sec_decays = {sec_decay};
    std::shared_ptr<CrossSectionCollection> secondary_cross_sections = std::make_shared<CrossSectionCollection>(ParticleType::NuF4, sec_decays);
    secondary_decay_inj_process->cross_sections = secondary_cross_sections;
    secondary_decay_phys_process->cross_sections = secondary_cross_sections;

    // Secondary physical distribution
    std::shared_ptr<const LI::geometry::Geometry> fid_vol = NULL;
    for(auto sector : earth_model->GetSectors()) {
      if(sector.name=="ccm_inner_argon") fid_vol = sector.geo;
    }
    std::shared_ptr<VertexPositionDistribution> secondary_pos_dist = std::make_shared<SecondaryPositionDistribution>(fid_vol);
    secondary_decay_inj_process->injection_distributions.push_back(secondary_pos_dist);

    secondary_injection_processes.push_back(secondary_decay_inj_process);
    secondary_physical_processes.push_back(secondary_decay_phys_process);

    // Put it all together!
    std::shared_ptr<InjectorBase> upper_injector = std::make_shared<InjectorBase>(events_to_inject, earth_model, primary_injection_process_upper_target, secondary_injection_processes, random);
    std::shared_ptr<InjectorBase> lower_injector = std::make_shared<InjectorBase>(events_to_inject, earth_model, primary_injection_process_lower_target, secondary_injection_processes, random);

    // Set stopping condition
    std::function<bool(std::shared_ptr<LI::dataclasses::InteractionTreeDatum>)> stopping_condition =
      [&] (std::shared_ptr<LI::dataclasses::InteractionTreeDatum> datum) {
        if(datum->depth() >=1) return true;
        return false;
    };
    upper_injector->SetStoppingCondition(stopping_condition);
    lower_injector->SetStoppingCondition(stopping_condition);

    std::shared_ptr<LeptonTreeWeighter> upper_weighter = std::make_shared<LeptonTreeWeighter>(std::vector<std::shared_ptr<InjectorBase>>{upper_injector}, earth_model, primary_physical_process_upper_target, secondary_physical_processes);
    std::shared_ptr<LeptonTreeWeighter> lower_weighter = std::make_shared<LeptonTreeWeighter>(std::vector<std::shared_ptr<InjectorBase>>{lower_injector}, earth_model, primary_physical_process_lower_target, secondary_physical_processes);
    
    // Upper injector
    std::ofstream myFile(args["output"].as<std::string>()+"_upper.csv");
    myFile << std::scientific << std::setprecision(4);
    myFile << "event interaction vertex_x vertex_y vertex_z ";
    myFile << "primary_ID ";
    myFile << "primary_p0 primary_p1 primary_p2 primary_p3 ";
    myFile << "secondary1_ID ";
    myFile << "secondary1_p0 secondary1_p1 secondary1_p2 secondary1_p3 ";
    myFile << "secondary2_ID ";
    myFile << "secondary2_p0 secondary2_p2 secondary2_p2 secondary2_p3 ";
    myFile << "weight ";
    myFile << std::endl;
    int i=0;
    while(*upper_injector) {
        InteractionTree tree = upper_injector->GenerateEvent();
        double weight = upper_weighter->EventWeight(tree);
        int j = 0;
        for(auto datum : tree.tree) {
          myFile << i << " " << j << " ";
          myFile << datum->record.interaction_vertex[0] << " ";
          myFile << datum->record.interaction_vertex[1] << " ";
          myFile << datum->record.interaction_vertex[2] << " ";
          myFile <<  datum->record.signature.primary_type << " ";
          myFile << datum->record.primary_momentum[0] << " ";
          myFile << datum->record.primary_momentum[1] << " ";
          myFile << datum->record.primary_momentum[2] << " ";
          myFile << datum->record.primary_momentum[3] << " ";
          myFile <<  datum->record.signature.secondary_types[0] << " ";
          myFile << datum->record.secondary_momenta[0][0] << " ";
          myFile << datum->record.secondary_momenta[0][1] << " ";
          myFile << datum->record.secondary_momenta[0][2] << " ";
          myFile << datum->record.secondary_momenta[0][3] << " ";
          myFile <<  datum->record.signature.secondary_types[1] << " ";
          myFile << datum->record.secondary_momenta[1][0] << " ";
          myFile << datum->record.secondary_momenta[1][1] << " ";
          myFile << datum->record.secondary_momenta[1][2] << " ";
          myFile << datum->record.secondary_momenta[1][3] << " ";
          myFile << weight << std::endl;
          ++j;
        }
        myFile << std::endl;
        if((++i)%int(events_to_inject/10.)==0) 
            std::cout << (int)(100*i/float(events_to_inject)) << "%" << std::endl;
    }
    myFile.close();
    
    // Lower injector
    myFile = std::ofstream(args["output"].as<std::string>()+"_lower.csv");
    myFile << std::scientific << std::setprecision(4);
    myFile << "event interaction vertex_x vertex_y vertex_z ";
    myFile << "primary_ID ";
    myFile << "primary_p0 primary_p1 primary_p2 primary_p3 ";
    myFile << "secondary1_ID ";
    myFile << "secondary1_p0 secondary1_p1 secondary1_p2 secondary1_p3 ";
    myFile << "secondary2_ID ";
    myFile << "secondary2_p0 secondary2_p2 secondary2_p2 secondary2_p3 ";
    myFile << "weight ";
    myFile << std::endl;
    i=0;
    while(*lower_injector) {
        InteractionTree tree = lower_injector->GenerateEvent();
        double weight = lower_weighter->EventWeight(tree);
        int j = 0;
        for(auto datum : tree.tree) {
          myFile << i << " " << j << " ";
          myFile << datum->record.interaction_vertex[0] << " ";
          myFile << datum->record.interaction_vertex[1] << " ";
          myFile << datum->record.interaction_vertex[2] << " ";
          myFile <<  datum->record.signature.primary_type << " ";
          myFile << datum->record.primary_momentum[0] << " ";
          myFile << datum->record.primary_momentum[1] << " ";
          myFile << datum->record.primary_momentum[2] << " ";
          myFile << datum->record.primary_momentum[3] << " ";
          myFile <<  datum->record.signature.secondary_types[0] << " ";
          myFile << datum->record.secondary_momenta[0][0] << " ";
          myFile << datum->record.secondary_momenta[0][1] << " ";
          myFile << datum->record.secondary_momenta[0][2] << " ";
          myFile << datum->record.secondary_momenta[0][3] << " ";
          myFile <<  datum->record.signature.secondary_types[1] << " ";
          myFile << datum->record.secondary_momenta[1][0] << " ";
          myFile << datum->record.secondary_momenta[1][1] << " ";
          myFile << datum->record.secondary_momenta[1][2] << " ";
          myFile << datum->record.secondary_momenta[1][3] << " ";
          myFile << weight << std::endl;
          ++j;
        }
        myFile << std::endl;
        if((++i)%int(events_to_inject/10.)==0) 
            std::cout << (int)(100*i/float(events_to_inject)) << "%" << std::endl;
    }
    myFile.close();
}
