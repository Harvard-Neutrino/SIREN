#include "phys-services/CrossSection.h"

#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Particle.h>
#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Constants.h>
#include "LeptonInjector/Weighter.h"

#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/EulerQuaternionConversions.h"
#include "earthmodel-service/Placement.h"

#include <string>
#include <iomanip>
#include <memory>
#include <chrono>
#include <ctime>
#include <argagg.hpp>
#include "date.h"

using namespace LeptonInjector;
bool z_samp = true;
bool in_invGeV = true;

template <class Precision>
std::string getISOCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    return date::format("%FT%TZ", date::floor<Precision>(now));
}

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

bool inFiducial(std::array<double,3> & int_vtx, earthmodel::ExtrPoly & fidVol) {
    earthmodel::Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    earthmodel::Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

bool inFiducial(std::array<double,3> & int_vtx, earthmodel::Sphere & fidVol) {
    earthmodel::Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    earthmodel::Vector3D dir(0,0,1);
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
            "flux_file", {"--flux-file"},
            "Path to NUMI flux file", 1,
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
            "FHC", {"--FHC"},
            "Inject according to NUMI FHC Flux?", 0,
        },
        {
            "RHC", {"--RHC"},
            "Inject according to NUMI RHC Flux?", 0,
        },
        {
            "LE", {"--LE"},
            "Inject according to NUMI LE Flux?", 0,
        },
        {
            "ME", {"--ME"},
            "Inject according to NUMI ME Flux?", 0,
        },
        {
            "nue", {"--nue"},
            "Inject nue events?", 0,
        },
        {
            "nuebar", {"--nuebar"},
            "Inject nuebar events?", 0,
        },
        {
            "numu", {"--numu"},
            "Inject numu events?", 0,
        },
        {
            "numubar", {"--numubar"},
            "Inject numubar events?", 0,
        },
        {
            "nutau", {"--nutau"},
            "Inject nutau events?", 0,
        },
        {
            "nutaubar", {"--nutaubar"},
            "Inject nutaubar events?", 0,
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

    std::string earth_file = args["earth_model"].as<std::string>("PREM_minerva");
    std::string materials_file = args["materials_model"].as<std::string>("Minerva");
    std::string flux_file;
    if(args["flux_file"]) {
        flux_file = args["flux_file"].as<std::string>();
    }
    else {
				std::cerr << "Please specify flux file!\n";
    }
    
    // default to MB best fit params
    double hnl_mass = args["hnl_mass"].as<double>(0.4);
    std::string mHNL = args["hnl_mass"].as<std::string>("0.4");
    double d_dipole = args["d_dipole"].as<double>(3e-7);

    // Decay parameters used to set the max range when injecting an HNL    
    double HNL_decay_width = std::pow(d_dipole,2)*std::pow(hnl_mass,3)/(4*Constants::pi); // in GeV; decay_width = d^2 m^3 / (4 * pi)
    double n_decay_lengths = 3.0;
    double max_distance = 240;

    // This should encompass Minerva, change to input argument if considering other detectors
    double disk_radius = 1.24; // in meters
    double endcap_length = 5; // in meters

    // Events to inject
    int events_to_inject = int(args["n_inj"].as<float>(float(1e5)));

    // Choose injection type (assumes only one is specified)
    Particle::ParticleType primary_type;
    if(args["nue"]) primary_type=Particle::ParticleType::NuE;
    else if(args["numu"]) primary_type=Particle::ParticleType::NuMu;
    else if(args["nutau"]) primary_type=Particle::ParticleType::NuTau;
    else if(args["nuebar"]) primary_type=Particle::ParticleType::NuEBar;
    else if(args["numubar"]) primary_type=Particle::ParticleType::NuMuBar;
    else if(args["nutaubar"]) primary_type=Particle::ParticleType::NuTauBar;
    else {std::cout << "Pick an input particle!\n"; exit(0);}

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
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass, d_dipole, DipoleFromTable::HelicityChannel::Flipping, z_samp, in_invGeV);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass, d_dipole, DipoleFromTable::HelicityChannel::Conserving, z_samp, in_invGeV);
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
    
    // Load the earth model
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->SetPath(path);
    earth_model->LoadMaterialModel(materials_file);
    earth_model->LoadEarthModel(earth_file);

    // Setup the primary type and mass
    std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type, 0);
    
    // Setup NUMI flux
    std::vector<double> params;
    if(args["LE"]) {
				if(args["FHC"]) {
						if(args["nue"]) params = {1.94e+00, 9.57e-01, 3.86e-01, 1.38e+01, 1.41e-01};
						if(args["numu"]) params = {2.06e+00, 6.52e-01, 3.36e-01, 7.50e+00, 1.19e-01};
						if(args["nuebar"]) params = {1.80e+00, 2.95e+00, 3.80e-01, 2.19e+01, 3.12e-01};
						if(args["numubar"]) params = {2.75e+00, 2.46e+00, 4.90e-01, 4.44e+00, 5.15e-02};
				}
				else if(args["RHC"]) {
						if(args["nue"]) params = {2.25e+00, 4.38e+00, 5.82e-01, 2.33e+02, 1.00e+00};
						if(args["numu"]) params = {3.75e+00, 3.04e+00, 5.53e-01, 1.50e+02, 3.11e-14};
						if(args["nuebar"]) params = {1.89e+00, 9.06e-01, 3.95e-01, 8.79e+00, 1.02e-01};
						if(args["numubar"]) params = {1.95e+00, 6.09e-01, 3.49e-01, 5.74e+00, 8.92e-02};
				}
    }
    else if(args["ME"] && args["numu"]) {
        params = {4.65e+00, 1.35e+00, 7.24e-02, 3.07e+00, 4.45e-03};
    }
    else {
        std::cout << "Please pick a valid flux! Exiting...\n";
        exit(0);
    }
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();
    std::shared_ptr<LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution> pdf = std::make_shared<LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution>(1.1*hnl_mass, 20, params[0], params[1], params[2], params[3], params[4]);
    
    // Setup tabulated flux
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(flux_file, true);
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf_gen = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(hnl_mass, 10, flux_file);

    // Change the flux units from cm^-2 to m^-2
    std::shared_ptr<LeptonInjector::WeightableDistribution> flux_units = std::make_shared<LeptonInjector::NormalizationConstant>(1e4);

    // Pick energy distribution
    std::shared_ptr<PrimaryEnergyDistribution> edist = tab_pdf_gen;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::FixedDirection>(earthmodel::Vector3D{0.0, 0.0, 1.0});

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to the decay distribution
    std::shared_ptr<RangeFunction> range_func = std::make_shared<LeptonInjector::DecayRangeFunction>(hnl_mass, HNL_decay_width, n_decay_lengths, max_distance);

    // Helicity distribution
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<LeptonInjector::PrimaryNeutrinoHelicityDistribution>();

    // Put it all together!
    std::shared_ptr<InjectorBase> injector = std::make_shared<RangedLeptonInjector>(events_to_inject, primary_injector, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length, helicity_distribution);

    std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions = {
        std::shared_ptr<WeightableDistribution>(tab_pdf),
        std::shared_ptr<WeightableDistribution>(flux_units),
        std::shared_ptr<WeightableDistribution>(ddist),
        std::shared_ptr<WeightableDistribution>(target_momentum_distribution),
        std::shared_ptr<WeightableDistribution>(helicity_distribution)
    };

    LeptonWeighter weighter(std::vector<std::shared_ptr<InjectorBase>>{injector}, earth_model, injector->GetCrossSections(), physical_distributions);
    
    // MINERvA Fiducial Volume
    std::vector<std::vector<double>> poly;
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
    std::vector<earthmodel::ExtrPoly::ZSection> zsecs;
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(0.125,offset,1));
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(4.1344,offset,1));
    earthmodel::Placement placement(earthmodel::Vector3D(0,0,0), earthmodel::QFromZXZr(0,0,0));
    earthmodel::ExtrPoly MINERvA_fiducial = earthmodel::ExtrPoly(placement, poly, zsecs);

    
    std::ofstream myFile(args["output"].as<std::string>()+".csv");
    myFile << std::scientific << std::setprecision(6);
    myFile << "intX intY intZ ";
    myFile << "p4gamma_0 p4gamma_1 p4gamma_2 p4gamma_3 ";
    myFile << "gamma_costh_hnlRest ";
    myFile << "decay_fid_weight prob_nopairprod simplified_weight fid ";
    if(!args["sparse"]) {
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
      myFile << "decay_length decay_ang_weight y target ";
    }
    myFile << std::endl << std::endl;
    int i=0;
    while(*injector) {
        LeptonInjector::InteractionRecord event = injector->GenerateEvent();
        LeptonInjector::DecayRecord decay;
        LeptonInjector::InteractionRecord pair_prod;
        double simplified_weight = 0;
        if(event.secondary_momenta.size() > 0) {
            
            injector->SampleSecondaryDecay(event, decay, HNL_decay_width, 0, 0, &MINERvA_fiducial, 0.1);
            injector->SamplePairProduction(decay, pair_prod);
            simplified_weight = weighter.SimplifiedEventWeight(event);
            
            myFile << event.interaction_vertex[0] << " ";
            myFile << event.interaction_vertex[1] << " ";
            myFile << event.interaction_vertex[2] << " ";

            myFile << decay.secondary_momenta[0][0] << " ";
            myFile << decay.secondary_momenta[0][1] << " ";
            myFile << decay.secondary_momenta[0][2] << " ";
            myFile << decay.secondary_momenta[0][3] << " ";
            myFile << decay.secondary_momenta[1][3]/decay.secondary_momenta[1][0] << " ";

            myFile << decay.decay_parameters[1] << " "; // decay fid weight
            myFile << pair_prod.interaction_parameters[0] << " "; // probability of no pair production
            myFile << simplified_weight << " ";
            myFile << int(inFiducial(pair_prod.interaction_vertex, MINERvA_fiducial)) << " "; // fid vol
            if(!args["sparse"]) {
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

              myFile << event.target_momentum[0] << " ";
              myFile << event.target_momentum[1] << " ";
              myFile << event.target_momentum[2] << " ";
              myFile << event.target_momentum[3] << " ";

              myFile << event.target_helicity << " ";

              myFile << event.secondary_momenta[0][0] << " ";
              myFile << event.secondary_momenta[0][1] << " ";
              myFile << event.secondary_momenta[0][2] << " ";
              myFile << event.secondary_momenta[0][3] << " ";

              myFile << event.secondary_helicity[0] << " ";

              myFile << event.secondary_momenta[1][0] << " ";
              myFile << event.secondary_momenta[1][1] << " ";
              myFile << event.secondary_momenta[1][2] << " ";
              myFile << event.secondary_momenta[1][3] << " ";

              myFile << event.secondary_helicity[1] << " ";
              myFile << decay.secondary_momenta[1][0] << " ";
              myFile << decay.secondary_momenta[1][1] << " ";
              myFile << decay.secondary_momenta[1][2] << " ";
              myFile << decay.secondary_momenta[1][3] << " ";

              myFile << decay.secondary_helicity[0] << " ";
              
              myFile << decay.decay_parameters[0] << " "; // decay length
              myFile << decay.decay_parameters[2] << " "; // decay ang weight
              myFile << event.interaction_parameters[1] << " "; // sampled y
              myFile << event.signature.target_type << " "; // target type
            }
            myFile << "\n\n";
        }
        if((++i)%int(events_to_inject/10.)==0) 
            std::cout << (int)(100*i/float(events_to_inject)) << "%" << std::endl;
    }
    myFile.close();
}
