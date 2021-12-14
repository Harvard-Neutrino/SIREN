#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Particle.h>
#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Constants.h>
#include <earthmodel-service/EarthModel.h>
//#include <earthmodel-service/EarthModelService.h>
#include <earthmodel-service/Geometry.h>
#include <string>
#include <memory>
#include <chrono>
#include <ctime>
#include <argagg.hpp>
#include "date.h"

template <class Precision>
std::string getISOCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    return date::format("%FT%TZ", date::floor<Precision>(now));
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
            "dune_atmo_path", {"--dune-atmo-path"},
            "Path to the DUNEAtmo source directory. Used for the injection cross section.", 1,
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
            "ice_type", {"--ice", "--ice-type"},
            "The type of ice model.", 1,
        },
        {
            "ice_angle", {"--ice-angle", "--ice-cap-angle"},
            "Angle of the ice cap in degrees if it exists.", 1,
        },
        {
            "depth", {"--depth", "--detector-depth"},
            "Depth of the detector origin in meters.", 1,
        },
        {
            "minE", {"--min-energy", "--minE"},
            "The minimum injection energy in GeV", 1,
        },
        {
            "maxE", {"--max-energy", "--maxE"},
            "The maximum injection energy in GeV", 1,
        },
        {
            "gamma", {"--gamma"},
            "The injection spectral index", 1,
        },
        {
            "minZenith", {"--min-zenith", "--minZen"},
            "The minimum injection zenith angle in degrees", 1,
        },
        {
            "maxZenith", {"--max-zenith", "--maxZen"},
            "The maximum injection zenith angle in degrees", 1,
        },
        {
            "minAzimuth", {"--min-azimuth", "--minAzi"},
            "The minimum injection azimuth angle in degrees", 1,
        },
        {
            "maxAzimuth", {"--max-azimuth", "--maxAzi"},
            "The maximum injection azimuth angle in degrees", 1,
        },
        {
            "ranged_radius", {"--ranged-radius"},
            "The ranged mode radius of injection in meters", 1,
        },
        {
            "ranged_length", {"--ranged-length", "--end-cap"},
            "The ranged mode end cap length in meters", 1,
        },
        {
            "volume_radius", {"--volume-radius", "--cylinder-radius"},
            "The volume mode injection cylinder radius in meters", 1,
        },
        {
            "volume_height", {"--volume-height", "--cylinder-height"},
            "The volume mode injection cylinder height in meters", 1,
        },
        {
            "n_ranged", {"--n-ranged"},
            "Number of ranged events per neutrino type per interaction", 1,
        },
        {
            "n_volume", {"--n-volume"},
            "Number of volume events per neutrino type per interaction", 1,
        },
        {
            "cc", {"--cc"},
            "Inject CC events?", 0,
        },
        {
            "nc", {"--nc"},
            "Inject NC events?", 0,
        },
        {
            "gr", {"--gr"},
            "Inject GR events?", 0,
        },
        {
            "gr_e", {"--gr-e"},
            "Inject GR e events?", 0,
        },
        {
            "gr_mu", {"--gr-mu"},
            "Inject GR mu events?", 0,
        },
        {
            "gr_tau", {"--gr-tau"},
            "Inject GR tau events?", 0,
        },
        {
            "gr_had", {"--gr-had"},
            "Inject GR had events?", 0,
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
            "nue_cc_diff_xs", {"--nue-cc-dsdxdy"},
            "Path to the nue CC differential cross section file", 1,
        },
        {
            "nue_cc_total_xs", {"--nue-cc-sigma"},
            "Path to the nue CC total cross section file", 1,
        },
        {
            "nuebar_cc_diff_xs", {"--nuebar-cc-dsdxdy"},
            "Path to the nuebar CC differential cross section file", 1,
        },
        {
            "nuebar_cc_total_xs", {"--nuebar-cc-sigma"},
            "Path to the nuebar CC total cross section file", 1,
        },
        {
            "numu_cc_diff_xs", {"--numu-cc-dsdxdy"},
            "Path to the numu CC differential cross section file", 1,
        },
        {
            "numu_cc_total_xs", {"--numu-cc-sigma"},
            "Path to the numu CC total cross section file", 1,
        },
        {
            "numubar_cc_diff_xs", {"--numubar-cc-dsdxdy"},
            "Path to the numubar CC differential cross section file", 1,
        },
        {
            "numubar_cc_total_xs", {"--numubar-cc-sigma"},
            "Path to the numubar CC total cross section file", 1,
        },
        {
            "nutau_cc_diff_xs", {"--nutau-cc-dsdxdy"},
            "Path to the nutau CC differential cross section file", 1,
        },
        {
            "nutau_cc_total_xs", {"--nutau-cc-sigma"},
            "Path to the nutau CC total cross section file", 1,
        },
        {
            "nutaubar_cc_diff_xs", {"--nutaubar-cc-dsdxdy"},
            "Path to the nutaubar CC differential cross section file", 1,
        },
        {
            "nutaubar_cc_total_xs", {"--nutaubar-cc-sigma"},
            "Path to the nutaubar CC total cross section file", 1,
        },
        {
            "nue_nc_diff_xs", {"--nue-nc-dsdxdy"},
            "Path to the nue NC differential cross section file", 1,
        },
        {
            "nue_nc_total_xs", {"--nue-nc-sigma"},
            "Path to the nue NC total cross section file", 1,
        },
        {
            "nuebar_nc_diff_xs", {"--nuebar-nc-dsdxdy"},
            "Path to the nuebar NC differential cross section file", 1,
        },
        {
            "nuebar_nc_total_xs", {"--nuebar-nc-sigma"},
            "Path to the nuebar NC total cross section file", 1,
        },
        {
            "numu_nc_diff_xs", {"--numu-nc-dsdxdy"},
            "Path to the numu NC differential cross section file", 1,
        },
        {
            "numu_nc_total_xs", {"--numu-nc-sigma"},
            "Path to the numu NC total cross section file", 1,
        },
        {
            "numubar_nc_diff_xs", {"--numubar-nc-dsdxdy"},
            "Path to the numubar NC differential cross section file", 1,
        },
        {
            "numubar_nc_total_xs", {"--numubar-nc-sigma"},
            "Path to the numubar NC total cross section file", 1,
        },
        {
            "nutau_nc_diff_xs", {"--nutau-nc-dsdxdy"},
            "Path to the nutau NC differential cross section file", 1,
        },
        {
            "nutau_nc_total_xs", {"--nutau-nc-sigma"},
            "Path to the nutau NC total cross section file", 1,
        },
        {
            "nutaubar_nc_diff_xs", {"--nutaubar-nc-dsdxdy"},
            "Path to the nutaubar NC differential cross section file", 1,
        },
        {
            "nutaubar_nc_total_xs", {"--nutaubar-nc-sigma"},
            "Path to the nutaubar NC total cross section file", 1,
        },
        {
            "nuebar_gr_diff_xs", {"--nuebar-gr-dsdxdy"},
            "Path to the nuebar GR differential cross section file", 1,
        },
        {
            "nuebar_gr_total_xs", {"--nuebar-gr-sigma"},
            "Path to the nuebar GR total cross section file", 1,
        },
        {
            "nuebar_gr_e_diff_xs", {"--nuebar-gr-e-dsdxdy"},
            "Path to the nuebar GR e differential cross section file", 1,
        },
        {
            "nuebar_gr_e_total_xs", {"--nuebar-gr-e-sigma"},
            "Path to the nuebar GR e total cross section file", 1,
        },
        {
            "nuebar_gr_mu_diff_xs", {"--nuebar-gr-mu-dsdxdy"},
            "Path to the nuebar GR mu differential cross section file", 1,
        },
        {
            "nuebar_gr_mu_total_xs", {"--nuebar-gr-mu-sigma"},
            "Path to the nuebar GR mu total cross section file", 1,
        },
        {
            "nuebar_gr_tau_diff_xs", {"--nuebar-gr-tau-dsdxdy"},
            "Path to the nuebar GR tau differential cross section file", 1,
        },
        {
            "nuebar_gr_tau_total_xs", {"--nuebar-gr-tau-sigma"},
            "Path to the nuebar GR tau total cross section file", 1,
        },
        {
            "nuebar_gr_had_diff_xs", {"--nuebar-gr-had-dsdxdy"},
            "Path to the nuebar GR had differential cross section file", 1,
        },
        {
            "nuebar_gr_had_total_xs", {"--nuebar-gr-had-sigma"},
            "Path to the nuebar GR had total cross section file", 1,
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

    std::string earth_model = args["earth_model"].as<std::string>("PREM_mmc");
    std::string materials_model = args["materials_model"].as<std::string>("Standard");
    std::string ice_type = args["ice_type"].as<std::string>("SimpleIceCap");
    double ice_angle = args["ice_angle"].as<double>(20.);
    double depth = args["depth"].as<double>(1480.);

    // define some parameters shared by the injectors
    int n_ranged_events = int(1e5);
    int n_volume_events = int(1e5);

    n_ranged_events = int(args["n_ranged"].as<float>(float(n_ranged_events)));
    n_volume_events = int(args["n_volume"].as<float>(float(n_volume_events)));

    double minE = args["minE"].as<double>(1e2)*LeptonInjector::Constants::GeV;
    double maxE = args["maxE"].as<double>(1e6)*LeptonInjector::Constants::GeV;
    double gamma = args["gamma"].as<double>(2.);
    double minZenith = args["minZenith"].as<double>(0.)*LeptonInjector::Constants::degrees;
    double maxZenith = args["maxZenith"].as<double>(180.)*LeptonInjector::Constants::degrees;
    double minAzimuth = args["minAzimuth"].as<double>(0.)*LeptonInjector::Constants::degrees;
    double maxAzimuth = args["maxAzimuth"].as<double>(360.)*LeptonInjector::Constants::degrees;
    double ranged_radius = args["ranged_radius"].as<double>(32)*LeptonInjector::Constants::m;
    double ranged_length = args["ranged_length"].as<double>(32)*LeptonInjector::Constants::m;
    double volume_radius = args["volume_radius"].as<double>(32)*LeptonInjector::Constants::m;
    double volume_height = args["volume_height"].as<double>(20)*LeptonInjector::Constants::m;

    std::vector<std::string> interactions;
    std::vector<std::string> possible_interactions = {"cc", "nc"};
    for(unsigned int i=0; i<possible_interactions.size(); ++i) {
        if(args[possible_interactions[i]]) {
            interactions.push_back(possible_interactions[i]);
        }
    }

    if(args["gr"]) {
        interactions.push_back("gr_e");
        interactions.push_back("gr_mu");
        interactions.push_back("gr_tau");
        interactions.push_back("gr_had");
    } else {
        std::vector<std::string> possible_interaction_suffixes = {"e", "mu", "tau", "had"};
        for(unsigned int i=0; i<possible_interaction_suffixes.size(); ++i) {
            if(args["gr_" + possible_interaction_suffixes[i]]) {
                interactions.push_back("gr_" + possible_interaction_suffixes[i]);
            }
        }
    }


    std::vector<std::string> neutrinos;
    std::vector<std::string> possible_neutrinos = {"nue", "nuebar", "numu", "numubar", "nutau", "nutaubar"};
    for(unsigned int i=0; i<possible_neutrinos.size(); ++i) {
        if(args[possible_neutrinos[i]]) {
            neutrinos.push_back(possible_neutrinos[i]);
        }
    }

    std::map<std::string, std::pair<LeptonInjector::Particle::ParticleType, LeptonInjector::Particle::ParticleType>> secondaries = {
        {"nue_cc", {LeptonInjector::Particle::EMinus, LeptonInjector::Particle::Hadrons}},
        {"nuebar_cc", {LeptonInjector::Particle::EPlus, LeptonInjector::Particle::Hadrons}},
        {"numu_cc", {LeptonInjector::Particle::MuMinus, LeptonInjector::Particle::Hadrons}},
        {"numubar_cc", {LeptonInjector::Particle::MuPlus, LeptonInjector::Particle::Hadrons}},
        {"nutau_cc", {LeptonInjector::Particle::TauMinus, LeptonInjector::Particle::Hadrons}},
        {"nutaubar_cc", {LeptonInjector::Particle::TauPlus, LeptonInjector::Particle::Hadrons}},
        {"nue_nc", {LeptonInjector::Particle::NuE, LeptonInjector::Particle::Hadrons}},
        {"nuebar_nc", {LeptonInjector::Particle::NuEBar, LeptonInjector::Particle::Hadrons}},
        {"numu_nc", {LeptonInjector::Particle::NuMu, LeptonInjector::Particle::Hadrons}},
        {"numubar_nc", {LeptonInjector::Particle::NuMuBar, LeptonInjector::Particle::Hadrons}},
        {"nutau_nc", {LeptonInjector::Particle::NuTau, LeptonInjector::Particle::Hadrons}},
        {"nutaubar_nc", {LeptonInjector::Particle::NuTauBar, LeptonInjector::Particle::Hadrons}},
        {"nuebar_gr_e", {LeptonInjector::Particle::NuEBar, LeptonInjector::Particle::EMinus}},
        {"nuebar_gr_mu", {LeptonInjector::Particle::NuMuBar, LeptonInjector::Particle::MuMinus}},
        {"nuebar_gr_tau", {LeptonInjector::Particle::NuTauBar, LeptonInjector::Particle::TauMinus}},
        {"nuebar_gr_had", {LeptonInjector::Particle::Hadrons, LeptonInjector::Particle::Hadrons}},
    };

    std::cout << getISOCurrentTimestamp<std::chrono::seconds>() << std::endl;

    std::string xs_base;
    if(args["dune_atmo_path"]) {
        xs_base = args["dune_atmo_path"].as<std::string>();
    }
    else {
        if (const char* env_p = getenv("GOLEMSOURCEPATH")){
            xs_base = std::string( env_p ) + "/DUNEAtmo/";
        }
        else {
            std::cerr << "WARNING no DUNEAtmo path specified and GOLEMSOURCEPATH not set! Assuming the default cross section information is in ../cross_sections/" << std::endl;
            xs_base = "../";
        }
    }
    if(xs_base[-1] != '/')
        xs_base += "/";
    xs_base += "cross_sections/csms_differential_v1.0/";

    // specify the final state particles, and construct the first injector
    std::string nu_cc_diff_xs = xs_base + "dsdxdy_nu_CC_iso.fits";
    std::string nu_cc_total_xs = xs_base + "sigma_nu_CC_iso.fits";
    std::string nubar_cc_diff_xs = xs_base + "dsdxdy_nubar_CC_iso.fits";
    std::string nubar_cc_total_xs = xs_base + "sigma_nubar_CC_iso.fits";
    std::string nu_nc_diff_xs = xs_base + "dsdxdy_nu_NC_iso.fits";
    std::string nu_nc_total_xs = xs_base + "sigma_nu_NC_iso.fits";
    std::string nubar_nc_diff_xs = xs_base + "dsdxdy_nubar_NC_iso.fits";
    std::string nubar_nc_total_xs = xs_base + "sigma_nubar_NC_iso.fits";

    std::string nuebar_gr_diff_xs = xs_base + "dsdxdy_nuebar_GR.fits";
    std::string nuebar_gr_total_xs = xs_base + "sigma_nuebar_GR.fits";

    if(args["nuebar_gr_diff_xs"]) {
        nuebar_gr_diff_xs = args["nuebar_gr_diff_xs"].as<std::string>();
    }
    if(args["nuebar_gr_total_xs"]) {
        nuebar_gr_total_xs = args["nuebar_gr_total_xs"].as<std::string>();
    }

    std::map<std::string, std::pair<std::string, std::string>> default_cross_sections = {
        {"nue_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"nuebar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"numu_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"numubar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"nutau_cc", {nu_cc_diff_xs, nu_cc_total_xs}},
        {"nutaubar_cc", {nubar_cc_diff_xs, nubar_cc_total_xs}},
        {"nue_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"nuebar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"numu_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"numubar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"nutau_nc", {nu_nc_diff_xs, nu_nc_total_xs}},
        {"nutaubar_nc", {nubar_nc_diff_xs, nubar_nc_total_xs}},
        {"nuebar_gr_e", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_mu", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_tau", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
        {"nuebar_gr_had", {nuebar_gr_diff_xs, nuebar_gr_total_xs}},
    };

    std::vector<std::string> possible_interaction_suffixes = {"e", "mu", "tau", "had"};
    for(unsigned int i=0; i<possible_interaction_suffixes.size(); ++i) {
        bool replace = false;
        std::string diff = nuebar_gr_diff_xs;
        std::string total = nuebar_gr_total_xs;
        std::string id = "nuebar_gr_" + possible_interaction_suffixes[i];
        if(args[id + "_diff_xs"]) {
            diff = args[id + "_diff_xs"].as<std::string>();
            replace = true;
        }
        if(args[id + "_total_xs"]) {
            total = args[id + "_total_xs"].as<std::string>();
            replace = true;
        }
        if(replace) {
            default_cross_sections[id] = {diff, total};
        }
    }

    std::vector<LeptonInjector::InjectorBase> injectors;
    for(std::string & interaction : interactions) {
        bool is_glashow = interaction.rfind("gr", 0) == 0;
        if(args[interaction]) {
            for(std::string & neutrino : neutrinos) {
                if(args[neutrino]) {
                    if(is_glashow and neutrino != "nuebar") {
                        continue;
                    }
                    std::string id = neutrino + "_" + interaction;

                    std::string diff_xs_id = id + "_diff_xs";
                    std::string diff_xs =
                        args[diff_xs_id].as<std::string>(default_cross_sections[id].first);

                    std::string total_xs_id = id + "_total_xs";
                    std::string total_xs =
                        args[total_xs_id].as<std::string>(default_cross_sections[id].second);

                    LeptonInjector::Particle::ParticleType first_particle = secondaries[id].first;
                    LeptonInjector::Particle::ParticleType second_particle = secondaries[id].second;

                    if(n_ranged_events > 0) {
                        std::cout << "loading ranged injector\n";
                        LeptonInjector::RangedInjector injector(n_ranged_events,
                                first_particle,
                                //second_particle,
                                diff_xs,
                                total_xs,
                                true); // is_ranged
                        injectors.push_back(injector);
                    }
                    if(n_volume_events > 0) {
                        LeptonInjector::Injector injector(n_volume_events,
                                first_particle,
                                second_particle,
                                diff_xs,
                                total_xs,
                                false); // is_ranged
                        injectors.push_back(injector);
                    }

                }
                else {
                    continue;
                }
            }
        }
        else {
            continue;
        }
    }

    if(injectors.size() == 0) {
        std::cerr << "Must have at least one neutrino and one interaction type specified!" << std::endl;
        return EXIT_FAILURE;
    }

    // build the Controller object. This will facilitate the simulation itself
    // We need to pass the first injector while building this Controller
    // Dimensions of DUNE module is 58.2 x 3.5 x 12 x 4
    LeptonInjector::Controller cont(injectors[0],
            minE, maxE,
            gamma,
            minAzimuth, maxAzimuth,
            minZenith, maxZenith,
            ranged_radius, // injection radius
            ranged_length, // injection length
            volume_radius, // cylinder radius
            volume_height); // cylinder height


    for(unsigned int i=1; i<injectors.size(); ++i) {
        cont.AddInjector(injectors[i]);
    }


    cont.setSeed(seed);


    std::cout << "Mat model: " << materials_model << std::endl;
    try {
        earthmodel::EarthModel earthModel;
        earthModel.SetPath(path);
        earthModel.LoadMaterialModel(materials_model);
        //earthModel.LoadConcentricShellsFromLegacyFile(earth_model, depth*LeptonInjector::Constants::m, ice_angle*LeptonInjector::Constants::degrees);
        try{
					earthModel.LoadEarthModel(earth_model);
        }
        catch(const std::exception& err){
					std::cout << err.what();
        }

        std::vector<earthmodel::EarthSector> sectors = earthModel.GetSectors();
        for(unsigned int i=0; i<sectors.size(); ++i) {
            std::cerr << "Sector " << sectors[i].name << std::endl;
            std::cerr << "id: " << sectors[i].material_id << std::endl;
            std::cerr << "level: " << sectors[i].level << std::endl;
            std::shared_ptr<const earthmodel::Geometry> geo = sectors[i].geo;
            std::cerr << "Geo placement: " << geo->GetPlacement() << std::endl;
            std::cerr << "density: " << sectors[i].density->Evaluate(geo->GetPlacement().GetPosition()) << std::endl;
        }

				/*
        earthmodel::EarthModelService old_earthModel(
                "DUNE",
                path,
                std::vector<std::string>({earth_model}),
                std::vector<std::string>({materials_model}),
                ice_type,
                ice_angle*LeptonInjector::Constants::degrees,
                depth*LeptonInjector::Constants::m);

        */

        cont.SetEarthModel(std::shared_ptr<earthmodel::EarthModel>(&earthModel));
        //cont.SetEarthModel(std::shared_ptr<earthmodel::EarthModelService>(&old_earthModel));

        cont.NameOutfile(output + ".h5");
        cont.NameLicFile(output + ".lic");

        // Run the program.
        cont.Execute();
    } catch(char const * s) {
        std::cout << "Failure!" << std::endl;
        std::cout << s << std::endl;
    }
}
