
#include <ios>
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
#include "LeptonInjector/Weighter.h"

#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/EulerQuaternionConversions.h"
#include "earthmodel-service/Placement.h"

//#define AUSTIN

using namespace LeptonInjector;

bool inFiducial(std::array<double,3> & int_vtx, earthmodel::ExtrPoly & fidVol) {
    earthmodel::Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    earthmodel::Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

TEST(ElasticScattering, Generation)
{
    using ParticleType = LeptonInjector::Particle::ParticleType;

#ifdef AUSTIN
    std::string material_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    std::string earth_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    std::string flux_file = "/home/austin/nu-dipole/fluxes/LE_FHC_numu.txt";
#else
    std::string material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    std::string earth_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    std::string flux_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/NUMI_Flux_Tables/ME_FHC_numu.txt";
#endif

    double max_distance = 240; // Maximum distance, set by distance from MiniBooNE to the decay pipe
    // This should encompass Minerva, should probably be smaller? Depends on how long Minerva is...
    double disk_radius = 1.4; // in meters
    double endcap_length = 5; // in meters


    // Events to inject
    unsigned int events_to_inject = 5e5;
    Particle::ParticleType primary_type = ParticleType::NuMu;

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::shared_ptr<ElasticScattering> es_xs = std::make_shared<ElasticScattering>();
    cross_sections.push_back(es_xs);

    // Load the earth model
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->LoadMaterialModel(material_file);
    earth_model->LoadEarthModel(earth_file);

    // Setup the primary type and mass
    //std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type, hnl_mass);
    std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type, 0);

    // Setup power law
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();

    // Setup tabulated flux
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(flux_file, true);
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf_gen = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(flux_file);

    // Change the flux units from cm^-2 to m^-2
    std::shared_ptr<LeptonInjector::WeightableDistribution> flux_units = std::make_shared<LeptonInjector::NormalizationConstant>(1e4);

    // Pick energy distribution
    std::shared_ptr<PrimaryEnergyDistribution> edist = tab_pdf_gen;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::FixedDirection>(earthmodel::Vector3D{0.0, 0.0, 1.0});

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to column depth
    std::shared_ptr<DepthFunction> depth_func = std::make_shared<LeptonInjector::LeptonDepthFunction>();

    // Helicity distribution
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<LeptonInjector::PrimaryNeutrinoHelicityDistribution>();

    // Put it all together!
    std::shared_ptr<InjectorBase> injector = std::make_shared<ColumnDepthLeptonInjector>(events_to_inject, primary_injector, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, depth_func, disk_radius, endcap_length, helicity_distribution);

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
    poly.push_back({0.0, 1.01758});
    poly.push_back({0.88125, 0.50879});
    poly.push_back({0.88125, -0.50879});
    poly.push_back({0.0, -1.01758});
    poly.push_back({-0.88125, -0.50879});
    poly.push_back({-0.88125, 0.50879});

    double offset[2];
    offset[0] = 0;
    offset[1] = 0;
    std::vector<earthmodel::ExtrPoly::ZSection> zsecs;
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(-2.0672,offset,1));
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(2.0672,offset,1));
    earthmodel::Placement placement(earthmodel::Vector3D(0,0,2.0672), earthmodel::QFromZXZr(0,0,0));
    earthmodel::ExtrPoly MINERvA_fiducial = earthmodel::ExtrPoly(placement, poly, zsecs);
    
    std::ofstream myFile("elasticscattering_test_events.csv");
    // myFile << std::fixed << std::setprecision(6);
    myFile << std::scientific << std::setprecision(16);
    myFile << "intX intY intZ ";
    myFile << "p4inu_0 p4inu_1 p4inu_2 p4inu_3 ";
    myFile << "helinu ";
    myFile << "p4iel_0 p4iel_1 p4iel_2 p4iel_3 ";
    myFile << "heliel ";
    myFile << "p4fnu_0 p4fnu_1 p4fnu_2 p4fnu_3 ";
    myFile << "helfnu ";
    myFile << "p4fel_0 p4fel_1 p4fel_2 p4fel_3 ";
    myFile << "helfel ";
    myFile << "simplified_weight interaction_prob y target fid\n";
    myFile << std::endl;
    int i = 0;
    while(*injector) {
        LeptonInjector::InteractionRecord event = injector->GenerateEvent();
        double simplified_weight, interaction_lengths, interaction_prob = 0;
        if(event.signature.target_type != LeptonInjector::Particle::ParticleType::unknown) {
            simplified_weight = weighter.SimplifiedEventWeight(event);
            interaction_prob = weighter.InteractionProbability(injector->InjectionBounds(event), event);
        }
        if(event.secondary_momenta.size() > 0) {
            myFile << event.interaction_vertex[0] << " ";
            myFile << event.interaction_vertex[1] << " ";
            myFile << event.interaction_vertex[2] << " ";

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

            myFile << simplified_weight << " ";
            myFile << interaction_prob << " ";
            myFile << event.interaction_parameters[0] << " "; // sampled y
            myFile << event.signature.target_type << " "; // target type
            myFile << int(inFiducial(event.interaction_vertex, MINERvA_fiducial)) << "\n"; // fid vol
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

