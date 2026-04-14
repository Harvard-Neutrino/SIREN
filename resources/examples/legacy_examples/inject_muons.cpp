#include <SIREN/Controller.h>
#include <SIREN/Particle.h>
#include <SIREN/SIREN.h>
#include <SIREN/Constants.h>
#include <string>

/*
This is an example constructed similarly to the python example.

It prepares a run with two injectors, one with CC muon events resulting from a neutrino interaction and one from a nubar interaction.

The output then creates two files:
    - "data_output_c.h5" which has all the data about the interactions
    - "confic_c.lic" which contains a 'snapshot' of the generation settings used to make the aforementioned data. This is to be loaded into LeptonWeighter to weight the events
*/

int main(void){

    // define some parameters shared by both of the injectors 
    int n_events = 25000;
    std::string diff_xs = "../test_xs.fits";
    std::string total_xs = "../test_xs_total.fits";
    bool is_ranged = true;

    // specify the final state particles, and construct the first injector 
    SIREN::Particle::ParticleType final_1 = SIREN::Particle::MuMinus;
    SIREN::Particle::ParticleType final_2 = SIREN::Particle::Hadrons;
    SIREN::Injector the_injector( n_events, final_1, final_2, diff_xs, total_xs, is_ranged);

    // specify more final state particles and contruct the second injector 
    SIREN::Particle::ParticleType final_3 = SIREN::Particle::MuPlus;
    SIREN::Injector next_injector( n_events, final_3, final_2, diff_xs, total_xs, is_ranged);

    // some global values shared by all the injectors 
    // units come from Constants.h 
    double minE = 1000.*SIREN::Constants::GeV;
    double maxE = 100000.*SIREN::Constants::GeV;
    double gamm = 2.;
    double minZenith = 80.*SIREN::Constants::deg;
    double maxZenith = 180.*SIREN::Constants::deg;
    double minAzimuth = 0.*SIREN::Constants::deg;
    double maxAzimuth = 180.*SIREN::Constants::deg;

    // build the Controller object. This will facilitate the simulation itself
    // We need to pass the first injector while building this Controller 
    SIREN::Controller cont(the_injector, minE, maxE, gamm, minAzimuth, maxAzimuth, minZenith, maxZenith, 1200.*SIREN::Constants::m, 1200.*SIREN::Constants::m, 1200.*SIREN::Constants::m, 1200.*SIREN::Constants::m);

    // Add the second injector, specify the output 
    cont.AddInjector(next_injector);
    cont.NameOutfile("./data_output_c.h5");
    cont.NameLicFile("./config_c.lic");

    // Run the program. 
    cont.Execute();
}
