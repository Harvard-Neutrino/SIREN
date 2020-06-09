#include <Controller.h>
#include <Particle.h>
#include <LeptonInjector.h>
#include <Constants.h>
#include <string>

int main(void){

    int n_events = 55000;
    std::string diff_xs = "../test_xs.fits";
    std::string total_xs = "../test_xs_total.fits";
    bool is_ranged = true;

    LeptonInjector::Particle::ParticleType final_1 = LeptonInjector::Particle::MuMinus;
    LeptonInjector::Particle::ParticleType final_2 = LeptonInjector::Particle::Hadrons;

    LeptonInjector::MinimalInjectionConfiguration the_injector( n_events, final_1, final_2, diff_xs, total_xs, is_ranged);

    // units come from Constants.h 
    double minE = 1000.*LeptonInjector::Constants::GeV;
    double maxE = 100000.*LeptonInjector::Constants::GeV;
    double gamm = 2.;
    double minZenith = 80.*LeptonInjector::Constants::deg;
    double maxZenith = 180.*LeptonInjector::Constants::deg;
    double minAzimuth = 0.*LeptonInjector::Constants::deg;
    double maxAzimuth = 180.*LeptonInjector::Constants::deg;

    LeptonInjector::Controller cont(the_injector, minE, maxE, gamm, minAzimuth, maxAzimuth, minZenith, maxZenith, 1200.*LeptonInjector::Constants::m, 1200.*LeptonInjector::Constants::m, 1200.*LeptonInjector::Constants::m, 1200.*LeptonInjector::Constants::m);

    cont.NameOutfile("./data_output_c.h5");
    cont.NameLicFile("./config.lic");
    cont.Execute();
}
