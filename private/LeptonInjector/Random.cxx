#include <LeptonInjector/Random.h>

namespace LeptonInjector{

    LI_random::LI_random(void){
        // default to boring seed 
        unsigned int seed   = 1;
        configuration       = std::default_random_engine(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0); 
    }

    LI_random::LI_random( unsigned int seed ){
        configuration       = std::default_random_engine(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0); 
    }

    double LI_random::Uniform(void){
        return( this->generator(configuration) );
    }

    double LI_random::Uniform( double from, double to ){
        if (to < from ){
            throw "'to' should be greater than 'from'";
        }

        double result = (from-to)*(this->generator(configuration)) + to;
        return( result );
    }

    void LI_random::set_seed( unsigned int new_seed) {
        this->configuration = std::default_random_engine(new_seed);
    }

} //end namespace LeptonInjector 
