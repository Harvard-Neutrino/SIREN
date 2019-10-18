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
        return( this->distribution(generator) );
    }

    double LI_random::Uniform( double from, double to ){
        if (to < from ){
            throw "'to' should be greater than 'from'";
        }

        double result = (from-to)*(this->distribution(generator)) + to;
        return( result );
    }

} //end namespace LeptonInjector 
