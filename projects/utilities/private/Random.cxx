#include "SIREN/utilities/Random.h"

#include <random>
#include <utility>

namespace siren {
namespace utilities {

    SIREN_random::SIREN_random(void){
        // default to boring seed
        seed   = 1;
        configuration       = std::default_random_engine(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0);
    }

    SIREN_random::SIREN_random( unsigned int _seed ){
        seed = _seed;
        configuration       = std::default_random_engine(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0);
    }

    // samples a number betwen the two specified values: (from, to)
    //      defaults to ( 0, 1)
    double SIREN_random::Uniform(double min, double max) {
        if(max < min)
            std::swap(min, max);
        double range = max - min;
        return range * (this->generator(configuration)) + min;
    }
    
    double SIREN_random::PowerLaw(double min, double max, double n) {
        if(max < min)
            std::swap(min, max);
        double range = max - min;
        double unif = range * (this->generator(configuration)) + min;
        
        double base = (pow(max,n+1) - pow(min,n+1)) * unif + pow(min,n+1) ; 
        double exp = 1/(n+1) ; 
        return pow(base, exp);
    }

    // reconfigures the generator with a new seed
    void SIREN_random::set_seed( unsigned int new_seed) {
        seed = new_seed;
        this->configuration = std::default_random_engine(seed);
    }

} // namespace utilities
} // namespace siren
