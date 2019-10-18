#ifndef LI_RANDOM
#define LI_RANDOM

// this implements a class to sample numbers just like in an i3 service
// Essentially, this will act as a wrapper object to dictate random number sampling 

#include <random>

namespace LeptonInjector {

    class LI_random{
        public:
            LI_random();
            LI_random( unsigned int seed );

            // this naming convention is used to
            double Uniform();
            double Uniform( double from, double to);

            // in case this is set up without a seed! 
            void set_seed(unsigned int new_seed);

        private: 
            std::default_random_engine configuration;
            std::uniform_real_distribution<double> generator;
    };

} //end namespace LeptonInjector

#endif
