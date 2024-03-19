#pragma once
#ifndef SIREN_Random_H
#define SIREN_Random_H

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// this implements a class to sample numbers just like in an i3 service

#include <random> // default_random_engine, uniform_real_distribution

namespace siren {
namespace utilities {

    class SIREN_random{
        public:
            SIREN_random();
            SIREN_random( unsigned int seed );

            // this naming convention is used to
            double Uniform( double from=0.0, double to=1.0);
            double PowerLaw(double min, double max, double n); 

            // in case this is set up without a seed!
            void set_seed(unsigned int new_seed);

        private:
            std::default_random_engine configuration;
            std::uniform_real_distribution<double> generator;
    };

} // namespace utilities
} // namespace siren

#endif // SIREN_Random_H

