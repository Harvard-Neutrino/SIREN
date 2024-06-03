#pragma once
#ifndef SIREN_Random_H
#define SIREN_Random_H

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// this implements a class to sample numbers just like in an i3 service

#include <random> // default_random_engine, uniform_real_distribution

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

namespace siren {
namespace utilities {

    class SIREN_random{
        public:
            SIREN_random();
            SIREN_random( unsigned int _seed );

            // this naming convention is used to
            double Uniform( double from=0.0, double to=1.0);
            double PowerLaw(double min, double max, double n); 

            // in case this is set up without a seed!
            void set_seed(unsigned int new_seed);

            template<class Archive>
            void save(Archive & archive, std::uint32_t const version) const {
            if(version == 0) {
                archive(::cereal::make_nvp("Seed", seed));
            } else {
                throw std::runtime_error("SIREN_random only supports version <= 0!");
            }
            };

            template<class Archive>
            void load(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(::cereal::make_nvp("Seed", seed));
                set_seed(seed);
            } else {
                throw std::runtime_error("SIREN_random only supports version <= 0!");
            }
            };

        private:
            unsigned int seed;
            std::default_random_engine configuration;
            std::uniform_real_distribution<double> generator;
    };

} // namespace utilities
} // namespace siren

#endif // SIREN_Random_H

