#pragma once
#ifndef SIREN_Random_H
#define SIREN_Random_H

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// this implements a class to sample numbers just like in an i3 service

#include <random> // mt19937_64, uniform_real_distribution
#include <string>                                       // for string
#include <cstdint>                                      // for uint64_t

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
    SIREN_random(uint64_t _seed);

    // this naming convention is used to
    double Uniform( double from=0.0, double to=1.0);
    double PowerLaw(double min, double max, double n);

    // in case this is set up without a seed!
    void set_seed(uint64_t new_seed);
    uint64_t get_seed() const;

    static uint64_t generate_seed();

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version <= 2) {
            if(version >= 2) {
                std::string format = EngineStateFormat();
                archive(::cereal::make_nvp("EngineStateFormat", format));
            }
            archive(::cereal::make_nvp("Seed", seed));
            if(version >= 1) {
                std::string state_string = WriteEngineState();
                archive(::cereal::make_nvp("EngineState", state_string));
            }
        } else {
            throw std::runtime_error("SIREN_random only supports version <= 2!");
        }
    };

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version <= 2) {
            if(version >= 2) {
                std::string format;
                archive(::cereal::make_nvp("EngineStateFormat", format));
                if(format != EngineStateFormat())
                    throw std::runtime_error("SIREN_random: incompatible engine state format; "
                        "RNG archives require the same C++ standard library");
            }
            uint64_t loaded_seed;
            archive(::cereal::make_nvp("Seed", loaded_seed));
            std::mt19937_64 loaded_configuration(loaded_seed);
            if(version >= 1) {
                std::string state_string;
                archive(::cereal::make_nvp("EngineState", state_string));
                loaded_configuration = ReadEngineState(state_string);
            }
            seed = loaded_seed;
            configuration = loaded_configuration;
            generator = std::uniform_real_distribution<double>(0.0, 1.0);
        } else {
            throw std::runtime_error("SIREN_random only supports version <= 2!");
        }
    };

private:
    static std::string EngineStateFormat();
    static std::mt19937_64 ReadEngineState(std::string const & state);
    std::string WriteEngineState() const;

    uint64_t seed;
    // 64-bit Mersenne Twister (period 2^19937-1). The period must far exceed
    // (number of seeds) x (draws per event) x (events per seed) so that distinct
    // seeds never collide into overlapping draw sequences. A given seed
    // deterministically fixes the entire sequence; version 1 serializes the
    // engine state. Version 2 identifies its library-specific encoding;
    // version 0 (seed only) restarts from the seed.
    std::mt19937_64 configuration;
    std::uniform_real_distribution<double> generator;
};

} // namespace utilities
} // namespace siren

CEREAL_CLASS_VERSION(siren::utilities::SIREN_random, 2);

#endif // SIREN_Random_H
