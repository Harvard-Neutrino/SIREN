#pragma once
#ifndef SIREN_Random_H
#define SIREN_Random_H

// Ben Smithers
// benjamin.smithers@mavs.uta.edu

// this implements a class to sample numbers just like in an i3 service

#include <random> // mt19937_64, uniform_real_distribution
#include <sstream>                                      // for engine-state (de)serialization
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
        if(version <= 1) {
            archive(::cereal::make_nvp("Seed", seed));
            // Version 1 also serializes the full engine state so a loaded
            // generator RESUMES the stream where it left off rather than
            // restarting from the seed. mt19937_64's stream insertion emits
            // its complete internal state (the state buffer plus position).
            if(version >= 1) {
                std::ostringstream state;
                state << configuration;
                std::string state_string = state.str();
                archive(::cereal::make_nvp("EngineState", state_string));
            }
        } else {
            throw std::runtime_error("SIREN_random only supports version <= 1!");
        }
    };

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version <= 1) {
            archive(::cereal::make_nvp("Seed", seed));
            // Version 0 archives carry only the seed, so re-seed a fresh engine
            // (the historical restart-from-seed behavior). Version 1 restores
            // the exact engine state, leaving the seed as a label; the
            // distribution is stateless, so resetting it to [0, 1) suffices.
            if(version >= 1) {
                std::string state_string;
                archive(::cereal::make_nvp("EngineState", state_string));
                std::istringstream state(state_string);
                state >> configuration;
                generator = std::uniform_real_distribution<double>(0.0, 1.0);
            } else {
                set_seed(seed);
            }
        } else {
            throw std::runtime_error("SIREN_random only supports version <= 1!");
        }
    };

private:
    uint64_t seed;
    // 64-bit Mersenne Twister (period 2^19937-1). The period must far exceed
    // (number of seeds) x (draws per event) x (events per seed) so that distinct
    // seeds never collide into overlapping draw sequences. A given seed
    // deterministically fixes the entire sequence; version 1 serializes the
    // engine state so a reloaded generator resumes mid-sequence, while version 0
    // (seed only) restarts from the seed.
    std::mt19937_64 configuration;
    std::uniform_real_distribution<double> generator;
};

} // namespace utilities
} // namespace siren

CEREAL_CLASS_VERSION(siren::utilities::SIREN_random, 1);

#endif // SIREN_Random_H
