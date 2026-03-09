#include "SIREN/utilities/Random.h"

#include <mutex>
#include <atomic>
#include <random>
#include <cstdint>                                      // for uint32_t
#include <utility>
#include <unistd.h>

namespace {
	std::mutex global_seed_lock;
}

namespace siren {
namespace utilities {

    SIREN_random::SIREN_random() {
        // default to boring seed
        seed   = generate_seed();
        configuration       = std::default_random_engine(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0);
    }

    SIREN_random::SIREN_random(uint64_t _seed) {
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
    void SIREN_random::set_seed(uint64_t new_seed) {
        seed = new_seed;
        this->configuration = std::default_random_engine(seed);
    }

    uint64_t SIREN_random::get_seed() const {
        return seed;
    }

    uint64_t SIREN_random::generate_seed() {
        std::atomic_thread_fence(std::memory_order_acquire);
        std::lock_guard<std::mutex> lg(global_seed_lock);
        std::hash<std::string> string_hash;
        std::stringstream s;
        s << time(0) << getpid() << gethostid();
        std::atomic_thread_fence(std::memory_order_release);
        uint64_t seed = string_hash(s.str());
        return seed;
    }

} // namespace utilities
} // namespace siren
