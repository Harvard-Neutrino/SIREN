#include "SIREN/utilities/Random.h"

#include <mutex>
#include <atomic>
#include <random>
#include <cstdint>                                      // for uint64_t
#include <utility>
#include <algorithm>
#include <charconv>
#include <locale>
#include <sstream>
#include <vector>
#include <unistd.h>

namespace {
	std::mutex global_seed_lock;

    std::size_t EngineStateWordCount() {
#if defined(_LIBCPP_VERSION)
        return std::mt19937_64::state_size;
#elif defined(__GLIBCXX__)
        return std::mt19937_64::state_size + 1;
#else
        throw std::runtime_error("SIREN_random: unsupported C++ standard library for RNG archives");
#endif
    }
}

namespace siren {
namespace utilities {

    std::string SIREN_random::EngineStateFormat() {
#if defined(_LIBCPP_VERSION)
        return "mt19937_64/libc++/1";
#elif defined(__GLIBCXX__)
        return "mt19937_64/libstdc++/1";
#else
        throw std::runtime_error("SIREN_random: unsupported C++ standard library for RNG archives");
#endif
    }

    std::string SIREN_random::WriteEngineState() const {
        std::ostringstream stream;
        stream.imbue(std::locale::classic());
        stream << configuration;
        return stream.str();
    }

    std::mt19937_64 SIREN_random::ReadEngineState(std::string const & text) {
        auto invalid = []() {
            return std::runtime_error("SIREN_random: invalid or incompatible engine state in archive; "
                "legacy RNG archives require the same C++ standard library");
        };
        std::istringstream tokens(text);
        tokens.imbue(std::locale::classic());
        std::vector<uint64_t> words;
        std::string token;
        auto count = EngineStateWordCount();
        while(tokens >> token) {
            uint64_t value;
            auto result = std::from_chars(token.data(), token.data() + token.size(), value);
            if(result.ec != std::errc{} || result.ptr != token.data() + token.size() || words.size() == count)
                throw invalid();
            words.push_back(value);
        }
        if(words.size() != count ||
                std::all_of(words.begin(), words.begin() + std::mt19937_64::state_size,
                    [](uint64_t value) { return value == 0; }))
            throw invalid();
        if(count > std::mt19937_64::state_size && words.back() > std::mt19937_64::state_size)
            throw invalid();

        std::istringstream stream(text);
        stream.imbue(std::locale::classic());
        std::mt19937_64 engine;
        if(!(stream >> engine)) throw invalid();
        stream >> std::ws;
        if(!stream.eof()) throw invalid();
        return engine;
    }

    SIREN_random::SIREN_random() {
        // default to a generated seed (hashed from time/pid/host)
        seed   = generate_seed();
        configuration       = std::mt19937_64(seed);
        generator           = std::uniform_real_distribution<double>( 0.0, 1.0);
    }

    SIREN_random::SIREN_random(uint64_t _seed) {
        seed = _seed;
        configuration       = std::mt19937_64(seed);
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
        this->configuration = std::mt19937_64(seed);
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
