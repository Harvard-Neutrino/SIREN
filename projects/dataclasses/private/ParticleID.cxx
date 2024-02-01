#include "LeptonInjector/dataclasses/ParticleID.h"

#include <tuple>
#include <mutex>
#include <atomic>
#include <unistd.h>

namespace LI {
namespace dataclasses {

namespace {
	std::mutex global_id_lock;
	std::atomic<int32_t> global_last_pid_(0);
	std::atomic<int32_t> global_minor_id_(0);
	std::atomic<uint64_t> global_major_id_(0);
}

bool ParticleID::operator<(ParticleID const & other) const {
    return std::tie(major_id, minor_id) < std::tie(other.major_id, other.minor_id);
}

bool ParticleID::operator==(ParticleID const & other) const {
    return std::tie(major_id, minor_id) == std::tie(other.major_id, other.minor_id);
}

// Adapted from https://github.com/icecube/icetray-public/blob/4436c3e10c23f95a8965c98fecccb7775a361fab/dataclasses/private/dataclasses/physics/I3Particle.cxx#L42-L93
ParticleID ParticleID::GenerateID() {
    int this_pid = getpid();

    int32_t last_pid = global_last_pid_.load(std::memory_order_relaxed);
    std::atomic_thread_fence(std::memory_order_acquire); //keep memory ops from wandering
    if (this_pid != last_pid) {
        std::lock_guard<std::mutex> lg(global_id_lock); //acquire the lock
                                                            //check whether another thread already updated this
        last_pid = global_last_pid_.load(std::memory_order_relaxed);
        if(this_pid != last_pid){
            std::atomic_thread_fence(std::memory_order_release);
            global_last_pid_.store(this_pid, std::memory_order_relaxed);
            global_major_id_.store(0, std::memory_order_relaxed); // this will cause a new major ID to be generated
            global_minor_id_.store(0, std::memory_order_relaxed); // reset the minor ID, too
        }
    }

    uint64_t old_major_id = global_major_id_.load(std::memory_order_relaxed);
    std::atomic_thread_fence(std::memory_order_acquire); //keep memory ops from wandering
    if(old_major_id == 0){
        std::lock_guard<std::mutex> lg(global_id_lock); //acquire the lock
        //check whether another thread already updated this
        old_major_id = global_major_id_.load(std::memory_order_relaxed);
        if(old_major_id == 0){
            std::hash<std::string> string_hash;
            std::stringstream s;
            s << time(0) << this_pid<<gethostid();
            std::atomic_thread_fence(std::memory_order_release);
            global_major_id_.store(string_hash(s.str()), std::memory_order_relaxed);
        }
    }
    ParticleID ID;

    ID.major_id = global_major_id_.load();
    ID.minor_id = global_minor_id_.fetch_add(1);

    return ID;
}

} // namespace dataclasses
} // namespace LI
