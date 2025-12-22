#include "SIREN/dataclasses/ParticleID.h"

#include <tuple>
#include <mutex>
#include <atomic>
#include <ostream>
#include <unistd.h>

#include "SIREN/utilities/StringManipulation.h"

std::ostream& operator<<(std::ostream& os, siren::dataclasses::ParticleID const & id) {
    os << to_repr(id);
    return os;
}

std::string to_str(siren::dataclasses::ParticleID const & id) {
    using siren::utilities::tab;
    std::stringstream ss;
    ss << "[ ParticleID (" << &id << ")\n";
    ss << tab << "IDSet: " << id.IsSet() << '\n';
    ss << tab << "MajorID: " << id.GetMajorID() << '\n';
    ss << tab << "MinorID: " << id.GetMinorID() << '\n';
    ss << ']';
    return ss.str();
}

std::string to_repr(siren::dataclasses::ParticleID const & id) {
    std::stringstream ss;
    ss << "ParticleID(";
    if(id.IsSet())
        ss << id.GetMajorID() << ", " << id.GetMinorID();
    else
        ss << "unset";
    ss << ")";
    return ss.str();
}

namespace siren {
namespace dataclasses {

namespace {
	std::mutex global_id_lock;
	std::atomic<int32_t> global_last_pid_(0);
	std::atomic<int32_t> global_minor_id_(0);
	std::atomic<uint64_t> global_major_id_(0);
}

ParticleID::ParticleID() : id_set(false), major_id(0), minor_id(0) {};

ParticleID::ParticleID(uint64_t major, int32_t minor) :
    id_set(true), major_id(major), minor_id(minor) {};

bool ParticleID::IsSet() const {
    return id_set;
}

ParticleID::operator bool() const {
    return id_set;
}

uint64_t ParticleID::GetMajorID() const {
    return major_id;
}

int32_t ParticleID::GetMinorID() const {
    return minor_id;
}

void ParticleID::SetID(uint64_t major, int32_t minor) {
    id_set = true;
    major_id = major;
    minor_id = minor;
}

bool ParticleID::operator<(ParticleID const & other) const {
    return std::tie(id_set, major_id, minor_id) < std::tie(id_set, other.major_id, other.minor_id);
}

bool ParticleID::operator==(ParticleID const & other) const {
    return std::tie(id_set, major_id, minor_id) == std::tie(id_set, other.major_id, other.minor_id);
}

bool ParticleID::operator!=(ParticleID const & other) const {
    return not (*this == other);
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
        if(old_major_id == 0) {
            std::hash<std::string> string_hash;
            std::stringstream s;
            s << time(0) << this_pid<<gethostid();
            std::atomic_thread_fence(std::memory_order_release);
            global_major_id_.store(string_hash(s.str()), std::memory_order_relaxed);
        }
    }
    ParticleID ID;

    ID.id_set = true;
    ID.major_id = global_major_id_.load();
    ID.minor_id = global_minor_id_.fetch_add(1);

    return ID;
}

} // namespace dataclasses
} // namespace siren
