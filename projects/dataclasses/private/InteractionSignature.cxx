#include <tuple>

#include "LeptonInjector/dataclasses/InteractionSignature.h"

namespace LI {
namespace dataclasses {

bool InteractionSignature::operator==(InteractionSignature const & other) const {
    return
        std::tie(primary_type, target_type, secondary_types)
        ==
        std::tie(other.primary_type, other.target_type, other.secondary_types);
    /*
    if(primary_type != other.primary_type or target_type != other.target_type) {
        return false;
    } else {
        std::map<LI::dataclasses::Particle::ParticleType, int> m0;
        for(auto p : secondary_types) {
            auto it = m0.find(p);
            if(it == m0.end()) {
                m0.insert({p, 1});
            } else {
                it->second += 1;
            }
        }
        std::map<LI::dataclasses::Particle::ParticleType, int> m1;
        for(auto p : other.secondary_types) {
            auto it = m1.find(p);
            if(it == m1.end()) {
                m1.insert({p, 1});
            } else {
                it->second += 1;
            }
        }
        return m0 == m1;
    }
    */
}

bool InteractionSignature::operator<(InteractionSignature const & other) const {
    return
        std::tie(primary_type, target_type, secondary_types)
        <
        std::tie(other.primary_type, other.target_type, other.secondary_types);
}

} // namespace dataclasses
} // namespace LI

std::ostream& operator<<(std::ostream& os, LI::dataclasses::InteractionSignature const& signature) {
    std::stringstream ss;
    ss << "InteractionSignature (" << &signature << ") ";
    os << ss.str() << '\n';


    os << "PrimaryType: " << signature.primary_type << "\n";
    os << "TargetType: " << signature.target_type << "\n";
    os << "SecondaryTypes:";
    for(auto secondary: signature.secondary_types) {
        os << " " << secondary;
    }
    os << std::endl;

    return os;
}

