#pragma once
#ifndef SIREN_Hadronization_H
#define SIREN_Hadronization_H

#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/dataclasses/InteractionRecord.h" // for InteractionSignature

#include "SIREN/utilities/Random.h" // for SIREN_random
#include "SIREN/geometry/Geometry.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class Hadronization {
friend cereal::access;
private:
public:
    Hadronization();
    virtual ~Hadronization() {};
    bool operator==(Hadronization const & other) const;
    virtual bool equal(Hadronization const & other) const = 0;
    
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const = 0;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const = 0;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const = 0;
    virtual double FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const = 0;

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {};
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {};

}; // class Hadronization

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::Hadronization, 0);


#endif // SIREN_Hadronization_H
