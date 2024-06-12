#pragma once
#ifndef SIREN_CharmHadronization_H
#define SIREN_CharmHadronization_H

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

#include "SIREN/interactions/Hadronization.h"   // for Hadronization
#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/dataclasses/InteractionSignature.h" // for InteractionSignature
#include "SIREN/dataclasses/InteractionRecord.h" // for InteractionSignature

#include "SIREN/utilities/Random.h" // for SIREN_random
#include "SIREN/geometry/Geometry.h"
#include "SIREN/utilities/Constants.h"            // for electronMass


namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class CharmHadronization : public Hadronization {
friend cereal::access;
private:
    const std::set<siren::dataclasses::Particle::ParticleType> primary_types = {siren::dataclasses::Particle::ParticleType::Charm, siren::dataclasses::Particle::ParticleType::CharmBar};
    
public:
    
    CharmHadronization();

    // virtual pybind11::object get_self();
    
    virtual bool equal(Hadronization const & other) const override;
    
    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & interaction, std::shared_ptr<siren::utilities::SIREN_random> random) const override;
    
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override; // Requires python-side implementation
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::Particle::ParticleType primary) const override; // Requires python-side implementation
    
    double FragmentationFraction(siren::dataclasses::Particle::ParticleType secondary) const override;

    static double getHadronMass(siren::dataclasses::ParticleType hadron_type);

public:
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Hadronization", cereal::virtual_base_class<Hadronization>(this)));
        } else {
            throw std::runtime_error("CharmHadronization only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<CharmHadronization> & construct, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Hadronization", cereal::virtual_base_class<Hadronization>(construct.ptr())));
        } else {
            throw std::runtime_error("CharmHadronization only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::CharmHadronization, 0);


#endif // SIREN_CharmHadronization_H
