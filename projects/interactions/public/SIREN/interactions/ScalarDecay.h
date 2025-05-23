#pragma once
#ifndef SIREN_ScalarDecay_H
#define SIREN_ScalarDecay_H

#include <set>                                    // for set
#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/interactions/Decay.h"   // for Decay
#include "SIREN/dataclasses/Particle.h"  // for Particle

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class ScalarDecay : public Decay {
friend cereal::access;
protected:
ScalarDecay() {};
private:
    double mass;
    double g; // coupling
    const std::set<siren::dataclasses::ParticleType> primary_types = {siren::dataclasses::ParticleType::Scalar};
public:
    ScalarDecay(double mass, double coupling, ChiralNature nature) : mass(mass), coupling(coupling) {};
    virtual bool equal(Decay const & other) const override;
    double GetMass() const {return mass;};
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override;
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const override;
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
public:
    virtual std::vector<std::string> DensityVariables() const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PrimaryTypes", primary_types));
            archive(::cereal::make_nvp("Mass", mass));
            archive(::cereal::make_nvp("Coupling", coupling));
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("ScalarDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<ScalarDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            std::set<siren::dataclasses::ParticleType> _primary_types;
            double _mass;
            double _coupling;

            archive(::cereal::make_nvp("PrimaryTypes", _primary_types));
            archive(::cereal::make_nvp("Mass", _mass));
            archive(::cereal::make_nvp("Coupling", _coupling));
            construct(_mass, _coupling, _primary_types);
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("ScalarDecay only supports version <= 0!");
        }
    }

};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::ScalarDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::ScalarDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::ScalarDecay);

#endif // SIREN_ScalarDecay_H
