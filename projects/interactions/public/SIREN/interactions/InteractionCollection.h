#pragma once
#ifndef SIREN_InteractionCollection_H
#define SIREN_InteractionCollection_H

#include <map>                                    // for map
#include <set>                                    // for set
#include <memory>                                 // for shared_ptr
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error

#include <map>
#include <set>
#include <memory>
#include <string>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"

namespace siren { namespace dataclasses { class InteractionRecord; } }

namespace siren {
namespace interactions {

class InteractionCollection {
private:
    siren::dataclasses::ParticleType primary_type;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<std::shared_ptr<Decay>> decays;
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>> cross_sections_by_target;
    std::set<siren::dataclasses::ParticleType> target_types;
    static const std::vector<std::shared_ptr<CrossSection>> empty;
    void InitializeTargetTypes();
public:
    InteractionCollection();
    virtual ~InteractionCollection() {};
    InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections);
    InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<Decay>> decays);
    InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::vector<std::shared_ptr<Decay>> decays);
    bool operator==(InteractionCollection const & other) const;
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSections() const {return cross_sections;}
    std::vector<std::shared_ptr<Decay>> const & GetDecays() const {return decays;}
    bool const HasCrossSections() const {return cross_sections.size() > 0;}
    bool const HasDecays() const {return decays.size() > 0;}
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSectionsForTarget(siren::dataclasses::ParticleType p) const;
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & GetCrossSectionsByTarget() const {
        return cross_sections_by_target;
    };
    std::set<siren::dataclasses::ParticleType> const & TargetTypes() const {
        return target_types;
    };
    double TotalDecayWidth(siren::dataclasses::InteractionRecord const & record) const;
    double TotalDecayLength(siren::dataclasses::InteractionRecord const & record) const;
    virtual bool MatchesPrimary(dataclasses::InteractionRecord const & record) const;
    std::map<siren::dataclasses::ParticleType, double> TotalCrossSectionByTarget(siren::dataclasses::InteractionRecord const & record) const;
    std::map<siren::dataclasses::ParticleType, double> TotalCrossSectionByTargetAllFinalStates(siren::dataclasses::InteractionRecord const & record) const;
public:
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::make_nvp("CrossSections", cross_sections));
            archive(cereal::make_nvp("Decays", decays));
        } else {
            throw std::runtime_error("InteractionCollection only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
            archive(cereal::make_nvp("TargetTypes", target_types));
            archive(cereal::make_nvp("CrossSections", cross_sections));
            archive(cereal::make_nvp("Decays", decays));
            InitializeTargetTypes();
        } else {
            throw std::runtime_error("InteractionCollection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::InteractionCollection, 0);

#endif // SIREN_InteractionCollection_H
