#pragma once
#ifndef LI_InteractionCollection_H
#define LI_InteractionCollection_H

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
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/dataclasses/Particle.h"  // for Particle

namespace SI { namespace interactions { class CrossSection; } }
namespace SI { namespace interactions { class Decay; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }

namespace SI {
namespace interactions {

class InteractionCollection {
private:
    SI::dataclasses::ParticleType primary_type;
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<std::shared_ptr<Decay>> decays;
    std::map<SI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>> cross_sections_by_target;
    std::set<SI::dataclasses::ParticleType> target_types;
    static const std::vector<std::shared_ptr<CrossSection>> empty;
    void InitializeTargetTypes();
public:
    InteractionCollection();
    virtual ~InteractionCollection() {};
    InteractionCollection(SI::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections);
    InteractionCollection(SI::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<Decay>> decays);
    InteractionCollection(SI::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::vector<std::shared_ptr<Decay>> decays);
    bool operator==(InteractionCollection const & other) const;
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSections() const {return cross_sections;}
    std::vector<std::shared_ptr<Decay>> const & GetDecays() const {return decays;}
    bool const HasCrossSections() const {return cross_sections.size() > 0;}
    bool const HasDecays() const {return decays.size() > 0;}
    std::vector<std::shared_ptr<CrossSection>> const & GetCrossSectionsForTarget(SI::dataclasses::ParticleType p) const;
    std::map<SI::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & GetCrossSectionsByTarget() const {
        return cross_sections_by_target;
    };
    std::set<SI::dataclasses::ParticleType> const & TargetTypes() const {
        return target_types;
    };
    double TotalDecayWidth(SI::dataclasses::InteractionRecord const & record) const;
    double TotalDecayLength(SI::dataclasses::InteractionRecord const & record) const;
    virtual bool MatchesPrimary(dataclasses::InteractionRecord const & record) const;
public:
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(cereal::make_nvp("PrimaryType", primary_type));
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
            archive(cereal::make_nvp("CrossSections", cross_sections));
            archive(cereal::make_nvp("Decays", decays));
        } else {
            throw std::runtime_error("InteractionCollection only supports version <= 0!");
        }
    }
};

} // namespace interactions
} // namespace SI

CEREAL_CLASS_VERSION(SI::interactions::InteractionCollection, 0);

#endif // LI_InteractionCollection_H
