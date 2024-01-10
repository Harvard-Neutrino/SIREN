#pragma once
#ifndef LI_DecayRangeLeptonInjector_H
#define LI_DecayRangeLeptonInjector_H

#include <memory>
#include <string>
#include <vector>
#include <cstdint>                                  // for uint32_t
#include <utility>
#include <stdexcept>                                // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/Decay.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/DecayRangePositionDistribution.h"
#include "LeptonInjector/injection/InjectorBase.h"  // for InjectorBase

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace injection { struct InjectionProcess; } }
namespace LI { namespace math { class Vector3D; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace injection {

class DecayRangeLeptonInjector : public InjectorBase {
friend cereal::access;
protected:
    std::shared_ptr<LI::distributions::DecayRangeFunction> range_func;
    double disk_radius;
    double endcap_length;
    std::shared_ptr<LI::distributions::DecayRangePositionDistribution> position_distribution;
    std::shared_ptr<LI::interactions::InteractionCollection> cross_sections;
    DecayRangeLeptonInjector();
public:
    DecayRangeLeptonInjector(unsigned int events_to_inject, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes, std::shared_ptr<LI::utilities::LI_random> random, std::shared_ptr<LI::distributions::DecayRangeFunction> range_func, double disk_radius, double endcap_length);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("RangeFunction", range_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<InjectorBase>(this));
        } else {
            throw std::runtime_error("DecayRangeLeptonInjector only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("RangeFunction", range_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<InjectorBase>(this));
        } else {
            throw std::runtime_error("DecayRangeLeptonInjector only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::DecayRangeLeptonInjector, 0);
CEREAL_REGISTER_TYPE(LI::injection::DecayRangeLeptonInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::injection::InjectorBase, LI::injection::DecayRangeLeptonInjector);

#endif // LI_DecayRangeLeptonInjector_H
