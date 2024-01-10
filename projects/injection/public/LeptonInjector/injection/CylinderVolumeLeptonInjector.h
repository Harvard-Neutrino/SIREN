#pragma once
#ifndef LI_CylinderVolumeLeptonInjector_H
#define LI_CylinderVolumeLeptonInjector_H

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

#include "LeptonInjector/crosssections/InteractionCollection.h"
#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/Decay.h"
#include "LeptonInjector/detector/EarthModel.h"
#include "LeptonInjector/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "LeptonInjector/geometry/Cylinder.h"       // for Cylinder
#include "LeptonInjector/injection/InjectorBase.h"  // for InjectorBase

namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace injection { struct InjectionProcess; } }
namespace LI { namespace math { class Vector3D; } }  // lines 40-40
namespace LI { namespace utilities { class LI_random; } }  // lines 44-44

namespace LI {
namespace math {
class Vector3D;
}

namespace utilities {
class LI_random;
}

namespace injection {

class CylinderVolumeLeptonInjector : public InjectorBase {
friend cereal::access;
protected:
    std::shared_ptr<LI::distributions::CylinderVolumePositionDistribution> position_distribution;
    std::shared_ptr<LI::crosssections::InteractionCollection> cross_sections;
    CylinderVolumeLeptonInjector();
public:
    CylinderVolumeLeptonInjector(unsigned int events_to_inject, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes, std::shared_ptr<LI::utilities::LI_random> random, LI::geometry::Cylinder cylinder);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<InjectorBase>(this));
        } else {
            throw std::runtime_error("CylinderVolumeLeptonInjector only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<InjectorBase>(this));
        } else {
            throw std::runtime_error("CylinderVolumeLeptonInjector only supports version <= 0!");
        }
    }
};


} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::CylinderVolumeLeptonInjector, 0);
CEREAL_REGISTER_TYPE(LI::injection::CylinderVolumeLeptonInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::injection::InjectorBase, LI::injection::CylinderVolumeLeptonInjector);

#endif // LI_CylinderVolumeLeptonInjector_H
