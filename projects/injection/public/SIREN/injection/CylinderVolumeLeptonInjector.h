#pragma once
#ifndef LI_CylinderVolumeSIREN_H
#define LI_CylinderVolumeSIREN_H

#include <tuple>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>                                  // for uint32_t
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

#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "SIREN/geometry/Cylinder.h"       // for Cylinder
#include "SIREN/injection/Injector.h"  // for Injector

namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace injection { class PrimaryInjectionProcess; } }
namespace SI { namespace math { class Vector3D; } }  // lines 40-40
namespace SI { namespace utilities { class LI_random; } }  // lines 44-44

namespace SI {
namespace math {
class Vector3D;
}

namespace utilities {
class LI_random;
}

namespace injection {

class CylinderVolumeSIREN : public Injector {
friend cereal::access;
protected:
    std::shared_ptr<SI::distributions::CylinderVolumePositionDistribution> position_distribution;
    std::shared_ptr<SI::interactions::InteractionCollection> interactions;
    CylinderVolumeSIREN();
public:
    CylinderVolumeSIREN(unsigned int events_to_inject, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<SI::utilities::LI_random> random, SI::geometry::Cylinder cylinder);
    std::string Name() const override;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const override;
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("CylinderVolumeSIREN only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("CylinderVolumeSIREN only supports version <= 0!");
        }
    }
};


} // namespace injection
} // namespace SI

CEREAL_CLASS_VERSION(SI::injection::CylinderVolumeSIREN, 0);
CEREAL_REGISTER_TYPE(SI::injection::CylinderVolumeSIREN);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::injection::Injector, SI::injection::CylinderVolumeSIREN);

#endif // LI_CylinderVolumeSIREN_H
