#pragma once
#ifndef LI_ColumnDepthSIREN_H
#define LI_ColumnDepthSIREN_H

#include <tuple>
#include <memory>
#include <string>
#include <vector>                                   // for vector
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

#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "SIREN/distributions/primary/vertex/DepthFunction.h"
#include "SIREN/injection/Injector.h"  // for Injector

namespace SI { namespace interactions { class InteractionCollection; } }
namespace SI { namespace dataclasses { class InteractionRecord; } }
namespace SI { namespace injection { class PrimaryInjectionProcess; } }
namespace SI { namespace math { class Vector3D; } }
namespace SI { namespace utilities { class LI_random; } }

namespace SI {
namespace injection {

class ColumnDepthSIREN : public Injector {
friend cereal::access;
protected:
    std::shared_ptr<SI::distributions::DepthFunction> depth_func;
    double disk_radius;
    double endcap_length;
    std::shared_ptr<SI::distributions::ColumnDepthPositionDistribution> position_distribution;
    std::shared_ptr<SI::interactions::InteractionCollection> interactions;
    ColumnDepthSIREN();
public:
    ColumnDepthSIREN(unsigned int events_to_inject, std::shared_ptr<SI::detector::DetectorModel> detector_model, std::shared_ptr<injection::PrimaryInjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<SI::utilities::LI_random> random, std::shared_ptr<SI::distributions::DepthFunction> depth_func, double disk_radius, double endcap_length);
    std::string Name() const override;
    virtual std::tuple<SI::math::Vector3D, SI::math::Vector3D> PrimaryInjectionBounds(SI::dataclasses::InteractionRecord const & interaction) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("DepthFunction", depth_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("ColumnDepthSIREN only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DepthFunction", depth_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("ColumnDepthSIREN only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace SI

CEREAL_CLASS_VERSION(SI::injection::ColumnDepthSIREN, 0);
CEREAL_REGISTER_TYPE(SI::injection::ColumnDepthSIREN);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::injection::Injector, SI::injection::ColumnDepthSIREN);

#endif // LI_ColumnDepthSIREN_H
