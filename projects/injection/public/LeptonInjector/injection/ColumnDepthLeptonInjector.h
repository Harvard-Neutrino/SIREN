#pragma once
#ifndef LI_ColumnDepthLeptonInjector_H
#define LI_ColumnDepthLeptonInjector_H

#include <memory>
#include <string>
#include <vector>                                   // for vector
#include <utility>
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

#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/Decay.h"
#include "LeptonInjector/detector/DetectorModel.h"
#include "LeptonInjector/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "LeptonInjector/distributions/primary/vertex/DepthFunction.h"
#include "LeptonInjector/injection/Injector.h"  // for Injector

namespace LI { namespace interactions { class InteractionCollection; } }
namespace LI { namespace dataclasses { struct InteractionRecord; } }
namespace LI { namespace injection { struct InjectionProcess; } }
namespace LI { namespace math { class Vector3D; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace injection {

class ColumnDepthLeptonInjector : public Injector {
friend cereal::access;
protected:
    std::shared_ptr<LI::distributions::DepthFunction> depth_func;
    double disk_radius;
    double endcap_length;
    std::shared_ptr<LI::distributions::ColumnDepthPositionDistribution> position_distribution;
    std::shared_ptr<LI::interactions::InteractionCollection> cross_sections;
    ColumnDepthLeptonInjector();
public:
    ColumnDepthLeptonInjector(unsigned int events_to_inject, std::shared_ptr<LI::detector::DetectorModel> earth_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes, std::shared_ptr<LI::utilities::LI_random> random, std::shared_ptr<LI::distributions::DepthFunction> depth_func, double disk_radius, double endcap_length);
    std::string Name() const override;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("DepthFunction", depth_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("ColumnDepthLeptonInjector only supports version <= 0!");
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
            throw std::runtime_error("ColumnDepthLeptonInjector only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::ColumnDepthLeptonInjector, 0);
CEREAL_REGISTER_TYPE(LI::injection::ColumnDepthLeptonInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::injection::Injector, LI::injection::ColumnDepthLeptonInjector);

#endif // LI_ColumnDepthLeptonInjector_H
