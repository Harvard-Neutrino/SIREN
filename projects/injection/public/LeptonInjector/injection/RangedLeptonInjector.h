#pragma once
#ifndef LI_RangedLeptonInjector_H
#define LI_RangedLeptonInjector_H

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
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/interactions/InteractionCollection.h"
#include "LeptonInjector/interactions/CrossSection.h"
#include "LeptonInjector/interactions/Decay.h"
#include "LeptonInjector/detector/DetectorModel.h"
#include "LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "LeptonInjector/distributions/primary/vertex/RangePositionDistribution.h"
#include "LeptonInjector/injection/Injector.h"  // for Injector

namespace LI { namespace dataclasses { class InteractionRecord; } }
namespace LI { namespace injection { class InjectionProcess; } }
namespace LI { namespace math { class Vector3D; } }
namespace LI { namespace utilities { class LI_random; } }

namespace LI {
namespace injection {

class RangedLeptonInjector : public Injector {
friend cereal::access;
protected:
    std::shared_ptr<LI::distributions::RangeFunction> range_func;
    double disk_radius;
    double endcap_length;
    std::shared_ptr<LI::distributions::RangePositionDistribution> position_distribution;
    std::shared_ptr<LI::interactions::InteractionCollection> interactions;
    RangedLeptonInjector();
public:
    RangedLeptonInjector(unsigned int events_to_inject, std::shared_ptr<LI::detector::DetectorModel> detector_model, std::shared_ptr<injection::InjectionProcess> primary_process, std::vector<std::shared_ptr<injection::SecondaryInjectionProcess>> secondary_processes, std::shared_ptr<LI::utilities::LI_random> random, std::shared_ptr<LI::distributions::RangeFunction> range_func, double disk_radius, double endcap_length);
    std::string Name() const override;
    virtual std::tuple<LI::math::Vector3D, LI::math::Vector3D> PrimaryInjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const override;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("RangeFunction", range_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("RangedLeptonInjector only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("RangeFunction", range_func));
            archive(::cereal::make_nvp("DiskRadius", disk_radius));
            archive(::cereal::make_nvp("EndcapLength", endcap_length));
            archive(::cereal::make_nvp("PositionDistribution", position_distribution));
            archive(cereal::virtual_base_class<Injector>(this));
        } else {
            throw std::runtime_error("RangedLeptonInjector only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::RangedLeptonInjector, 0);
CEREAL_REGISTER_TYPE(LI::injection::RangedLeptonInjector);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::injection::Injector, LI::injection::RangedLeptonInjector);

#endif // LI_RangedLeptonInjector_H
