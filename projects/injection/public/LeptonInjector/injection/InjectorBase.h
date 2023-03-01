#pragma once
#ifndef LI_LeptonInjector_H
#define LI_LeptonInjector_H

#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/detector/EarthModel.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"

#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/dataclasses/DecayRecord.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/type/PrimaryInjector.h"

namespace LI {
namespace math {
class Vector3D;
}

namespace utilities {
class LI_random;
}

namespace injection {

class InjectorBase {
friend cereal::access;
public:
    virtual ~InjectorBase() {};
protected:
    unsigned int events_to_inject = 0;
    unsigned int injected_events = 0;
    std::shared_ptr<LI::utilities::LI_random> random;
    std::shared_ptr<LI::distributions::PrimaryInjector> primary_injector;
    std::shared_ptr<LI::crosssections::CrossSectionCollection> cross_sections;
    std::shared_ptr<LI::detector::EarthModel> earth_model;
    std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> distributions;
    InjectorBase();
public:
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<distributions::PrimaryInjector> primary_injector, std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections, std::shared_ptr<LI::detector::EarthModel> earth_model, std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> distributions, std::shared_ptr<LI::utilities::LI_random> random);
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<distributions::PrimaryInjector> primary_injector, std::vector<std::shared_ptr<LI::crosssections::CrossSection>> cross_sections, std::shared_ptr<LI::detector::EarthModel> earth_model, std::shared_ptr<LI::utilities::LI_random> random);
    InjectorBase(unsigned int events_to_inject, std::shared_ptr<LI::crosssections::CrossSectionCollection> cross_sections);
    virtual LI::dataclasses::InteractionRecord NewRecord() const;
    void SetRandom(std::shared_ptr<LI::utilities::LI_random> random);
    virtual void SampleCrossSection(LI::dataclasses::InteractionRecord & record) const;
    virtual void SampleSecondaryDecay(LI::dataclasses::InteractionRecord const & interaction, LI::dataclasses::DecayRecord & decay, double width, double alpha_gen, double alpha_phys, LI::geometry::Geometry *fiducial, double buffer) const;
    virtual void SamplePairProduction(LI::dataclasses::DecayRecord const & decay, LI::dataclasses::InteractionRecord & pairprod) const;
    LI::dataclasses::InteractionRecord GenerateEvent();
    virtual std::string Name() const;
    virtual double GenerationProbability(LI::dataclasses::InteractionRecord const & record) const;
    virtual std::set<std::vector<std::string>> DensityVariables() const;
    virtual std::pair<LI::math::Vector3D, LI::math::Vector3D> InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const;
    virtual std::vector<std::shared_ptr<LI::distributions::InjectionDistribution>> GetInjectionDistributions() const;
    virtual std::shared_ptr<LI::detector::EarthModel> GetEarthModel() const;
    virtual std::shared_ptr<LI::crosssections::CrossSectionCollection> GetCrossSections() const;
    unsigned int InjectedEvents() const;
    unsigned int EventsToInject() const;
    operator bool() const;

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("PrimaryInjector", primary_injector));
            archive(::cereal::make_nvp("CrossSections", cross_sections));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("InjectionDistributions", distributions));
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("EventsToInject", events_to_inject));
            archive(::cereal::make_nvp("InjectedEvents", injected_events));
            archive(::cereal::make_nvp("PrimaryInjector", primary_injector));
            archive(::cereal::make_nvp("CrossSections", cross_sections));
            archive(::cereal::make_nvp("EarthModel", earth_model));
            archive(::cereal::make_nvp("InjectionDistributions", distributions));
        } else {
            throw std::runtime_error("InjectorBase only supports version <= 0!");
        }
    }
};

} // namespace injection
} // namespace LI

CEREAL_CLASS_VERSION(LI::injection::InjectorBase, 0);

#endif // LI_LeptonInjector_H

