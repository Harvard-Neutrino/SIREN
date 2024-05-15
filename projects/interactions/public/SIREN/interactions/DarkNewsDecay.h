#pragma once
#ifndef SIREN_DarkNewsDecay_H
#define SIREN_DarkNewsDecay_H

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
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include "SIREN/interactions/Decay.h"   // for Decay
#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

namespace siren { namespace dataclasses { class InteractionRecord; } }
namespace siren { namespace dataclasses { class CrossSectionDistributionRecord; } }
namespace siren { namespace dataclasses { struct InteractionSignature; } }
namespace siren { namespace utilities { class SIREN_random; } }

namespace siren {
namespace interactions {

class DarkNewsDecay : public Decay {
friend cereal::access;
public:
    
    DarkNewsDecay();
    
    virtual bool equal(Decay const & other) const override;
    
    virtual double TotalDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override;
    virtual double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const &) const override;
    virtual double DifferentialDecayWidth(dataclasses::InteractionRecord const &) const override;
    virtual void SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const;
    virtual void SampleFinalState(dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<siren::utilities::SIREN_random>) const override;
    
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override = 0; // Requires python-side implementation
    virtual std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary) const override = 0; // Requires python-side implementation
    
    virtual double FinalStateProbability(dataclasses::InteractionRecord const & record) const override;
public:
    virtual std::vector<std::string> DensityVariables() const override = 0; // Requires python-side implementation
    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version == 0) {
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(this)));
        } else {
            throw std::runtime_error("DarkNewsDecay only supports version <= 0!");
        }
    }
    template<typename Archive>
    void load_and_construct(Archive & archive, cereal::construct<DarkNewsDecay> & construct, std::uint32_t version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Decay", cereal::virtual_base_class<Decay>(construct.ptr())));
        } else {
            throw std::runtime_error("DarkNewsDecay only supports version <= 0!");
        }
    }

}; // class DarkNewsDecay

// Trampoline class for DarkNewsDecay
class pyDarkNewsDecay : public DarkNewsDecay,public Pybind11Trampoline<DarkNewsDecay, pyDarkNewsDecay> {
public:
    using DarkNewsDecay::DarkNewsDecay;
    pyDarkNewsDecay(DarkNewsDecay && parent) : DarkNewsDecay(std::move(parent)) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsDecay)))));
    }
    pyDarkNewsDecay(DarkNewsDecay const & parent) : DarkNewsDecay(parent) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsDecay)))));
    }
    pybind11::object self;

    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            std::cref(interaction)
        )
    }

    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            std::cref(interaction)
        )
    }

    double TotalDecayWidth(siren::dataclasses::ParticleType primary) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            std::cref(interaction)
        )
    }

    void SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleRecordFromDarkNews,
            "SampleRecordFromDarkNews",
            std::ref(record),
            random
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParent,
            "GetPossibleSignaturesFromParent",
            primary_type
        )
    }

    std::vector<std::string> DensityVariables() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<std::string>,
            DensityVariables,
            "DensityVariables"
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
        )
    }

    template<typename Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        archive(cereal::base_class<Pybind11Trampoline<DarkNewsDecay, pyDarkNewsDecay>>(this));
    }

    template<typename Archive>
    void load(Archive & archive, std::uint32_t version) {
        archive(cereal::base_class<Pybind11Trampoline<DarkNewsDecay, pyDarkNewsDecay>>(this));
    }
}; // class pyDarkNewsDecay

} // namespace interactions
} // namespace siren

CEREAL_CLASS_VERSION(siren::interactions::DarkNewsDecay, 0);
CEREAL_REGISTER_TYPE(siren::interactions::DarkNewsDecay);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::interactions::Decay, siren::interactions::DarkNewsDecay);

using DarkNewsDecayTrampolineType = Pybind11Trampoline<siren::interactions::DarkNewsDecay,siren::interactions::pyDarkNewsDecay>;
RegisterTrampolineCerealMethods(siren::interactions::DarkNewsDecay, siren::interactions::pyDarkNewsDecay, DarkNewsDecayTrampolineType)

#endif // SIREN_DarkNewsDecay_H
