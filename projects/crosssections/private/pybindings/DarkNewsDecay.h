#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/DarkNewsDecay.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

// Macro for defining pure virtual methods of PyDarkNewsDecay 
#define C_PYBIND11_OVERRIDE_PURE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        const BaseType * ref; \
        if(selfname) { \
            ref = selfname.cast<BaseType *>(); \
        } else { \
            ref = this; \
        } \
        do { \
            do { \
                auto *tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                pybind11::function override = \
                    tinfo ? pybind11::detail::get_type_override(static_cast<const BaseType *>(ref), tinfo, pyfuncname) : pybind11::function(); \
                if (override) { \
                    auto o = override(__VA_ARGS__); \
                    if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                        static pybind11::detail::override_caster_t<returnType> caster; \
                        return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                    } \
                    return pybind11::detail::cast_safe<returnType>(std::move(o)); \
                } \
            } while (false); \
            pybind11::pybind11_fail( \
                "Tried to call pure virtual function \"" PYBIND11_STRINGIFY(BaseType) "::" "cfuncname" "\""); \
        } while (false);

// Macro for defining virtual methods of PyDarkNewsDecay 
#define C_PYBIND11_OVERRIDE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        const BaseType * ref; \
        if(selfname) { \
            ref = selfname.cast<BaseType *>(); \
        } else { \
            ref = this; \
        } \
        do { \
            do { \
                auto *tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                pybind11::function override = \
                    tinfo ? pybind11::detail::get_type_override(static_cast<const BaseType *>(ref), tinfo, pyfuncname) : pybind11::function(); \
                if (override) { \
                    auto o = override(__VA_ARGS__); \
                    if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                        static pybind11::detail::override_caster_t<returnType> caster; \
                        return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                    } \
                    return pybind11::detail::cast_safe<returnType>(std::move(o)); \
                } \
            } while (false); \
            return BaseType::cfuncname(__VA_ARGS__); \
        } while (false);

namespace LI {
namespace crosssections {
// Trampoline class for CrossSection
class pyDecay : public Decay {
public:
    using Decay::Decay;
    pyDecay(Decay && parent) : Decay(std::move(parent)) {}
    pybind11::object self;

    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            interaction
        )
    }

    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            interaction
        )
    }

    double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            interaction
        )
    }

    void SampleFinalState(dataclasses::InteractionRecord & interaction, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            void,
            SampleFinalState,
            "SampleFinalState",
            interaction,
            random
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            record
        )
    }

    pybind11::object get_self() {
        return self;
    }
};
// Trampoline class for DarkNewsDecay
class pyDarkNewsDecay : public DarkNewsDecay {
public:
    using DarkNewsDecay::DarkNewsDecay;
    pyDarkNewsDecay(DarkNewsDecay && parent) : DarkNewsDecay(std::move(parent)) {}
    pybind11::object self;

    double TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            interaction
        )
    }

    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            interaction
        )
    }

    double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            interaction
        )
    }

    void SampleFinalState(dataclasses::InteractionRecord & interaction, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            void,
            SampleFinalState,
            "SampleFinalState",
            interaction,
            random
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            record
        )
    }

    pybind11::object get_self() override {
        return self;
    }
};
} // end crosssections namespace
} // end LI namespace

void register_DarkNewsDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    // Bindings for pyDarkNewsDecay
    class_<LI::crosssections::pyDarkNewsDecay> pyDarkNewsDecay(m, "pyDarkNewsDecay");

    pyDarkNewsDecay
        .def(init<>())
        .def("__eq__", [](const LI::crosssections::DarkNewsDecay &self, const LI::crosssections::DarkNewsDecay &other){ return self == other; })
        .def("equal", &LI::crosssections::DarkNewsDecay::equal)
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&DarkNewsDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&DarkNewsDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&DarkNewsDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&DarkNewsDecay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&DarkNewsDecay::DensityVariables)
        .def("FinalStateProbability",&DarkNewsDecay::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsDecay::SampleFinalState)
        .def("get_self", &pyDarkNewsDecay::get_self)
        .def(pybind11::pickle(
            [](const LI::crosssections::pyDarkNewsDecay & cpp_obj) {
                pybind11::object self;
                if(cpp_obj.self) {
                    self = pybind11::reinterpret_borrow<pybind11::object>(cpp_obj.self);
                } else {
                    auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsDecay));
                    pybind11::handle self_handle = get_object_handle(static_cast<const DarkNewsDecay *>(&cpp_obj), tinfo);
                    self = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
                }
                pybind11::dict d;
                if (pybind11::hasattr(self, "__dict__")) {
                    d = self.attr("__dict__");
                }
                return pybind11::make_tuple(d);
            },
            [](const pybind11::tuple &t) {
                if (t.size() != 1) {
                    throw std::runtime_error("Invalid state!");
                }
                auto cpp_state = std::unique_ptr<LI::crosssections::pyDarkNewsDecay>(new LI::crosssections::pyDarkNewsDecay);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;


    class_<DarkNewsDecay, std::shared_ptr<DarkNewsDecay>, CrossSection, LI::crosssections::pyDarkNewsDecay> DarkNewsDecay(m, "DarkNewsDecay");

    DarkNewsDecay
        .def(init<>())
        .def("__eq__", [](const LI::crosssections::DarkNewsDecay &self, const LI::crosssections::DarkNewsDecay &other){ return self == other; })
        .def("equal", &LI::crosssections::DarkNewsDecay::equal)
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&DarkNewsDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&DarkNewsDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&DarkNewsDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&DarkNewsDecay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&DarkNewsDecay::DensityVariables)
        .def("FinalStateProbability",&DarkNewsDecay::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsDecay::SampleFinalState)
        .def("get_self", &pyDarkNewsDecay::get_self)
        .def(pybind11::pickle(
            [](const LI::crosssections::DarkNewsDecay & cpp_obj) {
                pybind11::object self;
                if(dynamic_cast<LI::crosssections::pyDarkNewsDecay const *>(&cpp_obj) != nullptr and dynamic_cast<LI::crosssections::pyDarkNewsDecay const *>(&cpp_obj)->self) {
                    self = pybind11::reinterpret_borrow<pybind11::object>(dynamic_cast<LI::crosssections::pyDarkNewsDecay const *>(&cpp_obj)->self);
                } else {
                    auto *tinfo = pybind11::detail::get_type_info(typeid(LI::crosssections::DarkNewsDecay));
                    pybind11::handle self_handle = get_object_handle(static_cast<const LI::crosssections::DarkNewsDecay *>(&cpp_obj), tinfo);
                    self = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
                }
                pybind11::dict d;
                if (pybind11::hasattr(self, "__dict__")) {
                    d = self.attr("__dict__");
                }
                return pybind11::make_tuple(d);
            },
            [](const pybind11::tuple &t) {
                if (t.size() != 1) {
                    throw std::runtime_error("Invalid state!");
                }
                auto cpp_state = std::unique_ptr<LI::crosssections::DarkNewsDecay>(new LI::crosssections::pyDarkNewsDecay);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;
}