#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "../../public/LeptonInjector/interactions/DarkNewsDecay.h"
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
namespace interactions {
// Trampoline class for Decay
class pyDecay : public Decay {
public:
    using Decay::Decay;
    pyDecay(Decay && parent) : Decay(std::move(parent)) {}
    pybind11::object self;

    double TotalDecayLength(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayLength,
            "TotalDecayLength",
            interaction
        )
    }

    double TotalDecayLengthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            Decay,
            double,
            TotalDecayLengthForFinalState,
            "TotalDecayLengthForFinalState",
            interaction
        )
    }
    
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

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            void,
            SampleFinalState,
            "SampleFinalState",
            record,
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

    std::vector<std::string> DensityVariables() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            Decay,
            std::vector<std::string>,
            DensityVariables,
            "DensityVariables"
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

    pybind11::object get_representation() {
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
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            std::cref(interaction)
        )
    }

    double TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            std::cref(interaction)
        )
    }

    double TotalDecayWidth(LI::dataclasses::Particle::ParticleType primary) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

    double DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            std::cref(interaction)
        )
    }

    void SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleRecordFromDarkNews,
            "SampleRecordFromDarkNews",
            std::ref(record),
            random
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParent(LI::dataclasses::Particle::ParticleType primary_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParent,
            "GetPossibleSignaturesFromParent",
            primary_type
        )
    }

    std::vector<std::string> DensityVariables() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<std::string>,
            DensityVariables,
            "DensityVariables"
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
        )
    }

    pybind11::object get_representation() override {
        const DarkNewsDecay * ref;
        if(self) {
            ref = self.cast<DarkNewsDecay *>();
        } else {
            ref = this;
        }
        auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsDecay));
        pybind11::function override_func =
            tinfo ? pybind11::detail::get_type_override(static_cast<const DarkNewsDecay *>(ref), tinfo, "get_representation") : pybind11::function();
        if (override_func) {
            pybind11::object o = override_func();
            if(not pybind11::isinstance<pybind11::dict>(o)) {
                throw std::runtime_error("get_representation must return a dict");
            }
            return o;
        }

        pybind11::object _self;
        if(this->self) {
            self = pybind11::reinterpret_borrow<pybind11::object>(this->self);
        } else {
            auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsDecay));
            pybind11::handle self_handle = get_object_handle(static_cast<const DarkNewsDecay *>(this), tinfo);
            _self = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
        }
        pybind11::dict d;
        if (pybind11::hasattr(self, "__dict__")) {
            d = _self.attr("__dict__");
        }
        return d;
    }
};
} // end interactions namespace
} // end LI namespace

void register_DarkNewsDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    // Bindings for pyDarkNewsDecay
    class_<LI::interactions::pyDarkNewsDecay> pyDarkNewsDecay(m, "pyDarkNewsDecay");

    pyDarkNewsDecay
        .def(init<>())
        .def("__eq__", [](const LI::interactions::DarkNewsDecay &self, const LI::interactions::DarkNewsDecay &other){ return self == other; })
        .def("equal", &LI::interactions::DarkNewsDecay::equal)
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&DarkNewsDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&DarkNewsDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&DarkNewsDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&DarkNewsDecay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&DarkNewsDecay::DensityVariables)
        .def("FinalStateProbability",&DarkNewsDecay::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsDecay::SampleFinalState)
        .def("SampleRecordFromDarkNews",&DarkNewsDecay::SampleRecordFromDarkNews)
        .def("get_representation", &pyDarkNewsDecay::get_representation)
        .def(pybind11::pickle(
            [](LI::interactions::pyDarkNewsDecay & cpp_obj) {
                return pybind11::make_tuple(cpp_obj.get_representation());
            },
            [](const pybind11::tuple &t) {
                if (t.size() != 1) {
                    throw std::runtime_error("Invalid state!");
                }
                auto cpp_state = std::unique_ptr<LI::interactions::pyDarkNewsDecay>(new LI::interactions::pyDarkNewsDecay);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;


    class_<DarkNewsDecay, std::shared_ptr<DarkNewsDecay>, Decay, LI::interactions::pyDarkNewsDecay> DarkNewsDecay(m, "DarkNewsDecay");

    DarkNewsDecay
        .def(init<>())
        .def("__eq__", [](const LI::interactions::DarkNewsDecay &self, const LI::interactions::DarkNewsDecay &other){ return self == other; })
        .def("equal", &LI::interactions::DarkNewsDecay::equal)
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&DarkNewsDecay::TotalDecayWidth, const_))
        .def("TotalDecayWidthForFinalState",&DarkNewsDecay::TotalDecayWidthForFinalState)
        .def("DifferentialDecayWidth",&DarkNewsDecay::DifferentialDecayWidth)
        .def("GetPossibleSignatures",&DarkNewsDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&DarkNewsDecay::GetPossibleSignaturesFromParent)
        .def("DensityVariables",&DarkNewsDecay::DensityVariables)
        .def("FinalStateProbability",&DarkNewsDecay::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsDecay::SampleFinalState)
        .def("SampleRecordFromDarkNews",&DarkNewsDecay::SampleRecordFromDarkNews)
        .def("get_representation", &DarkNewsDecay::get_representation)
        .def(pybind11::pickle(
            [](LI::interactions::DarkNewsDecay & cpp_obj) {
                return pybind11::make_tuple(cpp_obj.get_representation());
            },
            [](const pybind11::tuple &t) {
                if (t.size() != 1) {
                    throw std::runtime_error("Invalid state!");
                }
                auto cpp_state = std::unique_ptr<LI::interactions::DarkNewsDecay>(new LI::interactions::pyDarkNewsDecay);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;
}
