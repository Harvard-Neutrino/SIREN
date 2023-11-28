#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/DarkNewsCrossSection.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

// Macro for defining pure virtual methods of PyDarkNewsCrossSection 
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

// Macro for defining virtual methods of PyDarkNewsCrossSection 
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
class pyCrossSection : public CrossSection {
public:
    using CrossSection::CrossSection;
    pyCrossSection(CrossSection && parent) : CrossSection(std::move(parent)) {}
    pybind11::object self;

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            interaction
        )
    }

    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            primary,
            energy,
            target
        )
    }

    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            interaction
        )
    }

    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            interaction
        )
    }

    void SampleFinalState(dataclasses::InteractionRecord & interaction, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            interaction,
            random
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            CrossSection,
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
// Trampoline class for DarkNewsCrossSection
class pyDarkNewsCrossSection : public DarkNewsCrossSection {
public:
    using DarkNewsCrossSection::DarkNewsCrossSection;
    pyDarkNewsCrossSection(DarkNewsCrossSection && parent) : DarkNewsCrossSection(std::move(parent)) {}
    pybind11::object self;

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            interaction
        )
    }

    double TotalCrossSection(LI::dataclasses::Particle::ParticleType primary, double energy, LI::dataclasses::Particle::ParticleType target) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            primary,
            energy,
            target
        )
    }

    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            interaction
        )
    }

    double DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary, LI::dataclasses::Particle::ParticleType target, double energy, double Q2) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            primary,
            target,
            energy,
            Q2
        )
    }

    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            interaction
        )
    }

    double Q2Min(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Min,
            "Q2Min",
            interaction
        )
    }

    double Q2Max(dataclasses::InteractionRecord const & interaction) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Max,
            "Q2Max",
            interaction
        )
    }

    void SampleFinalState(dataclasses::InteractionRecord & interaction, std::shared_ptr<LI::utilities::LI_random> random) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            interaction,
            random
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<LI::dataclasses::Particle::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override {
        C_PYBIND11_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<LI::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        C_PYBIND11_OVERRIDE(
            self,
            DarkNewsCrossSection,
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

void register_DarkNewsCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    // Bindings for pyDarkNewsCrossSection
    class_<LI::crosssections::pyDarkNewsCrossSection> pyDarkNewsCrossSection(m, "pyDarkNewsCrossSection");

    pyDarkNewsCrossSection
        .def(init<>())
        .def("__eq__", [](const LI::crosssections::DarkNewsCrossSection &self, const LI::crosssections::DarkNewsCrossSection &other){ return self == other; })
        .def("equal", &LI::crosssections::DarkNewsCrossSection::equal)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, LI::dataclasses::Particle::ParticleType, double, double>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DarkNewsCrossSection::InteractionThreshold)
        .def("Q2Min",&DarkNewsCrossSection::Q2Min)
        .def("Q2Max",&DarkNewsCrossSection::Q2Max)
        .def("GetPossibleTargets",&DarkNewsCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DarkNewsCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DarkNewsCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DarkNewsCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DarkNewsCrossSection::GetPossibleSignaturesFromParents)
        .def("DensityVariables",&DarkNewsCrossSection::DensityVariables)
        .def("FinalStateProbability",&DarkNewsCrossSection::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsCrossSection::SampleFinalState)
        .def("get_self", &pyDarkNewsCrossSection::get_self)
        .def(pybind11::pickle(
            [](const LI::crosssections::pyDarkNewsCrossSection & cpp_obj) {
                pybind11::object self;
                if(cpp_obj.self) {
                    self = pybind11::reinterpret_borrow<pybind11::object>(cpp_obj.self);
                } else {
                    auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsCrossSection));
                    pybind11::handle self_handle = get_object_handle(static_cast<const DarkNewsCrossSection *>(&cpp_obj), tinfo);
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
                auto cpp_state = std::unique_ptr<LI::crosssections::pyDarkNewsCrossSection>(new LI::crosssections::pyDarkNewsCrossSection);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;


    class_<DarkNewsCrossSection, std::shared_ptr<DarkNewsCrossSection>, CrossSection, LI::crosssections::pyDarkNewsCrossSection> DarkNewsCrossSection(m, "DarkNewsCrossSection");

    DarkNewsCrossSection
        .def(init<>())
        .def("__eq__", [](const LI::crosssections::DarkNewsCrossSection &self, const LI::crosssections::DarkNewsCrossSection &other){ return self == other; })
        .def("equal", &LI::crosssections::DarkNewsCrossSection::equal)
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, LI::dataclasses::Particle::ParticleType, double, double>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DarkNewsCrossSection::InteractionThreshold)
        .def("Q2Min",&DarkNewsCrossSection::Q2Min)
        .def("Q2Max",&DarkNewsCrossSection::Q2Max)
        .def("GetPossibleTargets",&DarkNewsCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DarkNewsCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DarkNewsCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DarkNewsCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DarkNewsCrossSection::GetPossibleSignaturesFromParents)
        .def("DensityVariables",&DarkNewsCrossSection::DensityVariables)
        .def("FinalStateProbability",&DarkNewsCrossSection::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsCrossSection::SampleFinalState)
        .def("get_self", &DarkNewsCrossSection::get_self)
        .def(pybind11::pickle(
            [](const LI::crosssections::DarkNewsCrossSection & cpp_obj) {
                pybind11::object self;
                if(dynamic_cast<LI::crosssections::pyDarkNewsCrossSection const *>(&cpp_obj) != nullptr and dynamic_cast<LI::crosssections::pyDarkNewsCrossSection const *>(&cpp_obj)->self) {
                    self = pybind11::reinterpret_borrow<pybind11::object>(dynamic_cast<LI::crosssections::pyDarkNewsCrossSection const *>(&cpp_obj)->self);
                } else {
                    auto *tinfo = pybind11::detail::get_type_info(typeid(LI::crosssections::DarkNewsCrossSection));
                    pybind11::handle self_handle = get_object_handle(static_cast<const LI::crosssections::DarkNewsCrossSection *>(&cpp_obj), tinfo);
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
                auto cpp_state = std::unique_ptr<LI::crosssections::DarkNewsCrossSection>(new LI::crosssections::pyDarkNewsCrossSection);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;
}

