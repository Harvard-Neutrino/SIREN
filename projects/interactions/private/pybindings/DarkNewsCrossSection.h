#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/DarkNewsCrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../utilities/public/SIREN/utilities/Pybind11Trampoline.h"

namespace siren {
namespace interactions {
// Trampoline class for CrossSection
class pyCrossSection : public CrossSection {
public:
    using CrossSection::CrossSection;
    pyCrossSection(CrossSection && parent) : CrossSection(std::move(parent)) {}
    pybind11::object self;

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            interaction
        )
    }

    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            CrossSection,
            double,
            TotalCrossSectionAllFinalStates,
            "TotalCrossSectionAllFinalStates",
            record
        )
    }

    double DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            interaction
        )
    }

    double InteractionThreshold(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            interaction
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::LI_random> random) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            record,
            random
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE_PURE(
            self,
            CrossSection,
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
// Trampoline class for DarkNewsCrossSection
class pyDarkNewsCrossSection : public DarkNewsCrossSection {
public:
    using DarkNewsCrossSection::DarkNewsCrossSection;
    pyDarkNewsCrossSection(DarkNewsCrossSection && parent) : DarkNewsCrossSection(std::move(parent)) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
    }
    pyDarkNewsCrossSection(DarkNewsCrossSection const & parent) : DarkNewsCrossSection(parent) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsCrossSection)))));
    }
    pybind11::object self;

    double TotalCrossSectionAllFinalStates(siren::dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            CrossSection,
            double,
            TotalCrossSectionAllFinalStates,
            "TotalCrossSectionAllFinalStates",
            std::cref(record)
        )
    }

    double TotalCrossSection(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TotalCrossSection,
            "TotalCrossSection",
            std::cref(interaction)
        )
    }

    double TotalCrossSection(siren::dataclasses::ParticleType primary, double energy, siren::dataclasses::ParticleType target) const override {
        SELF_OVERRIDE_PURE(
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
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            DifferentialCrossSection,
            "DifferentialCrossSection",
            std::cref(interaction)
        )
    }

    double DifferentialCrossSection(siren::dataclasses::ParticleType primary, siren::dataclasses::ParticleType target, double energy, double Q2) const override {
        SELF_OVERRIDE(
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
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            InteractionThreshold,
            "InteractionThreshold",
            std::cref(interaction)
        )
    }

    double Q2Min(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Min,
            "Q2Min",
            std::cref(interaction)
        )
    }

    double Q2Max(dataclasses::InteractionRecord const & interaction) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            Q2Max,
            "Q2Max",
            std::cref(interaction)
        )
    }

    double TargetMass(dataclasses::ParticleType const & target_type) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            TargetMass,
            "TargetMass",
            std::cref(target_type)
        )
    }

    std::vector<double> SecondaryMasses(std::vector<dataclasses::ParticleType> const & secondary_types) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryMasses,
            "SecondaryMasses",
            std::cref(secondary_types)
        )
    }

    std::vector<double> SecondaryHelicities(dataclasses::InteractionRecord const & record) const override{
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            std::vector<double>,
            SecondaryHelicities,
            "SecondaryHelicities",
            std::cref(record)
        )
    }

    void SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::LI_random> random) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargets() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargets,
            "GetPossibleTargets"
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(siren::dataclasses::ParticleType primary_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossibleTargetsFromPrimary,
            "GetPossibleTargetsFromPrimary",
            primary_type
        )
    }

    std::vector<siren::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::ParticleType>,
            GetPossiblePrimaries,
            "GetPossiblePrimaries"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(siren::dataclasses::ParticleType primary_type, siren::dataclasses::ParticleType target_type) const override {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsCrossSection,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParents,
            "GetPossibleSignaturesFromParents",
            primary_type,
            target_type
        )
    }

    double FinalStateProbability(dataclasses::InteractionRecord const & record) const override {
        SELF_OVERRIDE(
            self,
            DarkNewsCrossSection,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
        )
    }

    pybind11::object get_representation() override {
        const DarkNewsCrossSection * ref;
        if(self) {
            ref = self.cast<DarkNewsCrossSection *>();
        } else {
            ref = this;
        }
        auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsCrossSection));
        pybind11::function override_func =
            tinfo ? pybind11::detail::get_type_override(static_cast<const DarkNewsCrossSection *>(ref), tinfo, "get_representation") : pybind11::function();
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
            auto *tinfo = pybind11::detail::get_type_info(typeid(DarkNewsCrossSection));
            pybind11::handle self_handle = get_object_handle(static_cast<const DarkNewsCrossSection *>(this), tinfo);
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

void register_DarkNewsCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<DarkNewsCrossSection, std::shared_ptr<DarkNewsCrossSection>, CrossSection, siren::interactions::pyDarkNewsCrossSection> DarkNewsCrossSection(m, "DarkNewsCrossSection");

    DarkNewsCrossSection
        .def(init<>())
        .def("__eq__", [](const siren::interactions::DarkNewsCrossSection &self, const siren::interactions::DarkNewsCrossSection &other){ return self == other; })
        .def_readwrite("m_ups",&DarkNewsCrossSection::m_ups)
        .def_readwrite("m_target",&DarkNewsCrossSection::m_target)
        .def("equal", &siren::interactions::DarkNewsCrossSection::equal)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType>(&DarkNewsCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::ParticleType, siren::dataclasses::ParticleType, double, double>(&DarkNewsCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DarkNewsCrossSection::InteractionThreshold)
        .def("Q2Min",&DarkNewsCrossSection::Q2Min)
        .def("Q2Max",&DarkNewsCrossSection::Q2Max)
        .def("TargetMass",&DarkNewsCrossSection::TargetMass)
        .def("SecondaryMasses",&DarkNewsCrossSection::SecondaryMasses)
        .def("SecondaryHelicities",&DarkNewsCrossSection::SecondaryHelicities)
        .def("GetPossibleTargets",&DarkNewsCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DarkNewsCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DarkNewsCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DarkNewsCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DarkNewsCrossSection::GetPossibleSignaturesFromParents)
        .def("DensityVariables",&DarkNewsCrossSection::DensityVariables)
        .def("FinalStateProbability",&DarkNewsCrossSection::FinalStateProbability)
        .def("SampleFinalState",&DarkNewsCrossSection::SampleFinalState)
        .def("get_representation", &DarkNewsCrossSection::get_representation)
        .def(pybind11::pickle(
            [](siren::interactions::DarkNewsCrossSection & cpp_obj) {
                return pybind11::make_tuple(cpp_obj.get_representation());
            },
            [](const pybind11::tuple &t) {
                if (t.size() != 1) {
                    throw std::runtime_error("Invalid state!");
                }
                auto cpp_state = std::unique_ptr<siren::interactions::DarkNewsCrossSection>(new siren::interactions::pyDarkNewsCrossSection);
                auto py_state = t[0].cast<pybind11::dict>();
                return std::make_pair(std::move(cpp_state), py_state);
            })
        )
        ;
}

