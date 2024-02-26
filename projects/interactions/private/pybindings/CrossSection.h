#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

using namespace pybind11;
using namespace LI::interactions;
class PyCrossSection : public LI::interactions::CrossSection {
public:
    using CrossSection::CrossSection;

    bool equal(const CrossSection& other) const override {
        PYBIND11_OVERRIDE_PURE_NAME(
            bool,
            CrossSection,
            "_equal",
            equal,
            other
        );
    }

    double TotalCrossSection(LI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            TotalCrossSection,
            record
        );
    }

    double TotalCrossSectionAllFinalStates(LI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE(
            double,
            CrossSection,
            TotalCrossSectionAllFinalStates,
            record
        );
    }

    double DifferentialCrossSection(LI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            DifferentialCrossSection,
            record
        );
    }

    double InteractionThreshold(LI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            InteractionRecord,
            record
        );
    }

    void SampleFinalState(LI::dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<LI::utilities::LI_random> rand) const override {
        PYBIND11_OVERRIDE_PURE(
            void,
            CrossSection,
            SampleFinalState,
            record,
            rand
        );
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargets() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<LI::dataclasses::Particle::ParticleType>,
            CrossSection,
            GetPossibleTargets
        );
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<LI::dataclasses::Particle::ParticleType>,
            CrossSection,
            GetPossibleTargetsFromPrimary,
            primary_type
        );
    }

    std::vector<LI::dataclasses::Particle::ParticleType> GetPossiblePrimaries() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<LI::dataclasses::Particle::ParticleType>,
            CrossSection,
            GetPossiblePrimaries
        );
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<LI::dataclasses::InteractionSignature>,
            CrossSection,
            GetPossibleSignatures
        );
    }

    std::vector<LI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<LI::dataclasses::InteractionSignature>,
            CrossSection,
            GetPossibleSignaturesFromParents,
            primary_type,
            target_type
        );
    }

    double FinalStateProbability(LI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            FinalStateProbability,
            record
        );
    }

    std::vector<std::string> DensityVariables() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<std::string>,
            CrossSection,
            DensityVariables
        );
    }
};

void register_CrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<CrossSection, std::shared_ptr<CrossSection>, PyCrossSection>(m, "CrossSection")
        .def(init<>())
        .def("__eq__", [](const CrossSection &self, const CrossSection &other){ return self == other; })
        .def("equal", &CrossSection::equal)
        .def("TotalCrossSection", (double (CrossSection::*)(LI::dataclasses::InteractionRecord const &) const)(&CrossSection::TotalCrossSection))
        .def("TotalCrossSection", (double (CrossSection::*)(LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType) const)(&CrossSection::TotalCrossSection))
        .def("DifferentialCrossSection", &CrossSection::DifferentialCrossSection)
        .def("InteractionThreshold", &CrossSection::InteractionThreshold)
        .def("SampleFinalState", (void (CrossSection::*)(LI::dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<LI::utilities::LI_random>) const)(&CrossSection::SampleFinalState))
        .def("GetPossibleTargets", &CrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary", &CrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries", &CrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures", &CrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents", &CrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability", &CrossSection::FinalStateProbability)
        .def("DensityVariables", &CrossSection::DensityVariables)
        ;
}

