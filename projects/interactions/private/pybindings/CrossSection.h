#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

using namespace pybind11;
using namespace SI::interactions;
class PyCrossSection : public SI::interactions::CrossSection {
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

    double TotalCrossSection(SI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            TotalCrossSection,
            record
        );
    }

    double TotalCrossSectionAllFinalStates(SI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_NAME(
            double, // Return type (ret_type)
            CrossSection,      // Parent class (cname)
            "TotalCrossSectionAllFinalStates",   // Name of method in Python (name)
            TotalCrossSectionAllFinalStates,    // Name of function in C++ (fn)
            record
        );
    }

    double DifferentialCrossSection(SI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            DifferentialCrossSection,
            record
        );
    }

    double InteractionThreshold(SI::dataclasses::InteractionRecord const & record) const override {
        PYBIND11_OVERRIDE_PURE(
            double,
            CrossSection,
            InteractionRecord,
            record
        );
    }

    void SampleFinalState(SI::dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<SI::utilities::LI_random> rand) const override {
        PYBIND11_OVERRIDE_PURE(
            void,
            CrossSection,
            SampleFinalState,
            record,
            rand
        );
    }

    std::vector<SI::dataclasses::ParticleType> GetPossibleTargets() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SI::dataclasses::ParticleType>,
            CrossSection,
            GetPossibleTargets
        );
    }

    std::vector<SI::dataclasses::ParticleType> GetPossibleTargetsFromPrimary(SI::dataclasses::ParticleType primary_type) const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SI::dataclasses::ParticleType>,
            CrossSection,
            GetPossibleTargetsFromPrimary,
            primary_type
        );
    }

    std::vector<SI::dataclasses::ParticleType> GetPossiblePrimaries() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SI::dataclasses::ParticleType>,
            CrossSection,
            GetPossiblePrimaries
        );
    }

    std::vector<SI::dataclasses::InteractionSignature> GetPossibleSignatures() const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SI::dataclasses::InteractionSignature>,
            CrossSection,
            GetPossibleSignatures
        );
    }

    std::vector<SI::dataclasses::InteractionSignature> GetPossibleSignaturesFromParents(SI::dataclasses::ParticleType primary_type, SI::dataclasses::ParticleType target_type) const override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<SI::dataclasses::InteractionSignature>,
            CrossSection,
            GetPossibleSignaturesFromParents,
            primary_type,
            target_type
        );
    }

    double FinalStateProbability(SI::dataclasses::InteractionRecord const & record) const override {
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
    using namespace SI::interactions;

    class_<CrossSection, std::shared_ptr<CrossSection>, PyCrossSection>(m, "CrossSection")
        .def(init<>())
        .def("__eq__", [](const CrossSection &self, const CrossSection &other){ return self == other; })
        .def("equal", &CrossSection::equal)
        .def("TotalCrossSection", (double (CrossSection::*)(SI::dataclasses::InteractionRecord const &) const)(&CrossSection::TotalCrossSection))
        .def("TotalCrossSectionAllFinalStates", (double (CrossSection::*)(SI::dataclasses::InteractionRecord const &) const)(&CrossSection::TotalCrossSectionAllFinalStates))
        .def("TotalCrossSection", (double (CrossSection::*)(SI::dataclasses::ParticleType, double, SI::dataclasses::ParticleType) const)(&CrossSection::TotalCrossSection))
        .def("DifferentialCrossSection", &CrossSection::DifferentialCrossSection)
        .def("InteractionThreshold", &CrossSection::InteractionThreshold)
        .def("SampleFinalState", (void (CrossSection::*)(SI::dataclasses::CrossSectionDistributionRecord &, std::shared_ptr<SI::utilities::LI_random>) const)(&CrossSection::SampleFinalState))
        .def("GetPossibleTargets", &CrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary", &CrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries", &CrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures", &CrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents", &CrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability", &CrossSection::FinalStateProbability)
        .def("DensityVariables", &CrossSection::DensityVariables)
        ;
}

