#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#include "pyDarkNewsCrossSection.h"

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/DarkNewsCrossSection.h"

void register_DarkNewsCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<DarkNewsCrossSection, std::shared_ptr<DarkNewsCrossSection>, siren::interactions::pyDarkNewsCrossSection, CrossSection> DarkNewsCrossSection(m, "DarkNewsCrossSection");

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
        ;

    RegisterTrampolinePickleMethods(DarkNewsCrossSection,pyDarkNewsCrossSection);
}

