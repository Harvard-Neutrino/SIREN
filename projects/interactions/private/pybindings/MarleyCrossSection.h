#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/MarleyCrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_MarleyCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<MarleyCrossSection, std::shared_ptr<MarleyCrossSection>, CrossSection> marleycrosssection(m, "MarleyCrossSection");

    marleycrosssection

        .def(init<std::string, std::string, std::string, std::string, std::string>())
        .def(self == self)
        .def("TotalCrossSection",&MarleyCrossSection::TotalCrossSection)
        .def("TotalCrossSectionAllFinalStates",&MarleyCrossSection::TotalCrossSectionAllFinalStates)
        .def("DifferentialCrossSection",&MarleyCrossSection::DifferentialCrossSection)
        .def("InteractionThreshold",&MarleyCrossSection::InteractionThreshold)
        .def("GetPossibleTargets",&MarleyCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&MarleyCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&MarleyCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&MarleyCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&MarleyCrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&MarleyCrossSection::FinalStateProbability)
        .def("DensityVariables",&MarleyCrossSection::DensityVariables);
}
