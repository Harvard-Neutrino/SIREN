#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/TrivialCrossSection.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_TrivialCrossSection(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<TrivialCrossSection, std::shared_ptr<TrivialCrossSection>, CrossSection> trivialcrosssection(m, "TrivialCrossSection");

    trivialcrosssection
        .def(init<double, std::vector<siren::dataclasses::ParticleType>, std::vector<siren::dataclasses::ParticleType>>(),
                arg("cross_section_cm2"), arg("primary_types"), arg("target_types"))
        .def(init<std::map<siren::dataclasses::ParticleType, std::pair<std::vector<double>, std::vector<double>>>, std::vector<siren::dataclasses::ParticleType>>(),
                arg("tables"), arg("target_types"))
        .def("_equal", &TrivialCrossSection::equal)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&TrivialCrossSection::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double>(&TrivialCrossSection::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&TrivialCrossSection::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&TrivialCrossSection::InteractionThreshold)
        .def("GetPossibleTargets",&TrivialCrossSection::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&TrivialCrossSection::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&TrivialCrossSection::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&TrivialCrossSection::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&TrivialCrossSection::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&TrivialCrossSection::FinalStateProbability)
        .def("DensityVariables",&TrivialCrossSection::DensityVariables)
        ;
}
