#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/CharmDISFromSpline.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_DMesonELoss(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<DMesonELoss, std::shared_ptr<DMesonELoss>, CrossSection> dmesoneloss(m, "DMesonELoss");

    dmesoneloss
        .def(init<>())
        .def(self == self)
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&DMesonELoss::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::Particle::ParticleType, double>(&DMesonELoss::TotalCrossSection, const_))
        // .def("DifferentialCrossSection",overload_cast<siren::dataclasses::CrossSectionDistributionRecord const &>(&DMesonELoss::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DMesonELoss::InteractionThreshold)
        .def("GetPossibleTargets",&DMesonELoss::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DMesonELoss::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DMesonELoss::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DMesonELoss::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DMesonELoss::GetPossibleSignaturesFromParents)
        .def("FinalStateProbability",&DMesonELoss::FinalStateProbability);
}

