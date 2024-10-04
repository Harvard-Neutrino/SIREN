#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/HNLDipoleFromTable.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_HNLDipoleFromTable(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<HNLDipoleFromTable, std::shared_ptr<HNLDipoleFromTable>, CrossSection> HNLDipoleFromTable(m, "HNLDipoleFromTable");

    HNLDipoleFromTable
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel>())
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel, bool, bool>())
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel, bool, bool, bool>())
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel, std::set<siren::dataclasses::ParticleType>>())
        .def(self == self)
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel, bool, bool, std::set<siren::dataclasses::ParticleType>>())
        .def(init<double, double, HNLDipoleFromTable::HelicityChannel, bool, bool, bool, std::set<siren::dataclasses::ParticleType>>())
        .def("TotalCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDipoleFromTable::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType>(&HNLDipoleFromTable::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::InteractionRecord const &>(&HNLDipoleFromTable::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType, double, double>(&HNLDipoleFromTable::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<siren::dataclasses::ParticleType, double, siren::dataclasses::ParticleType, double, double, double>(&HNLDipoleFromTable::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&HNLDipoleFromTable::InteractionThreshold)
        .def("GetPossibleTargets",&HNLDipoleFromTable::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&HNLDipoleFromTable::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&HNLDipoleFromTable::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&HNLDipoleFromTable::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&HNLDipoleFromTable::GetPossibleSignaturesFromParents)
        .def("AddDifferentialCrossSectionFile",&HNLDipoleFromTable::AddDifferentialCrossSectionFile)
        .def("AddTotalCrossSectionFile",&HNLDipoleFromTable::AddTotalCrossSectionFile)
        .def("AddDifferentialCrossSection",&HNLDipoleFromTable::AddDifferentialCrossSection)
        .def("AddTotalCrossSection",&HNLDipoleFromTable::AddTotalCrossSection)
        .def("DensityVariables",&HNLDipoleFromTable::DensityVariables)
        .def("FinalStateProbability",&HNLDipoleFromTable::FinalStateProbability)
        ;

    enum_<HNLDipoleFromTable::HelicityChannel>(HNLDipoleFromTable, "HelicityChannel")
        .value("Conserving",HNLDipoleFromTable::HelicityChannel::Conserving)
        .value("Flipping",HNLDipoleFromTable::HelicityChannel::Flipping)
        .export_values()
        ;
}

