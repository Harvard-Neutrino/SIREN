#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/interactions/CrossSection.h"
#include "../../public/LeptonInjector/interactions/DipoleFromTable.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_DipoleFromTable(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::interactions;

    class_<DipoleFromTable, std::shared_ptr<DipoleFromTable>, CrossSection> dipolefromtable(m, "DipoleFromTable");

    dipolefromtable
        .def(init<double, double, DipoleFromTable::HelicityChannel>())
        .def(init<double, double, DipoleFromTable::HelicityChannel, bool, bool>())
        .def(init<double, double, DipoleFromTable::HelicityChannel, bool, bool, bool>())
        .def(init<double, double, DipoleFromTable::HelicityChannel, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def(self == self)
        .def(init<double, double, DipoleFromTable::HelicityChannel, bool, bool, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def(init<double, double, DipoleFromTable::HelicityChannel, bool, bool, bool, std::set<LI::dataclasses::Particle::ParticleType>>())
        .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DipoleFromTable::TotalCrossSection, const_))
        .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType>(&DipoleFromTable::TotalCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DipoleFromTable::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType, double, double>(&DipoleFromTable::DifferentialCrossSection, const_))
        .def("DifferentialCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType, double, double, double>(&DipoleFromTable::DifferentialCrossSection, const_))
        .def("InteractionThreshold",&DipoleFromTable::InteractionThreshold)
        .def("GetPossibleTargets",&DipoleFromTable::GetPossibleTargets)
        .def("GetPossibleTargetsFromPrimary",&DipoleFromTable::GetPossibleTargetsFromPrimary)
        .def("GetPossiblePrimaries",&DipoleFromTable::GetPossiblePrimaries)
        .def("GetPossibleSignatures",&DipoleFromTable::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParents",&DipoleFromTable::GetPossibleSignaturesFromParents)
        .def("AddDifferentialCrossSectionFile",&DipoleFromTable::AddDifferentialCrossSectionFile)
        .def("AddTotalCrossSectionFile",&DipoleFromTable::AddTotalCrossSectionFile)
        .def("AddDifferentialCrossSection",&DipoleFromTable::AddDifferentialCrossSection)
        .def("AddTotalCrossSection",&DipoleFromTable::AddTotalCrossSection)
        .def("DensityVariables",&DipoleFromTable::DensityVariables)
        .def("FinalStateProbability",&DipoleFromTable::FinalStateProbability)
        ;

    enum_<DipoleFromTable::HelicityChannel>(dipolefromtable, "HelicityChannel")
        .value("Conserving",DipoleFromTable::HelicityChannel::Conserving)
        .value("Flipping",DipoleFromTable::HelicityChannel::Flipping)
        .export_values()
        ;
}

