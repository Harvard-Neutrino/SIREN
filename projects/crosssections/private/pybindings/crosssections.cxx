
#include <vector>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/Decay.h"
#include "../../public/LeptonInjector/crosssections/DipoleFromTable.h"
#include "../../public/LeptonInjector/crosssections/NeutrissimoDecay.h"
#include "../../public/LeptonInjector/crosssections/CrossSectionCollection.h"
#include "../../public/LeptonInjector/crosssections/DISFromSpline.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(crosssections,m) {
  using namespace LI::crosssections;

  class_<CrossSection, std::shared_ptr<CrossSection>>(m, "CrossSection");

  class_<Decay, std::shared_ptr<Decay>>(m, "Decay")
    .def("TotalDecayLength",&Decay::TotalDecayLength)
    .def("TotalDecayLengthForFinalState",&Decay::TotalDecayLengthForFinalState);

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
    .def("FinalStateProbability",&DipoleFromTable::FinalStateProbability);

  enum_<DipoleFromTable::HelicityChannel>(dipolefromtable, "HelicityChannel")
    .value("Conserving",DipoleFromTable::HelicityChannel::Conserving)
    .value("Flipping",DipoleFromTable::HelicityChannel::Flipping)
    .export_values();
    
    
  class_<DISFromSpline, std::shared_ptr<DISFromSpline>, CrossSection> disfromspline(m, "DISFromSpline");
   
  disfromspline
    .def(init<>())
    .def(init<std::vector<char>, std::vector<char>, int, double, double, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def(init<std::vector<char>, std::vector<char>, int, double, double, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
    .def(init<std::string, std::string, int, double, double, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def(init<std::string, std::string, std::set<LI::dataclasses::Particle::ParticleType>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def(init<std::string, std::string, int, double, double, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
    .def(init<std::string, std::string, std::vector<LI::dataclasses::Particle::ParticleType>, std::vector<LI::dataclasses::Particle::ParticleType>>())
    .def(self == self)
    .def("TotalCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DISFromSpline::TotalCrossSection, const_))
    .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double>(&DISFromSpline::TotalCrossSection, const_))
    .def("TotalCrossSection",overload_cast<LI::dataclasses::Particle::ParticleType, double, LI::dataclasses::Particle::ParticleType>(&DISFromSpline::TotalCrossSection, const_))
    .def("DifferentialCrossSection",overload_cast<LI::dataclasses::InteractionRecord const &>(&DISFromSpline::DifferentialCrossSection, const_))
    .def("DifferentialCrossSection",overload_cast<double, double, double, double>(&DISFromSpline::DifferentialCrossSection, const_))
    .def("InteractionThreshold",&DISFromSpline::InteractionThreshold)
    .def("SampleFinalState",&DISFromSpline::SampleFinalState)
    .def("GetPossibleTargets",&DISFromSpline::GetPossibleTargets)
    .def("GetPossibleTargetsFromPrimary",&DISFromSpline::GetPossibleTargetsFromPrimary)
    .def("GetPossiblePrimaries",&DISFromSpline::GetPossiblePrimaries)
    .def("GetPossibleSignatures",&DISFromSpline::GetPossibleSignatures)
    .def("GetPossibleSignaturesFromParents",&DISFromSpline::GetPossibleSignaturesFromParents)
    .def("FinalStateProbability",&DISFromSpline::FinalStateProbability);
      

  class_<NeutrissimoDecay, std::shared_ptr<NeutrissimoDecay>, Decay> neutrissimodecay(m, "NeutrissimoDecay");

  neutrissimodecay
    .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature>())
    .def(init<double, std::vector<double>, NeutrissimoDecay::ChiralNature, std::set<LI::dataclasses::Particle::ParticleType> const &>())
    .def(init<double, double, NeutrissimoDecay::ChiralNature>())
    .def(init<double, double, NeutrissimoDecay::ChiralNature, std::set<LI::dataclasses::Particle::ParticleType> const &>())
    .def(self == self)
    .def("GetHNLMass",&NeutrissimoDecay::GetHNLMass)
    .def("TotalDecayWidth",overload_cast<LI::dataclasses::InteractionRecord const &>(&NeutrissimoDecay::TotalDecayWidth, const_))
    .def("TotalDecayWidth",overload_cast<LI::dataclasses::Particle::ParticleType>(&NeutrissimoDecay::TotalDecayWidth, const_))
    .def("TotalDecayWidthForFinalState",&NeutrissimoDecay::TotalDecayWidthForFinalState)
    .def("DifferentialDecayWidth",&NeutrissimoDecay::DifferentialDecayWidth)
    .def("GetPossibleSignatures",&NeutrissimoDecay::GetPossibleSignatures)
    .def("GetPossibleSignaturesFromParent",&NeutrissimoDecay::GetPossibleSignaturesFromParent)
    .def("FinalStateProbability",&NeutrissimoDecay::FinalStateProbability);

  enum_<NeutrissimoDecay::ChiralNature>(neutrissimodecay, "ChiralNature")
    .value("Dirac",NeutrissimoDecay::ChiralNature::Dirac)
    .value("Majorana",NeutrissimoDecay::ChiralNature::Majorana)
    .export_values();

  class_<CrossSectionCollection, std::shared_ptr<CrossSectionCollection>>(m, "CrossSectionCollection")
    .def(init<>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<Decay>>>())
    .def(init<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>, std::vector<std::shared_ptr<Decay>>>())
    .def(self == self)
    .def("GetDecays",&CrossSectionCollection::GetDecays)
    .def("HasCrossSections",&CrossSectionCollection::HasCrossSections)
    .def("HasDecays",&CrossSectionCollection::HasDecays)
    .def("GetCrossSectionsForTarget",&CrossSectionCollection::GetCrossSectionsForTarget)
    .def("GetCrossSectionsByTarget",&CrossSectionCollection::GetCrossSectionsByTarget)
    .def("TargetTypes",&CrossSectionCollection::TargetTypes)
    .def("TotalDecayWidth",&CrossSectionCollection::TotalDecayWidth)
    .def("TotalDecayLength",&CrossSectionCollection::TotalDecayLength)
    .def("MatchesPrimary",&CrossSectionCollection::MatchesPrimary);

}
