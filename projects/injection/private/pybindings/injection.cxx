
#include <vector>

#include "../../public/LeptonInjector/injection/Process.h"
#include "../../public/LeptonInjector/injection/InjectorBase.h"
#include "../../public/LeptonInjector/injection/ColumnDepthLeptonInjector.h"
#include "../../public/LeptonInjector/injection/CylinderVolumeLeptonInjector.h"
#include "../../public/LeptonInjector/injection/DecayRangeLeptonInjector.h"
#include "../../public/LeptonInjector/injection/RangedLeptonInjector.h"
#include "../../public/LeptonInjector/injection/TreeWeighter.h"
#include "../../public/LeptonInjector/injection/Weighter.h"

#include "../../../utilities/public/LeptonInjector/utilities/Random.h"
#include "../../../detector/public/LeptonInjector/detector/EarthModel.h"
#include "../../../crosssections/public/LeptonInjector/crosssections/CrossSectionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(injection,m) {
  using namespace LI::injection;

  // Process

  class_<Process, std::shared_ptr<Process>>(m, "Process")
    .def_readwrite("primary_type",&Process::primary_type)
    .def_readwrite("cross_sections",&Process::cross_sections)
    .def("SetCrossSections",&Process::SetCrossSections);

  class_<InjectionProcess, std::shared_ptr<InjectionProcess>, Process>(m, "InjectionProcess")
    .def(init<>())
    .def("AddInjectionDistribution",&InjectionProcess::AddInjectionDistribution)
    .def_readwrite("injection_distributions",&InjectionProcess::injection_distributions);

  class_<PhysicalProcess, std::shared_ptr<PhysicalProcess>, Process>(m, "PhysicalProcess")
    .def(init<>())
    .def("AddPhysicalDistribution",&PhysicalProcess::AddPhysicalDistribution)
    .def_readwrite("physical_distributions",&PhysicalProcess::physical_distributions);

  // Injection

  class_<InjectorBase, std::shared_ptr<InjectorBase>>(m, "InjectorBase")
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::utilities::LI_random>>())
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::shared_ptr<LI::utilities::LI_random>>())
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>>())
    .def("SetStoppingCondition",&InjectorBase::SetStoppingCondition)
    .def("SetPrimaryProcess",&InjectorBase::SetPrimaryProcess)
    .def("AddSecondaryProcess",&InjectorBase::AddSecondaryProcess)
    .def("GetPrimaryProcess",&InjectorBase::GetPrimaryProcess)
    .def("GetSecondaryProcesses",&InjectorBase::GetSecondaryProcesses)
    .def("GetSecondaryProcessMap",&InjectorBase::GetSecondaryProcessMap)
    .def("NewRecord",&InjectorBase::NewRecord)
    .def("SetRandom",&InjectorBase::SetRandom)
    .def("GenerateEvent",&InjectorBase::GenerateEvent)
    .def("DensityVariables",&InjectorBase::DensityVariables)
    .def("Name",&InjectorBase::Name)
    .def("GetInjectionDistributions",&InjectorBase::GetInjectionDistributions)
    .def("GetEarthModel",&InjectorBase::GetEarthModel)
    .def("GetCrossSections",&InjectorBase::GetCrossSections)
    .def("InjectedEvents",&InjectorBase::InjectedEvents)
    .def("EventsToInject",&InjectorBase::EventsToInject);

  class_<RangedLeptonInjector, std::shared_ptr<RangedLeptonInjector>, InjectorBase>(m, "RangedLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::RangeFunction>, double, double>())
    .def("Name",&RangedLeptonInjector::Name);

  class_<DecayRangeLeptonInjector, std::shared_ptr<DecayRangeLeptonInjector>, InjectorBase>(m, "DecayRangeLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::DecayRangeFunction>, double, double>())
    .def("Name",&DecayRangeLeptonInjector::Name);

  class_<ColumnDepthLeptonInjector, std::shared_ptr<ColumnDepthLeptonInjector>, InjectorBase>(m, "ColumnDepthLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::DepthFunction>, double, double>())
    .def("Name",&ColumnDepthLeptonInjector::Name)
    .def("InjectionBounds",&ColumnDepthLeptonInjector::InjectionBounds);

  class_<CylinderVolumeLeptonInjector, std::shared_ptr<CylinderVolumeLeptonInjector>, InjectorBase>(m, "CylinderVolumeLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, LI::geometry::Cylinder>())
    .def("Name",&CylinderVolumeLeptonInjector::Name);

  // Weighter classes

  class_<LeptonProcessWeighter, std::shared_ptr<LeptonProcessWeighter>>(m, "LeptonProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<InjectionProcess>, std::shared_ptr<LI::detector::EarthModel>>())
    .def("InteractionProbability",&LeptonProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&LeptonProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",&LeptonProcessWeighter::PhysicalProbability)
    .def("GenerationProbability",&LeptonProcessWeighter::GenerationProbability)
    .def("EventWeight",&LeptonProcessWeighter::EventWeight);

  class_<LeptonTreeWeighter, std::shared_ptr<LeptonTreeWeighter>>(m, "LeptonTreeWeighter")
    .def(init<std::vector<std::shared_ptr<InjectorBase>>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<PhysicalProcess>, std::vector<std::shared_ptr<PhysicalProcess>>>())
    .def(init<std::vector<std::shared_ptr<InjectorBase>>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<PhysicalProcess>>())
    .def("EventWeight",&LeptonTreeWeighter::EventWeight);

  class_<LeptonWeighter, std::shared_ptr<LeptonWeighter>>(m, "LeptonWeighter")
    .def(init<std::vector<std::shared_ptr<InjectorBase>>, std::shared_ptr<LI::detector::EarthModel>, std::shared_ptr<LI::crosssections::CrossSectionCollection>, std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>>>())
    .def("EventWeight",&LeptonWeighter::EventWeight)
    .def("SimplifiedEventWeight",&LeptonWeighter::SimplifiedEventWeight);




}
