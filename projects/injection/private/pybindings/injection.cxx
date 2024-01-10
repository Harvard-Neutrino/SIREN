
#include <vector>

#include "../../public/LeptonInjector/injection/Process.h"
#include "../../public/LeptonInjector/injection/Injector.h"
#include "../../public/LeptonInjector/injection/ColumnDepthLeptonInjector.h"
#include "../../public/LeptonInjector/injection/CylinderVolumeLeptonInjector.h"
#include "../../public/LeptonInjector/injection/DecayRangeLeptonInjector.h"
#include "../../public/LeptonInjector/injection/RangedLeptonInjector.h"
#include "../../public/LeptonInjector/injection/TreeWeighter.h"
#include "../../public/LeptonInjector/injection/Weighter.h"
#include "../../public/LeptonInjector/injection/WeightingUtils.h"

#include "../../../distributions/public/LeptonInjector/distributions/primary/vertex/DepthFunction.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"
#include "../../../detector/public/LeptonInjector/detector/DetectorModel.h"
#include "../../../interactions/public/LeptonInjector/interactions/InteractionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(injection,m) {
  using namespace LI::injection;

  // Utils function
    
  m.def("CrossSectionProbability", &CrossSectionProbability);
    
  // Process

  class_<Process, std::shared_ptr<Process>>(m, "Process")
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("cross_sections", &Process::GetInteractions, &Process::SetInteractions);

  class_<InjectionProcess, std::shared_ptr<InjectionProcess>, Process>(m, "InjectionProcess")
    .def(init<>())
    .def("AddInjectionDistribution",&InjectionProcess::AddInjectionDistribution)
    .def("GetInjectionDistributions",&InjectionProcess::GetInjectionDistributions);

  class_<PhysicalProcess, std::shared_ptr<PhysicalProcess>, Process>(m, "PhysicalProcess")
    .def(init<>())
    .def("AddPhysicalDistribution",&PhysicalProcess::AddPhysicalDistribution)
    .def("GetPhysicalDistributions",&PhysicalProcess::GetPhysicalDistributions);

  // Injection

  class_<Injector, std::shared_ptr<Injector>>(m, "Injector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<LI::utilities::LI_random>>())
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::shared_ptr<LI::utilities::LI_random>>())
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>>())
    .def("SetStoppingCondition",&Injector::SetStoppingCondition)
    .def("SetPrimaryProcess",&Injector::SetPrimaryProcess)
    .def("AddSecondaryProcess",&Injector::AddSecondaryProcess)
    .def("GetPrimaryProcess",&Injector::GetPrimaryProcess)
    .def("GetSecondaryProcesses",&Injector::GetSecondaryProcesses)
    .def("GetSecondaryProcessMap",&Injector::GetSecondaryProcessMap)
    .def("NewRecord",&Injector::NewRecord)
    .def("SetRandom",&Injector::SetRandom)
    .def("GenerateEvent",&Injector::GenerateEvent)
    .def("DensityVariables",&Injector::DensityVariables)
    .def("Name",&Injector::Name)
    .def("GetInjectionDistributions",&Injector::GetInjectionDistributions)
    .def("GetDetectorModel",&Injector::GetDetectorModel)
    .def("GetInteractions",&Injector::GetInteractions)
    .def("InjectedEvents",&Injector::InjectedEvents)
    .def("EventsToInject",&Injector::EventsToInject);

  class_<RangedLeptonInjector, std::shared_ptr<RangedLeptonInjector>, Injector>(m, "RangedLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::RangeFunction>, double, double>())
    .def("Name",&RangedLeptonInjector::Name);

  class_<DecayRangeLeptonInjector, std::shared_ptr<DecayRangeLeptonInjector>, Injector>(m, "DecayRangeLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::DecayRangeFunction>, double, double>())
    .def("Name",&DecayRangeLeptonInjector::Name);

  class_<ColumnDepthLeptonInjector, std::shared_ptr<ColumnDepthLeptonInjector>, Injector>(m, "ColumnDepthLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::distributions::DepthFunction>, double, double>())
    .def("Name",&ColumnDepthLeptonInjector::Name)
    .def("InjectionBounds",&ColumnDepthLeptonInjector::InjectionBounds);

  class_<CylinderVolumeLeptonInjector, std::shared_ptr<CylinderVolumeLeptonInjector>, Injector>(m, "CylinderVolumeLeptonInjector")
    .def(init<unsigned int, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<InjectionProcess>, std::vector<std::shared_ptr<InjectionProcess>>, std::shared_ptr<LI::utilities::LI_random>, LI::geometry::Cylinder>())
    .def("Name",&CylinderVolumeLeptonInjector::Name);

  // Weighter classes

  class_<LeptonProcessWeighter, std::shared_ptr<LeptonProcessWeighter>>(m, "LeptonProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<InjectionProcess>, std::shared_ptr<LI::detector::DetectorModel>>())
    .def("InteractionProbability",&LeptonProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&LeptonProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",&LeptonProcessWeighter::PhysicalProbability)
    .def("GenerationProbability",&LeptonProcessWeighter::GenerationProbability)
    .def("EventWeight",&LeptonProcessWeighter::EventWeight);

  class_<LeptonTreeWeighter, std::shared_ptr<LeptonTreeWeighter>>(m, "LeptonTreeWeighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>, std::vector<std::shared_ptr<PhysicalProcess>>>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>>())
    .def("EventWeight",&LeptonTreeWeighter::EventWeight);

  class_<LeptonWeighter, std::shared_ptr<LeptonWeighter>>(m, "LeptonWeighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<LI::detector::DetectorModel>, std::shared_ptr<LI::interactions::InteractionCollection>, std::vector<std::shared_ptr<LI::distributions::WeightableDistribution>>>())
    .def("EventWeight",&LeptonWeighter::EventWeight)
    .def("SimplifiedEventWeight",&LeptonWeighter::SimplifiedEventWeight);




}
