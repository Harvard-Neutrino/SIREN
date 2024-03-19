
#include <vector>

#include "../../public/SIREN/injection/Process.h"
#include "../../public/SIREN/injection/Injector.h"
//#include "../../public/SIREN/injection/ColumnDepthSIREN.h"
//#include "../../public/SIREN/injection/CylinderVolumeSIREN.h"
//#include "../../public/SIREN/injection/DecayRangeSIREN.h"
//#include "../../public/SIREN/injection/RangedSIREN.h"
#include "../../public/SIREN/injection/TreeWeighter.h"
#include "../../public/SIREN/injection/Weighter.h"
#include "../../public/SIREN/injection/WeightingUtils.h"

#include "../../../distributions/public/SIREN/distributions/primary/vertex/DepthFunction.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../detector/public/SIREN/detector/DetectorModel.h"
#include "../../../interactions/public/SIREN/interactions/InteractionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(injection,m) {
  using namespace siren::injection;

  // Utils function
    
  m.def("CrossSectionProbability", &CrossSectionProbability);
    
  // Process

  class_<Process, std::shared_ptr<Process>>(m, "Process")
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, &Process::SetInteractions);

  class_<PhysicalProcess, std::shared_ptr<PhysicalProcess>, Process>(m, "PhysicalProcess")
    .def(init<>())
    .def("AddPhysicalDistribution",&PhysicalProcess::AddPhysicalDistribution)
    .def("GetPhysicalDistributions",&PhysicalProcess::GetPhysicalDistributions);

  class_<PrimaryInjectionProcess, std::shared_ptr<PrimaryInjectionProcess>, Process>(m, "PrimaryInjectionProcess")
    .def(init<>())
    .def("AddPrimaryInjectionDistribution",&PrimaryInjectionProcess::AddPrimaryInjectionDistribution)
    .def("GetPrimaryInjectionDistributions",&PrimaryInjectionProcess::GetPrimaryInjectionDistributions);

  class_<SecondaryInjectionProcess, std::shared_ptr<SecondaryInjectionProcess>, Process>(m, "SecondaryInjectionProcess")
    .def(init<>())
    .def("AddSecondaryInjectionDistribution",&SecondaryInjectionProcess::AddSecondaryInjectionDistribution)
    .def("GetSecondaryInjectionDistributions",&SecondaryInjectionProcess::GetSecondaryInjectionDistributions);


  // Injection

  class_<Injector, std::shared_ptr<Injector>>(m, "Injector")
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>>())
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
    .def("GetPrimaryInjectionDistributions",&Injector::GetPrimaryInjectionDistributions)
    .def("GetDetectorModel",&Injector::GetDetectorModel)
    .def("GetInteractions",&Injector::GetInteractions)
    .def("InjectedEvents",&Injector::InjectedEvents)
    .def("EventsToInject",&Injector::EventsToInject);

//  class_<RangedSIREN, std::shared_ptr<RangedSIREN>, Injector>(m, "RangedSIREN")
//    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::distributions::RangeFunction>, double, double>())
//    .def("Name",&RangedSIREN::Name);

//  class_<DecayRangeSIREN, std::shared_ptr<DecayRangeSIREN>, Injector>(m, "DecayRangeSIREN")
//    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::distributions::DecayRangeFunction>, double, double>())
//    .def("Name",&DecayRangeSIREN::Name);
//
//  class_<ColumnDepthSIREN, std::shared_ptr<ColumnDepthSIREN>, Injector>(m, "ColumnDepthSIREN")
//    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::distributions::DepthFunction>, double, double>())
//    .def("Name",&ColumnDepthSIREN::Name)
//    .def("PrimaryInjectionBounds",&ColumnDepthSIREN::PrimaryInjectionBounds)
//    .def("SecondaryInjectionBounds",&ColumnDepthSIREN::SecondaryInjectionBounds);
//
//  class_<CylinderVolumeSIREN, std::shared_ptr<CylinderVolumeSIREN>, Injector>(m, "CylinderVolumeSIREN")
//    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>, siren::geometry::Cylinder>())
//    .def("Name",&CylinderVolumeSIREN::Name);

  // Weighter classes

  class_<PrimaryProcessWeighter, std::shared_ptr<PrimaryProcessWeighter>>(m, "PrimaryProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<PrimaryInjectionProcess>, std::shared_ptr<siren::detector::DetectorModel>>())
    .def("InteractionProbability",&PrimaryProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&PrimaryProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",&PrimaryProcessWeighter::PhysicalProbability)
    .def("GenerationProbability",&PrimaryProcessWeighter::GenerationProbability)
    .def("EventWeight",&PrimaryProcessWeighter::EventWeight);

  class_<SecondaryProcessWeighter, std::shared_ptr<SecondaryProcessWeighter>>(m, "SecondaryProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<SecondaryInjectionProcess>, std::shared_ptr<siren::detector::DetectorModel>>())
    .def("InteractionProbability",&SecondaryProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&SecondaryProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",&SecondaryProcessWeighter::PhysicalProbability)
    .def("GenerationProbability",&SecondaryProcessWeighter::GenerationProbability)
    .def("EventWeight",&SecondaryProcessWeighter::EventWeight);

  class_<LeptonTreeWeighter, std::shared_ptr<LeptonTreeWeighter>>(m, "LeptonTreeWeighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>, std::vector<std::shared_ptr<PhysicalProcess>>>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>>())
    .def("EventWeight",&LeptonTreeWeighter::EventWeight);

  class_<LeptonWeighter, std::shared_ptr<LeptonWeighter>>(m, "LeptonWeighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<siren::interactions::InteractionCollection>, std::vector<std::shared_ptr<siren::distributions::WeightableDistribution>>>())
    .def("EventWeight",&LeptonWeighter::EventWeight)
    .def("SimplifiedEventWeight",&LeptonWeighter::SimplifiedEventWeight);




}
