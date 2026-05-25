
#include <vector>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/injection/Process.h"
#include "../../public/SIREN/injection/Injector.h"
#include "../../public/SIREN/injection/Weighter.h"
#include "../../public/SIREN/injection/WeightingUtils.h"
#include "../../public/SIREN/injection/PhaseSpaceChannel.h"
#include "../../public/SIREN/injection/Isotropic2BodyChannel.h"
#include "../../public/SIREN/injection/DetectorDirected2BodyChannel.h"
#include "../../public/SIREN/injection/DetectorDirected3BodyChannel.h"
#include "../../public/SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "../../public/SIREN/injection/PhysicalChannelAdapters.h"
#include "../../public/SIREN/injection/TwoBodyKinematics.h"

#include "../../../geometry/public/SIREN/geometry/Geometry.h"

#include "../../../distributions/public/SIREN/distributions/primary/vertex/DepthFunction.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../detector/public/SIREN/detector/DetectorModel.h"
#include "../../../interactions/public/SIREN/interactions/InteractionCollection.h"

#include "../../../interactions/public/SIREN/interactions/pyDarkNewsCrossSection.h"
#include "../../../interactions/public/SIREN/interactions/pyDarkNewsDecay.h"

#include "../../../serialization/public/SIREN/serialization/ByteString.h"

#include "SIREN/dataclasses/serializable.h"
#include "SIREN/detector/serializable.h"
#include "SIREN/distributions/serializable.h"
#include "SIREN/geometry/serializable.h"
#include "SIREN/injection/serializable.h"
#include "SIREN/interactions/serializable.h"
#include "SIREN/math/serializable.h"
#include "SIREN/utilities/serializable.h"

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>);

using namespace pybind11;

PYBIND11_MODULE(injection,m) {
  using namespace siren::injection;

  // Utils function

  m.def("CrossSectionProbability", &CrossSectionProbability);

  // Phase space channels

  enum_<PhaseSpaceConvention>(m, "PhaseSpaceConvention")
    .value("RestFrameSolidAngle", PhaseSpaceConvention::RestFrameSolidAngle)
    .value("LabFrameSolidAngle", PhaseSpaceConvention::LabFrameSolidAngle)
    .value("Recursive2Body", PhaseSpaceConvention::Recursive2Body)
    .value("Dalitz", PhaseSpaceConvention::Dalitz)
    .value("HelicityAngles", PhaseSpaceConvention::HelicityAngles)
    .value("BjorkenXY", PhaseSpaceConvention::BjorkenXY)
    .value("MandelstamST", PhaseSpaceConvention::MandelstamST)
    .value("Custom", PhaseSpaceConvention::Custom);

  m.def("PhaseSpaceConventionName", &PhaseSpaceConventionName);

  class_<PhaseSpaceChannel, std::shared_ptr<PhaseSpaceChannel>>(m, "PhaseSpaceChannel")
    .def("Sample", &PhaseSpaceChannel::Sample)
    .def("Density", &PhaseSpaceChannel::Density)
    .def("Name", &PhaseSpaceChannel::Name)
    .def("Convention", &PhaseSpaceChannel::Convention)
    ;

  class_<MultiChannelPhaseSpace, std::shared_ptr<MultiChannelPhaseSpace>>(m, "MultiChannelPhaseSpace")
    .def(init<>())
    .def_readwrite("channels", &MultiChannelPhaseSpace::channels)
    .def_readwrite("weights", &MultiChannelPhaseSpace::weights)
    .def("Sample", &MultiChannelPhaseSpace::Sample)
    .def("Density", &MultiChannelPhaseSpace::Density)
    .def("CommonConvention", &MultiChannelPhaseSpace::CommonConvention)
    .def("ValidateConventions", &MultiChannelPhaseSpace::ValidateConventions)
    .def("ValidateChannels", &MultiChannelPhaseSpace::ValidateChannels,
         arg("random"), arg("detector_model"), arg("template_record"),
         arg("samples_per_channel") = 100)
    ;

  class_<Isotropic2BodyChannel, std::shared_ptr<Isotropic2BodyChannel>, PhaseSpaceChannel>(m, "Isotropic2BodyChannel")
    .def(init<int>(), arg("daughter_index") = 0)
    ;

  enum_<DetectorDirected2BodyChannel::Mode>(m, "DirectedMode")
    .value("Cone", DetectorDirected2BodyChannel::Mode::Cone)
    .value("Volume", DetectorDirected2BodyChannel::Mode::Volume);

  class_<DetectorDirected2BodyChannel, std::shared_ptr<DetectorDirected2BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected2BodyChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int, DetectorDirected2BodyChannel::Mode>(),
         arg("target"), arg("daughter_index") = 0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume)
    .def("SetVolume", &DetectorDirected2BodyChannel::SetVolume)
    ;

  enum_<DetectorDirected3BodyChannel::InvariantMassMode>(m, "InvariantMassMode")
    .value("Uniform", DetectorDirected3BodyChannel::InvariantMassMode::Uniform)
    .value("BreitWigner", DetectorDirected3BodyChannel::InvariantMassMode::BreitWigner)
    .value("PowerLaw", DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw);

  class_<DetectorDirected3BodyChannel, std::shared_ptr<DetectorDirected3BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected3BodyChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int, int, int, int,
              DetectorDirected3BodyChannel::InvariantMassMode, double, double, double, double,
              DetectorDirected2BodyChannel::Mode>(),
         arg("target"), arg("spectator_index") = 0,
         arg("pair_first_index") = 1, arg("pair_second_index") = 2,
         arg("directed_pair_index") = 1,
         arg("mass_mode") = DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
         arg("resonance_mass") = 0.0, arg("resonance_width") = 0.0,
         arg("power_law_nu") = 0.8, arg("power_law_offset") = 0.0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume)
    .def("SetVolume", &DetectorDirected3BodyChannel::SetVolume)
    ;

  enum_<DetectorDirectedScatteringChannel::Variable>(m, "ScatteringVariable")
    .value("Q2", DetectorDirectedScatteringChannel::Variable::Q2)
    .value("BjorkenY", DetectorDirectedScatteringChannel::Variable::BjorkenY)
    .value("RecoilY", DetectorDirectedScatteringChannel::Variable::RecoilY);

  class_<DetectorDirectedScatteringChannel, std::shared_ptr<DetectorDirectedScatteringChannel>, PhaseSpaceChannel>(m, "DetectorDirectedScatteringChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int,
              DetectorDirectedScatteringChannel::Variable,
              DetectorDirected2BodyChannel::Mode>(),
         arg("target"), arg("directed_index") = 0,
         arg("variable") = DetectorDirectedScatteringChannel::Variable::Q2,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume)
    .def("SetVolume", &DetectorDirectedScatteringChannel::SetVolume)
    ;

  // Physical channel adapters

  class_<PhysicalDecayChannel, std::shared_ptr<PhysicalDecayChannel>, PhaseSpaceChannel>(m, "PhysicalDecayChannel")
    .def(init<std::shared_ptr<siren::interactions::Decay>>())
    .def(init<std::shared_ptr<siren::interactions::Decay>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::Decay>, PhaseSpaceConvention>())
    .def("GetDecay", &PhysicalDecayChannel::GetDecay)
    ;

  class_<PhysicalCrossSectionChannel, std::shared_ptr<PhysicalCrossSectionChannel>, PhaseSpaceChannel>(m, "PhysicalCrossSectionChannel")
    .def(init<std::shared_ptr<siren::interactions::CrossSection>>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>, PhaseSpaceConvention>())
    .def("GetCrossSection", &PhysicalCrossSectionChannel::GetCrossSection)
    ;

  // Two-body kinematics utilities

  m.def("TwoBodyRestMomentum", &TwoBodyRestMomentum);
  m.def("TwoBodyRestEnergy", &TwoBodyRestEnergy);
  m.def("Kallen", &Kallen);

  // Process

  class_<Process, std::shared_ptr<Process>>(m, "Process")
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, &Process::SetInteractions)
    ;

  class_<PhysicalProcess, std::shared_ptr<PhysicalProcess>, Process>(m, "PhysicalProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>())
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, &Process::SetInteractions)
    .def_property("distributions", &PhysicalProcess::GetPhysicalDistributions, &PhysicalProcess::SetPhysicalDistributions)
    .def("SetPhaseSpace", &PhysicalProcess::SetPhaseSpace)
    .def("GetPhaseSpace", &PhysicalProcess::GetPhaseSpace)
    .def("HasPhaseSpace", overload_cast<siren::dataclasses::InteractionSignature const &>(&PhysicalProcess::HasPhaseSpace, const_))
    .def("HasAnyPhaseSpace", &PhysicalProcess::HasAnyPhaseSpace)
    ;

  class_<PrimaryInjectionProcess, std::shared_ptr<PrimaryInjectionProcess>, PhysicalProcess>(m, "PrimaryInjectionProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>())
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, &Process::SetInteractions)
    .def_property("distributions", &PrimaryInjectionProcess::GetPrimaryInjectionDistributions, &PrimaryInjectionProcess::SetPrimaryInjectionDistributions)
    ;

  class_<SecondaryInjectionProcess, std::shared_ptr<SecondaryInjectionProcess>, PhysicalProcess>(m, "SecondaryInjectionProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>())
    .def_property("secondary_type", &SecondaryInjectionProcess::GetSecondaryType, &SecondaryInjectionProcess::SetSecondaryType)
    .def_property("interactions", &Process::GetInteractions, &Process::SetInteractions)
    .def_property("distributions", &SecondaryInjectionProcess::GetSecondaryInjectionDistributions, &SecondaryInjectionProcess::SetSecondaryInjectionDistributions)
    ;

  // Injection

  class_<Injector, std::shared_ptr<Injector>>(m, "Injector")
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def(init<unsigned int, std::string, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>>())
    .def("SetStoppingCondition",&Injector::SetStoppingCondition)
    .def("GetStoppingCondition",&Injector::GetStoppingCondition)
    .def("SetPrimaryProcess",&Injector::SetPrimaryProcess)
    .def("AddSecondaryProcess",&Injector::AddSecondaryProcess)
    .def("GetPrimaryProcess",&Injector::GetPrimaryProcess)
    .def("GetSecondaryProcesses",&Injector::GetSecondaryProcesses)
    .def("GetSecondaryProcessMap",&Injector::GetSecondaryProcessMap)
    .def("NewRecord",&Injector::NewRecord)
    .def("SetRandom",&Injector::SetRandom)
    .def("GenerateEvent",&Injector::GenerateEvent, pybind11::return_value_policy::move)
    .def("DensityVariables",&Injector::DensityVariables)
    .def("Name",&Injector::Name)
    .def("GetPrimaryInjectionDistributions",&Injector::GetPrimaryInjectionDistributions)
    .def("GetDetectorModel",&Injector::GetDetectorModel)
    .def("SetDetectorModel",&Injector::SetDetectorModel)
    .def("GetInteractions",&Injector::GetInteractions)
    .def("InjectedEvents",&Injector::InjectedEvents)
    .def("InjectionAttempts",&Injector::InjectionAttempts)
    .def("EventsToInject",&Injector::EventsToInject)
    .def("__len__", &Injector::EventsToInject)
    .def("ResetInjectedEvents",&Injector::ResetInjectedEvents)
    .def("SaveInjector",&Injector::SaveInjector)
    .def("LoadInjector",&Injector::LoadInjector)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<Injector>),
        &(siren::serialization::pickle_load<Injector>)
    ))
    ;

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
    .def("EventWeight",&PrimaryProcessWeighter::EventWeight)
    ;

  class_<SecondaryProcessWeighter, std::shared_ptr<SecondaryProcessWeighter>>(m, "SecondaryProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<SecondaryInjectionProcess>, std::shared_ptr<siren::detector::DetectorModel>>())
    .def("InteractionProbability",&SecondaryProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&SecondaryProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",&SecondaryProcessWeighter::PhysicalProbability)
    .def("GenerationProbability",&SecondaryProcessWeighter::GenerationProbability)
    .def("EventWeight",&SecondaryProcessWeighter::EventWeight)
    ;

  class_<Weighter, std::shared_ptr<Weighter>>(m, "Weighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>, std::vector<std::shared_ptr<PhysicalProcess>>>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::string>())
    .def("EventWeight",&Weighter::EventWeight)
    .def("GetInteractionProbabilities",&Weighter::GetInteractionProbabilities, arg("tree"), arg("i_inj")=0)
    .def("GetSurvivalProbabilities",&Weighter::GetSurvivalProbabilities, arg("tree"), arg("i_inj")=0)
    .def("SaveWeighter",&Weighter::SaveWeighter)
    .def("LoadWeighter",&Weighter::LoadWeighter)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<Weighter>),
        &(siren::serialization::pickle_load<Weighter>)
    ))
    ;
}
