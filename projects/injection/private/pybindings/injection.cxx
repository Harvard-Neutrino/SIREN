
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
#include "../../public/SIREN/injection/DetectorDirectedAngularSectorChannel.h"
#include "../../public/SIREN/injection/DetectorDirected3BodyChannel.h"
#include "../../public/SIREN/injection/DetectorDirectedScatteringChannel.h"
#include "../../public/SIREN/injection/PhysicalChannelAdapters.h"
#include "../../public/SIREN/injection/TwoBodyKinematics.h"
#include "../../public/SIREN/injection/InvariantMassMapping.h"

#include "../../../geometry/public/SIREN/geometry/Geometry.h"

#include "../../../distributions/public/SIREN/distributions/primary/vertex/DepthFunction.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../utilities/public/SIREN/utilities/Errors.h"
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
  m.def("CrossSectionProbabilityWithPhaseSpace", &CrossSectionProbabilityWithPhaseSpace);
  m.def("ChannelSelectionProbability", &ChannelSelectionProbability);

  // Vertex weighting mode
  using VWM = siren::dataclasses::VertexWeightingMode;

  enum_<VWM::BoundSource>(m, "BoundSource")
    .value("Geometry", VWM::BoundSource::Geometry)
    .value("Distribution", VWM::BoundSource::Distribution)
    .value("None", VWM::BoundSource::None);

  class_<VWM>(m, "VertexWeightingMode")
    .def(init<>())
    .def_readwrite("compute_interaction_probability", &VWM::compute_interaction_probability)
    .def_readwrite("compute_position_probability", &VWM::compute_position_probability)
    .def_readwrite("bound_source", &VWM::bound_source)
    .def("__eq__", &VWM::operator==)
    .def("__ne__", &VWM::operator!=)
    .def_static("Propagated", &VWM::Propagated)
    .def_static("Fixed", &VWM::Fixed)
    .def_static("ExternalBounds", &VWM::ExternalBounds)
    ;

  // Phase space channels

  // New topology/measure enums
  enum_<PhaseSpaceTopology>(m, "PhaseSpaceTopology")
    .value("Decay2Body", PhaseSpaceTopology::Decay2Body)
    .value("Decay3Body", PhaseSpaceTopology::Decay3Body)
    .value("DecayNBody", PhaseSpaceTopology::DecayNBody)
    .value("Scatter2to2", PhaseSpaceTopology::Scatter2to2)
    .value("Scatter2to3", PhaseSpaceTopology::Scatter2to3)
    .value("Unspecified", PhaseSpaceTopology::Unspecified);

  enum_<siren::utilities::FailureReason>(m, "FailureReason")
    .value("Unspecified", siren::utilities::FailureReason::Unspecified)
    .value("NoPathThroughVolume", siren::utilities::FailureReason::NoPathThroughVolume)
    .value("NoTargetsOnPath", siren::utilities::FailureReason::NoTargetsOnPath)
    .value("NoColumnDepthSolution", siren::utilities::FailureReason::NoColumnDepthSolution)
    .value("KinematicallyForbidden", siren::utilities::FailureReason::KinematicallyForbidden)
    .value("UnregisteredSecondaryType", siren::utilities::FailureReason::UnregisteredSecondaryType)
    .value("PrimaryVertexFailure", siren::utilities::FailureReason::PrimaryVertexFailure)
    .value("TopLevelCatch", siren::utilities::FailureReason::TopLevelCatch);

  enum_<PhaseSpaceMeasure::Type>(m, "PhaseSpaceMeasureType")
    .value("SolidAngleRest", PhaseSpaceMeasure::Type::SolidAngleRest)
    .value("SolidAngleLab", PhaseSpaceMeasure::Type::SolidAngleLab)
    .value("Recursive2Body", PhaseSpaceMeasure::Type::Recursive2Body)
    .value("DalitzPair", PhaseSpaceMeasure::Type::DalitzPair)
    .value("HelicityAngles", PhaseSpaceMeasure::Type::HelicityAngles)
    .value("MandelstamQ2", PhaseSpaceMeasure::Type::MandelstamQ2)
    .value("BjorkenXY", PhaseSpaceMeasure::Type::BjorkenXY)
    .value("Unspecified", PhaseSpaceMeasure::Type::Unspecified);

  class_<PhaseSpaceMeasure>(m, "PhaseSpaceMeasure")
    .def(init<>())
    .def_readwrite("type", &PhaseSpaceMeasure::type)
    .def_readwrite("spectator", &PhaseSpaceMeasure::spectator)
    .def_readwrite("pair_first", &PhaseSpaceMeasure::pair_first)
    .def_readwrite("pair_second", &PhaseSpaceMeasure::pair_second)
    .def("__eq__", &PhaseSpaceMeasure::operator==)
    .def("__ne__", &PhaseSpaceMeasure::operator!=)
    .def("__hash__", [](PhaseSpaceMeasure const & m) {
        size_t h = std::hash<int>()(static_cast<int>(m.type));
        h ^= std::hash<int>()(m.spectator) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(m.pair_first) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(m.pair_second) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    })
    .def_static("SolidAngleRest", &PhaseSpaceMeasure::SolidAngleRest)
    .def_static("SolidAngleLab", &PhaseSpaceMeasure::SolidAngleLab)
    .def_static("Recursive2Body", &PhaseSpaceMeasure::Recursive2Body,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2)
    .def_static("DalitzPair", &PhaseSpaceMeasure::DalitzPair,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2)
    .def_static("HelicityAngles", &PhaseSpaceMeasure::HelicityAngles,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2)
    .def_static("MandelstamQ2", &PhaseSpaceMeasure::MandelstamQ2)
    .def_static("BjorkenXY", &PhaseSpaceMeasure::BjorkenXY)
    .def_static("Unspecified", &PhaseSpaceMeasure::Unspecified)
    ;

  m.def("PhaseSpaceTopologyName", &PhaseSpaceTopologyName);
  m.def("PhaseSpaceMeasureName", &PhaseSpaceMeasureName);

  // Convert a sampling density between phase-space measures within a topology,
  // applying the analytic Jacobian (the same conversion the mixture uses).
  m.def("ConvertDensity",
        (double (*)(double, PhaseSpaceMeasure const &, PhaseSpaceMeasure const &,
                    PhaseSpaceTopology, siren::dataclasses::InteractionRecord const &))
            &siren::injection::ConvertDensity,
        arg("density"), arg("from_measure"), arg("to_measure"),
        arg("topology"), arg("record"));

  // Legacy convention enum (kept for backward compatibility)
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
    .def("Topology", &PhaseSpaceChannel::Topology)
    .def("Measure", &PhaseSpaceChannel::Measure)
    .def("Convention", &PhaseSpaceChannel::Convention)
    ;

  class_<MultiChannelPhaseSpace, std::shared_ptr<MultiChannelPhaseSpace>> multi_channel_phase_space(m, "MultiChannelPhaseSpace");

  // Severity-tagged compatibility diagnostic returned by ValidateChannelsDetailed,
  // bound as a nested type so the ValidateChannelsDetailed return type crosses.
  {
    class_<MultiChannelPhaseSpace::ChannelDiagnostic> channel_diagnostic(
        multi_channel_phase_space, "ChannelDiagnostic");
    enum_<MultiChannelPhaseSpace::ChannelDiagnostic::Severity>(channel_diagnostic, "Severity")
      .value("Info", MultiChannelPhaseSpace::ChannelDiagnostic::Severity::Info)
      .value("Fatal", MultiChannelPhaseSpace::ChannelDiagnostic::Severity::Fatal);
    channel_diagnostic
      .def_readonly("severity", &MultiChannelPhaseSpace::ChannelDiagnostic::severity)
      .def_readonly("message", &MultiChannelPhaseSpace::ChannelDiagnostic::message);
  }

  multi_channel_phase_space
    .def(init<>())
    .def(init<std::vector<std::shared_ptr<PhaseSpaceChannel>>, std::vector<double>, bool>(),
         arg("channels"), arg("weights") = std::vector<double>{},
         arg("allow_incompatible") = false)
    .def_readwrite("channels", &MultiChannelPhaseSpace::channels)
    .def_readwrite("weights", &MultiChannelPhaseSpace::weights)
    .def("Normalize", &MultiChannelPhaseSpace::Normalize)
    .def("Sample", &MultiChannelPhaseSpace::Sample)
    .def("Density", &MultiChannelPhaseSpace::Density)
    .def("DensityBreakdown", &MultiChannelPhaseSpace::DensityBreakdown,
         arg("detector_model"), arg("record"))
    .def("Accumulate", &MultiChannelPhaseSpace::Accumulate,
         arg("detector_model"), arg("record"), arg("weight"),
         arg("discount_fallback") = true, arg("recurse") = true)
    .def("UpdateWeights", &MultiChannelPhaseSpace::UpdateWeights,
         arg("update_rule"), arg("damping"), arg("min_weight"),
         arg("recurse") = true, arg("failure_mode") = "throughput")
    .def("ResetAccumulators", &MultiChannelPhaseSpace::ResetAccumulators,
         arg("recurse") = true)
    .def("AccumulateSelection", &MultiChannelPhaseSpace::AccumulateSelection,
         arg("detector_model"), arg("record"), arg("failed"))
    .def_readwrite("kp_accumulator", &MultiChannelPhaseSpace::kp_accumulator_)
    .def_readwrite("kp_count", &MultiChannelPhaseSpace::kp_count_)
    .def_readwrite("kp_succ_select", &MultiChannelPhaseSpace::kp_succ_select_)
    .def_readwrite("kp_fail_select", &MultiChannelPhaseSpace::kp_fail_select_)
    .def("CommonTopology", &MultiChannelPhaseSpace::CommonTopology)
    .def("CommonMeasure", &MultiChannelPhaseSpace::CommonMeasure)
    .def("CommonConvention", &MultiChannelPhaseSpace::CommonConvention)
    .def("ValidateChannels", &MultiChannelPhaseSpace::ValidateChannels)
    .def("ValidateChannelsDetailed", &MultiChannelPhaseSpace::ValidateChannelsDetailed)
    .def("ValidateChannelDensities", &MultiChannelPhaseSpace::ValidateChannelDensities,
         arg("random"), arg("detector_model"), arg("template_record"),
         arg("samples_per_channel") = 100)
    ;

  // A sub-mixture wrapped as one channel: encapsulates a set of (e.g.
  // geometric) channels whose inner weights remain part of the optimization.
  class_<NestedMixtureChannel, std::shared_ptr<NestedMixtureChannel>, PhaseSpaceChannel>(m, "NestedMixtureChannel")
    .def(init<std::shared_ptr<MultiChannelPhaseSpace>>(), arg("mixture"))
    .def_readwrite("mixture", &NestedMixtureChannel::mixture)
    .def_readwrite("label", &NestedMixtureChannel::label)
    .def("DirectingActive", &NestedMixtureChannel::DirectingActive, arg("record"))
    ;

  class_<Isotropic2BodyChannel, std::shared_ptr<Isotropic2BodyChannel>, PhaseSpaceChannel>(m, "Isotropic2BodyChannel")
    .def(init<int>(), arg("daughter_index") = 0)
    ;

  enum_<DetectorDirected2BodyChannel::Mode>(m, "DirectedMode")
    .value("Cone", DetectorDirected2BodyChannel::Mode::Cone)
    .value("Volume", DetectorDirected2BodyChannel::Mode::Volume);

  class_<DetectorDirected2BodyChannel, std::shared_ptr<DetectorDirected2BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected2BodyChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int, DetectorDirected2BodyChannel::Mode, double>(),
         arg("target"), arg("daughter_index") = 0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("volume") = -1.0)
    .def("SetVolume", &DetectorDirected2BodyChannel::SetVolume)
    .def("DirectingActive", &DetectorDirected2BodyChannel::DirectingActive, arg("record"))
    ;

  class_<DetectorDirectedAngularSectorChannel, std::shared_ptr<DetectorDirectedAngularSectorChannel>, PhaseSpaceChannel>(m, "DetectorDirectedAngularSectorChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, double, double, double, double, int>(),
         arg("target"), arg("u_lo"), arg("u_hi"), arg("phi_lo"), arg("phi_hi"),
         arg("daughter_index") = 0)
    .def("DirectingActive", &DetectorDirectedAngularSectorChannel::DirectingActive, arg("record"))
    ;

  // 1-D importance maps: one object provides BOTH the draw (Forward) and
  // its own normalized density (Density), so a model/channel routing both
  // through a shared instance cannot let sampling and density drift apart
  // (Contract C1).  Models consume these; they do not subclass in Python,
  // so no trampoline is needed.
  class_<Mapping1D, std::shared_ptr<Mapping1D>>(m, "Mapping1D")
    .def("Forward", &Mapping1D::Forward, arg("r"))
    .def("Inverse", &Mapping1D::Inverse, arg("x"))
    .def("Density", &Mapping1D::Density, arg("x"))
    .def("Accumulate", &Mapping1D::Accumulate, arg("x"), arg("weight"))
    .def("Refine", &Mapping1D::Refine);

  class_<BreitWignerMapping, std::shared_ptr<BreitWignerMapping>, Mapping1D>(m, "BreitWignerMapping")
    .def(init<double, double, double, double>(),
         arg("mass"), arg("width"), arg("s_min"), arg("s_max"));

  class_<PowerLawMapping, std::shared_ptr<PowerLawMapping>, Mapping1D>(m, "PowerLawMapping")
    .def(init<double, double, double, double>(),
         arg("nu"), arg("m2"), arg("s_min"), arg("s_max"));

  class_<TabulatedMapping, std::shared_ptr<TabulatedMapping>, Mapping1D>(m, "TabulatedMapping")
    .def(init<std::vector<double>, std::vector<double>, double, double>(),
         arg("s_nodes"), arg("cdf_nodes"), arg("s_min"), arg("s_max"));

  class_<PropagatorMapping, std::shared_ptr<PropagatorMapping>, Mapping1D>(m, "PropagatorMapping")
    .def(init<double, double, double>(),
         arg("m2"), arg("x_min"), arg("x_max"));

  class_<UniformMapping, std::shared_ptr<UniformMapping>, Mapping1D>(m, "UniformMapping")
    .def(init<double, double>(),
         arg("s_min"), arg("s_max"));

  class_<LogMapping, std::shared_ptr<LogMapping>, Mapping1D>(m, "LogMapping")
    .def(init<double, double>(),
         arg("x_min"), arg("x_max"));

  class_<ExponentialMapping, std::shared_ptr<ExponentialMapping>, Mapping1D>(m, "ExponentialMapping")
    .def(init<double, double, double>(),
         arg("tau"), arg("x_min"), arg("x_max"));

  class_<GaussianMapping, std::shared_ptr<GaussianMapping>, Mapping1D>(m, "GaussianMapping")
    .def(init<double, double, double, double>(),
         arg("mu"), arg("sigma"), arg("x_min"), arg("x_max"));

  class_<AdaptiveMapping, std::shared_ptr<AdaptiveMapping>, Mapping1D>(m, "AdaptiveMapping")
    .def(init<double, double, int, double, double>(),
         arg("x_min"), arg("x_max"), arg("n_bins") = 32,
         arg("damping") = 0.5, arg("floor_frac") = 1e-3)
    .def_readonly("p", &AdaptiveMapping::p)
    .def_readonly("n_bins", &AdaptiveMapping::n_bins);

  enum_<DetectorDirected3BodyChannel::InvariantMassMode>(m, "InvariantMassMode")
    .value("Uniform", DetectorDirected3BodyChannel::InvariantMassMode::Uniform)
    .value("BreitWigner", DetectorDirected3BodyChannel::InvariantMassMode::BreitWigner)
    .value("PowerLaw", DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw)
    .value("Tabulated", DetectorDirected3BodyChannel::InvariantMassMode::Tabulated);

  class_<DetectorDirected3BodyChannel, std::shared_ptr<DetectorDirected3BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected3BodyChannel")
    // Direct mode: specify which daughter to bias, others inferred.
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int,
              DetectorDirected3BodyChannel::InvariantMassMode, double, double, double, double,
              DetectorDirected2BodyChannel::Mode, PhaseSpaceTopology,
              std::vector<double>, std::vector<double>>(),
         arg("target"), arg("directed_index"),
         arg("mass_mode") = DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
         arg("resonance_mass") = 0.0, arg("resonance_width") = 0.0,
         arg("power_law_nu") = 0.8, arg("power_law_offset") = 0.0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("topology") = PhaseSpaceTopology::Decay3Body,
         arg("mass_cdf_nodes") = std::vector<double>{},
         arg("mass_cdf_values") = std::vector<double>{})
    // Recursive mode: specify all indices (backward compatible).
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int, int, int, int,
              DetectorDirected3BodyChannel::InvariantMassMode, double, double, double, double,
              DetectorDirected2BodyChannel::Mode, PhaseSpaceTopology,
              std::vector<double>, std::vector<double>>(),
         arg("target"), arg("spectator_index"),
         arg("pair_first_index"), arg("pair_second_index"),
         arg("directed_pair_index"),
         arg("mass_mode") = DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
         arg("resonance_mass") = 0.0, arg("resonance_width") = 0.0,
         arg("power_law_nu") = 0.8, arg("power_law_offset") = 0.0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("topology") = PhaseSpaceTopology::Decay3Body,
         arg("mass_cdf_nodes") = std::vector<double>{},
         arg("mass_cdf_values") = std::vector<double>{})
    .def("SetVolume", &DetectorDirected3BodyChannel::SetVolume)
    .def("DirectingActive", &DetectorDirected3BodyChannel::DirectingActive, arg("record"))
    ;

  enum_<DetectorDirectedScatteringChannel::Variable>(m, "ScatteringVariable")
    .value("Q2", DetectorDirectedScatteringChannel::Variable::Q2)
    .value("BjorkenY", DetectorDirectedScatteringChannel::Variable::BjorkenY)
    .value("RecoilY", DetectorDirectedScatteringChannel::Variable::RecoilY);

  enum_<DetectorDirectedScatteringChannel::Q2Mode>(m, "ScatteringQ2Mode")
    .value("Geometry", DetectorDirectedScatteringChannel::Q2Mode::Geometry)
    .value("Propagator", DetectorDirectedScatteringChannel::Q2Mode::Propagator)
    .value("Tabulated", DetectorDirectedScatteringChannel::Q2Mode::Tabulated);

  class_<DetectorDirectedScatteringChannel, std::shared_ptr<DetectorDirectedScatteringChannel>, PhaseSpaceChannel>(m, "DetectorDirectedScatteringChannel")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int,
              DetectorDirectedScatteringChannel::Variable,
              DetectorDirected2BodyChannel::Mode,
              DetectorDirectedScatteringChannel::Q2Mode, double,
              std::vector<double>, std::vector<double>>(),
         arg("target"), arg("directed_index") = 0,
         arg("variable") = DetectorDirectedScatteringChannel::Variable::Q2,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("q2_mode") = DetectorDirectedScatteringChannel::Q2Mode::Geometry,
         arg("mediator_mass") = 0.0,
         arg("q2_cdf_nodes") = std::vector<double>{},
         arg("q2_cdf_values") = std::vector<double>{})
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
    .def("SetWeightingMode", &PhysicalProcess::SetWeightingMode)
    .def("GetWeightingMode", &PhysicalProcess::GetWeightingMode)
    .def_property("weighting_mode", &PhysicalProcess::GetWeightingMode, &PhysicalProcess::SetWeightingMode)
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

  class_<FailureLedger, std::shared_ptr<FailureLedger>>(m, "FailureLedger")
    .def(init<>())
    .def("Clear", &FailureLedger::Clear)
    .def("entries", [](FailureLedger const & ledger) {
        pybind11::dict out;
        for(auto const & item : ledger.entries) {
            pybind11::tuple key = pybind11::make_tuple(
                item.first.depth, item.first.parent_pdg, item.first.reason);
            pybind11::tuple value = pybind11::make_tuple(
                item.second.count, item.second.exemplar);
            out[key] = value;
        }
        return out;
    });

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
    .def("FailedEvents",&Injector::FailedEvents)
    .def("GetFailureCounts",&Injector::GetFailureCounts)
    .def("GetLastFailureReason",&Injector::GetLastFailureReason)
    .def("GetLastFailedTree",&Injector::GetLastFailedTree, pybind11::return_value_policy::reference_internal)
    .def("GetFailureLedger",&Injector::GetFailureLedger, pybind11::return_value_policy::reference_internal)
    .def("ResetInjectedEvents",overload_cast<unsigned int>(&Injector::ResetInjectedEvents))
    .def("ResetInjectedEvents",overload_cast<>(&Injector::ResetInjectedEvents))
    .def("GetPhaseSpaces",&Injector::GetPhaseSpaces)
    .def("AccumulateEventToMixtures",&Injector::AccumulateEventToMixtures,
         arg("tree"), arg("weight"), arg("discount_fallback") = true,
         arg("recurse") = false)
    .def("AccumulateSelectionToMixtures",&Injector::AccumulateSelectionToMixtures,
         arg("tree"), arg("failed"))
    .def("PrimaryInjectionBounds",&Injector::PrimaryInjectionBounds)
    .def("SecondaryInjectionBounds",&Injector::SecondaryInjectionBounds)
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

  class_<VertexWeightFactors, std::shared_ptr<VertexWeightFactors>>(m, "VertexWeightFactors")
    .def_readonly("depth", &VertexWeightFactors::depth)
    .def_readonly("primary_pdg", &VertexWeightFactors::primary_pdg)
    .def_readonly("generation", &VertexWeightFactors::generation)
    .def_readonly("physical", &VertexWeightFactors::physical)
    .def_readonly("interaction_prob", &VertexWeightFactors::interaction_prob)
    .def_readonly("position_prob", &VertexWeightFactors::position_prob)
    .def_readonly("channel_densities", &VertexWeightFactors::channel_densities)
    .def_readonly("cancelled", &VertexWeightFactors::cancelled)
    .def_readonly("flags", &VertexWeightFactors::flags);

  class_<EventWeightBreakdown, std::shared_ptr<EventWeightBreakdown>>(m, "EventWeightBreakdown")
    .def_readonly("total", &EventWeightBreakdown::total)
    .def_readonly("vertices", &EventWeightBreakdown::vertices);

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
    .def("EventWeightWithBreakdown",&Weighter::EventWeightWithBreakdown)
    .def("GetInjectors",&Weighter::GetInjectors)
    .def("GetDetectorModel",&Weighter::GetDetectorModel)
    .def("GetPrimaryPhysicalProcess",&Weighter::GetPrimaryPhysicalProcess)
    .def("GetSecondaryPhysicalProcesses",&Weighter::GetSecondaryPhysicalProcesses)
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
