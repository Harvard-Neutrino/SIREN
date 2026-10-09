
#include <cmath>
#include <stdexcept>
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
#include "../../public/SIREN/injection/PhysicalChannelAdapters.h"
#include "../../public/SIREN/injection/TwoBodyKinematics.h"

#include "../../../geometry/public/SIREN/geometry/Geometry.h"

#include "../../../distributions/public/SIREN/distributions/primary/vertex/DepthFunction.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../utilities/public/SIREN/utilities/Errors.h"
#include "../../../detector/public/SIREN/detector/DetectorModel.h"
#include "../../../interactions/public/SIREN/interactions/InteractionCollection.h"
#include "../../../interactions/public/SIREN/interactions/CharmMesonDecay.h"

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

  m.def("CrossSectionProbability",
        overload_cast<
            std::shared_ptr<siren::detector::DetectorModel const>,
            std::shared_ptr<siren::interactions::InteractionCollection const>,
            siren::dataclasses::InteractionRecord const &>(
            &CrossSectionProbability));
  m.def("CrossSectionProbabilityWithPhaseSpace",
        overload_cast<
            std::shared_ptr<siren::detector::DetectorModel const>,
            std::shared_ptr<siren::interactions::InteractionCollection const>,
            siren::dataclasses::InteractionRecord const &,
            MultiChannelPhaseSpace const &>(
            &CrossSectionProbabilityWithPhaseSpace));
  m.def("ChannelSelectionProbability", &ChannelSelectionProbability);
  m.def("FixedVertexChannelSelectionProbability", &FixedVertexChannelSelectionProbability,
        "Channel-selection factor a Fixed vertex charges when multiple channels "
        "compete: exactly 1.0 for a single candidate signature, else "
        "selected_rate/total_rate. Throws WeightCalculationError on a "
        "non-positive total rate or non-finite result.");

  // Vertex weighting mode
  using VWM = siren::dataclasses::VertexWeightingMode;

  enum_<VWM::BoundSource>(m, "BoundSource")
    .value("Geometry", VWM::BoundSource::Geometry)
    .value("Distribution", VWM::BoundSource::Distribution)
    // Exposed as "Unbounded": the C++ enumerator is None, but None is a Python
    // keyword and unreachable through attribute access.
    .value("Unbounded", VWM::BoundSource::None);

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

  enum_<PhaseSpaceTopology>(m, "PhaseSpaceTopology")
    .value("Decay2Body", PhaseSpaceTopology::Decay2Body,
           "Two-body final state from a single parent.")
    .value("Decay3Body", PhaseSpaceTopology::Decay3Body,
           "Three-body final state from a single parent.")
    .value("DecayNBody", PhaseSpaceTopology::DecayNBody,
           "N-body decay final state.")
    .value("Scatter2to2", PhaseSpaceTopology::Scatter2to2,
           "2->2 scattering (two incoming, two outgoing).")
    .value("Scatter2to3", PhaseSpaceTopology::Scatter2to3,
           "2->3 scattering.")
    .value("Unspecified", PhaseSpaceTopology::Unspecified,
           "Topology not declared; blocks mixing with typed channels.");

  enum_<siren::utilities::FailureReason>(m, "FailureReason")
    .value("Unspecified", siren::utilities::FailureReason::Unspecified)
    .value("NoPathThroughVolume", siren::utilities::FailureReason::NoPathThroughVolume)
    .value("NoTargetsOnPath", siren::utilities::FailureReason::NoTargetsOnPath)
    .value("NoColumnDepthSolution", siren::utilities::FailureReason::NoColumnDepthSolution)
    .value("KinematicallyForbidden", siren::utilities::FailureReason::KinematicallyForbidden);

  enum_<PhaseSpaceMeasure::Type>(m, "PhaseSpaceMeasureType")
    .value("CosThetaRest", PhaseSpaceMeasure::Type::CosThetaRest)
    .value("SolidAngleRest", PhaseSpaceMeasure::Type::SolidAngleRest,
           "Rest-frame solid angle.")
    .value("SolidAngleLab", PhaseSpaceMeasure::Type::SolidAngleLab,
           "Lab-frame solid angle.")
    .value("Recursive2Body", PhaseSpaceMeasure::Type::Recursive2Body,
           "Recursive two-body decomposition with spectator/pair indices.")
    .value("DalitzPair", PhaseSpaceMeasure::Type::DalitzPair,
           "Dalitz variables over a chosen pair.")
    .value("HelicityAngles", PhaseSpaceMeasure::Type::HelicityAngles,
           "Helicity-frame angles.")
    .value("MandelstamQ2", PhaseSpaceMeasure::Type::MandelstamQ2,
           "Momentum-transfer Q^2 with a uniform azimuth integrated out.")
    .value("MandelstamQ2Phi", PhaseSpaceMeasure::Type::MandelstamQ2Phi,
           "Momentum-transfer Q^2 and explicit beam-axis azimuth.")
    .value("BjorkenXY", PhaseSpaceMeasure::Type::BjorkenXY,
           "Bjorken x,y with a uniform azimuth integrated out.")
    .value("BjorkenXYPhi", PhaseSpaceMeasure::Type::BjorkenXYPhi,
           "Bjorken x,y and explicit beam-axis azimuth.")
    .value("FixedMassY", PhaseSpaceMeasure::Type::FixedMassY,
           "Fixed-mass y with a uniform azimuth integrated out.")
    .value("FixedMassYPhi", PhaseSpaceMeasure::Type::FixedMassYPhi,
           "Fixed-mass y and explicit beam-axis azimuth.")
    .value("MandelstamQ2Y", PhaseSpaceMeasure::Type::MandelstamQ2Y,
           "Momentum-transfer Q^2 and y with a uniform azimuth integrated out.")
    .value("MandelstamQ2YPhi", PhaseSpaceMeasure::Type::MandelstamQ2YPhi,
           "Momentum-transfer Q^2, y, and explicit beam-axis azimuth.")
    .value("Unspecified", PhaseSpaceMeasure::Type::Unspecified,
           "No measure declared.");

  class_<PhaseSpaceMeasure>(m, "PhaseSpaceMeasure")
    .def(init<>())
    .def_readwrite("type", &PhaseSpaceMeasure::type)
    .def_readwrite("spectator", &PhaseSpaceMeasure::spectator)
    .def_readwrite("pair_first", &PhaseSpaceMeasure::pair_first)
    .def_readwrite("pair_second", &PhaseSpaceMeasure::pair_second)
    .def("__eq__", &PhaseSpaceMeasure::operator==)
    .def("__ne__", &PhaseSpaceMeasure::operator!=)
    .def("__hash__", [](PhaseSpaceMeasure const & m) {
        // Mirror operator==: the factorization indices only participate in
        // equality (and therefore in the hash) for the index-relevant types.
        size_t h = std::hash<int>()(static_cast<int>(m.type));
        bool indices_relevant =
            m.type == PhaseSpaceMeasure::Type::Recursive2Body ||
            m.type == PhaseSpaceMeasure::Type::DalitzPair ||
            m.type == PhaseSpaceMeasure::Type::HelicityAngles;
        if (indices_relevant) {
            h ^= std::hash<int>()(m.spectator) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int>()(m.pair_first) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int>()(m.pair_second) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    })
    .def_static("CosThetaRest", &PhaseSpaceMeasure::CosThetaRest,
         "Rest-frame cos(theta) measure with uniform azimuth integrated out.")
    .def_static("SolidAngleRest", &PhaseSpaceMeasure::SolidAngleRest,
         "Rest-frame solid angle measure.")
    .def_static("SolidAngleLab", &PhaseSpaceMeasure::SolidAngleLab,
         arg("daughter_index") = 0,
         "Lab-frame solid angle measure of the indexed daughter.")
    .def_static("Recursive2Body", &PhaseSpaceMeasure::Recursive2Body,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2,
         "Recursive two-body decomposition measure with spectator/pair indices.")
    .def_static("DalitzPair", &PhaseSpaceMeasure::DalitzPair,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2,
         "Dalitz-variable measure over the chosen pair.")
    .def_static("HelicityAngles", &PhaseSpaceMeasure::HelicityAngles,
         arg("spectator") = 0, arg("pair_first") = 1, arg("pair_second") = 2,
         "Helicity-frame angle measure.")
    .def_static("MandelstamQ2", &PhaseSpaceMeasure::MandelstamQ2,
         "Momentum-transfer Q^2 measure with uniform azimuth integrated out.")
    .def_static("MandelstamQ2Phi", &PhaseSpaceMeasure::MandelstamQ2Phi,
         "Momentum-transfer Q^2 and explicit beam-axis azimuth measure.")
    .def_static("BjorkenXY", &PhaseSpaceMeasure::BjorkenXY,
         "Bjorken x,y measure with uniform azimuth integrated out.")
    .def_static("BjorkenXYPhi", &PhaseSpaceMeasure::BjorkenXYPhi,
         "Bjorken x,y and explicit beam-axis azimuth measure.")
    .def_static("FixedMassY", &PhaseSpaceMeasure::FixedMassY,
         "Fixed-mass y measure with uniform azimuth integrated out.")
    .def_static("FixedMassYPhi", &PhaseSpaceMeasure::FixedMassYPhi,
         "Fixed-mass y and explicit beam-axis azimuth measure.")
    .def_static("MandelstamQ2Y", &PhaseSpaceMeasure::MandelstamQ2Y,
         "Momentum-transfer Q^2 and y measure with uniform azimuth integrated out.")
    .def_static("MandelstamQ2YPhi", &PhaseSpaceMeasure::MandelstamQ2YPhi,
         "Momentum-transfer Q^2, y, and explicit beam-axis azimuth measure.")
    .def_static("Unspecified", &PhaseSpaceMeasure::Unspecified,
         "No measure declared.")
    ;

  class_<PhaseSpaceConvention>(m, "PhaseSpaceConvention")
    .def(init<>())
    .def(init([](PhaseSpaceTopology topology, PhaseSpaceMeasure measure) {
             return PhaseSpaceConvention{topology, measure};
         }),
         arg("topology"), arg("measure"))
    .def_readwrite("topology", &PhaseSpaceConvention::topology)
    .def_readwrite("measure", &PhaseSpaceConvention::measure)
    .def("__eq__", &PhaseSpaceConvention::operator==)
    .def("__ne__", &PhaseSpaceConvention::operator!=)
    ;

  m.def("PhaseSpaceTopologyName", &PhaseSpaceTopologyName);
  m.def("PhaseSpaceMeasureName", &PhaseSpaceMeasureName);
  m.def("PhaseSpaceDensityConvertible", &PhaseSpaceDensityConvertible,
        arg("topology"), arg("from_measure"), arg("to_measure"),
        "Return whether a density can be converted pointwise in the requested direction.");

  // Convert a sampling density between phase-space measures within a topology,
  // applying the analytic Jacobian (the same conversion the mixture uses).
  m.def("ConvertDensity",
        (double (*)(double, PhaseSpaceMeasure const &, PhaseSpaceMeasure const &,
                    PhaseSpaceTopology, siren::dataclasses::InteractionRecord const &))
            &siren::injection::ConvertDensity,
        arg("density"), arg("from_measure"), arg("to_measure"),
        arg("topology"), arg("record"));

  class_<PhaseSpaceChannel, std::shared_ptr<PhaseSpaceChannel>>(m, "PhaseSpaceChannel")
    .def("Sample", &PhaseSpaceChannel::Sample)
    .def("Density", &PhaseSpaceChannel::Density)
    .def("Name", &PhaseSpaceChannel::Name)
    .def("Topology", &PhaseSpaceChannel::Topology)
    .def("Measure", &PhaseSpaceChannel::Measure)
    ;

  class_<MultiChannelPhaseSpace, std::shared_ptr<MultiChannelPhaseSpace>> multi_channel_phase_space(m, "MultiChannelPhaseSpace",
      R"pbdoc(
      Weighted mixture of PhaseSpaceChannel objects sampled and evaluated as one
      combined density g(x) = sum_i alpha_i g_i(x). Each channel proposes candidate
      kinematics and every channel's density is evaluated at whatever point was
      drawn, so the combined density stays consistent regardless of which channel
      produced the sample.
      )pbdoc");

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
    .def_readwrite("channels", &MultiChannelPhaseSpace::channels,
         "The list of PhaseSpaceChannel objects making up the mixture.")
    .def_readwrite("weights", &MultiChannelPhaseSpace::weights,
         "Per-channel mixture weights alpha_i; call Normalize() after assigning unnormalized values.")
    .def("Normalize", &MultiChannelPhaseSpace::Normalize,
         "Rescale weights in place so they sum to one.")
    .def("Sample", &MultiChannelPhaseSpace::Sample,
         "Pick a channel by weight and draw kinematics from it into record.")
    .def("Density", &MultiChannelPhaseSpace::Density,
         "Evaluate the combined mixture density g(x) = sum_i alpha_i g_i(x) at record.")
    .def("DensityIn",
         overload_cast<
             std::shared_ptr<siren::detector::DetectorModel const>,
             siren::dataclasses::InteractionRecord const &,
             PhaseSpaceMeasure const &>(
             &MultiChannelPhaseSpace::DensityIn, const_),
         arg("detector_model"), arg("record"), arg("measure"),
         "Evaluate the combined mixture density in an explicitly requested measure.")
    .def("DensityBreakdown", &MultiChannelPhaseSpace::DensityBreakdown,
         arg("detector_model"), arg("record"),
         "Per-channel density contributions at record, for diagnosing which "
         "channel dominates the mixture at a given point.")
    .def("CommonTopology", &MultiChannelPhaseSpace::CommonTopology)
    .def("CommonMeasure", &MultiChannelPhaseSpace::CommonMeasure)
    .def("CommonConvention", &MultiChannelPhaseSpace::CommonConvention)
    .def("ValidateChannelsDetailed", &MultiChannelPhaseSpace::ValidateChannelsDetailed)
    .def("ValidateChannelDensities", &MultiChannelPhaseSpace::ValidateChannelDensities,
         arg("random"), arg("detector_model"), arg("template_record"),
         arg("samples_per_channel") = 100)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<MultiChannelPhaseSpace>),
        &(siren::serialization::pickle_load<MultiChannelPhaseSpace>)
    ))
    ;

  // A mixture used as one channel of another mixture.
  class_<NestedMixtureChannel, std::shared_ptr<NestedMixtureChannel>, PhaseSpaceChannel>(m, "NestedMixtureChannel")
    .def(init<std::shared_ptr<MultiChannelPhaseSpace>>(), arg("mixture"))
    .def_readwrite("mixture", &NestedMixtureChannel::mixture)
    .def_readwrite("label", &NestedMixtureChannel::label)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<NestedMixtureChannel>),
        &(siren::serialization::pickle_load<NestedMixtureChannel>)
    ))
    ;

  class_<Isotropic2BodyChannel, std::shared_ptr<Isotropic2BodyChannel>, PhaseSpaceChannel>(m, "Isotropic2BodyChannel",
      "Samples a two-body decay isotropically in the parent rest frame. "
      "Topology Decay2Body, measure SolidAngleRest.")
    .def(init<int>(), arg("daughter_index") = 0)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<Isotropic2BodyChannel>),
        &(siren::serialization::pickle_load<Isotropic2BodyChannel>)
    ))
    ;

  class_<PhysicalDecayChannel, std::shared_ptr<PhysicalDecayChannel>, PhaseSpaceChannel>(m, "PhysicalDecayChannel",
      "Samples a Decay's final state with the model's own sampler and reports "
      "its FinalStateProbability. The topology and measure are the model's "
      "declaration for the signature, or the convention passed explicitly.")
    .def(init<std::shared_ptr<siren::interactions::Decay>>())
    .def(init<std::shared_ptr<siren::interactions::Decay>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::Decay>,
              siren::dataclasses::InteractionSignature const &,
              PhaseSpaceConvention const &>(),
         arg("decay"), arg("signature"), arg("convention"))
    .def("GetDecay", &PhysicalDecayChannel::GetDecay)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<PhysicalDecayChannel>),
        &(siren::serialization::pickle_load<PhysicalDecayChannel>)
    ))
    ;

  class_<PhysicalCrossSectionChannel, std::shared_ptr<PhysicalCrossSectionChannel>, PhaseSpaceChannel>(m, "PhysicalCrossSectionChannel",
      "Samples a CrossSection's final state with the model's own sampler and "
      "reports its FinalStateProbability. The topology and measure are the "
      "model's declaration for the signature, or the convention passed explicitly.")
    .def(init<std::shared_ptr<siren::interactions::CrossSection>>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>,
              siren::dataclasses::InteractionSignature const &,
              PhaseSpaceConvention const &>(),
         arg("cross_section"), arg("signature"), arg("convention"))
    .def("GetCrossSection", &PhysicalCrossSectionChannel::GetCrossSection)
    .def(pybind11::pickle(
        &(siren::serialization::pickle_save<PhysicalCrossSectionChannel>),
        &(siren::serialization::pickle_load<PhysicalCrossSectionChannel>)
    ))
    ;

  // Two-body kinematics utilities

  m.def("TwoBodyRestMomentum", &TwoBodyRestMomentum);
  m.def("TwoBodyRestEnergy", &TwoBodyRestEnergy);
  m.def("Kallen", &Kallen);

  // Process

  // The keep_alive policies below tie python-defined distributions, cross
  // sections, and decays to the process, injector, or weighter consuming
  // them; a python-defined object held only by C++ shared_ptrs loses its
  // python half to garbage collection and virtual calls then fail.

  class_<Process, std::shared_ptr<Process>>(m, "Process")
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, cpp_function(&Process::SetInteractions, keep_alive<1, 2>()))
    ;

  class_<PhysicalProcess, std::shared_ptr<PhysicalProcess>, Process>(m, "PhysicalProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>(), keep_alive<1, 3>())
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, cpp_function(&Process::SetInteractions, keep_alive<1, 2>()))
    .def_property("distributions", &PhysicalProcess::GetPhysicalDistributions, cpp_function(&PhysicalProcess::SetPhysicalDistributions, keep_alive<1, 2>()))
    .def("SetPhaseSpace", &PhysicalProcess::SetPhaseSpace)
    .def("GetPhaseSpace", &PhysicalProcess::GetPhaseSpace)
    .def("GetPhaseSpaceMap", &PhysicalProcess::GetPhaseSpaceMap)
    .def("HasPhaseSpace", overload_cast<siren::dataclasses::InteractionSignature const &>(&PhysicalProcess::HasPhaseSpace, const_))
    .def("HasAnyPhaseSpace", &PhysicalProcess::HasAnyPhaseSpace)
    .def("SetWeightingMode", &PhysicalProcess::SetWeightingMode)
    .def("GetWeightingMode", &PhysicalProcess::GetWeightingMode)
    .def_property("weighting_mode", &PhysicalProcess::GetWeightingMode, &PhysicalProcess::SetWeightingMode)
    ;

  class_<PrimaryInjectionProcess, std::shared_ptr<PrimaryInjectionProcess>, PhysicalProcess>(m, "PrimaryInjectionProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>(), keep_alive<1, 3>())
    .def_property("primary_type", &Process::GetPrimaryType, &Process::SetPrimaryType)
    .def_property("interactions", &Process::GetInteractions, cpp_function(&Process::SetInteractions, keep_alive<1, 2>()))
    .def_property("distributions", &PrimaryInjectionProcess::GetPrimaryInjectionDistributions, cpp_function(&PrimaryInjectionProcess::SetPrimaryInjectionDistributions, keep_alive<1, 2>()))
    ;

  class_<SecondaryInjectionProcess, std::shared_ptr<SecondaryInjectionProcess>, PhysicalProcess>(m, "SecondaryInjectionProcess")
    .def(init<>())
    .def(init<siren::dataclasses::ParticleType, std::shared_ptr<siren::interactions::InteractionCollection>>(), keep_alive<1, 3>())
    .def_property("secondary_type", &SecondaryInjectionProcess::GetSecondaryType, &SecondaryInjectionProcess::SetSecondaryType)
    .def_property("interactions", &Process::GetInteractions, cpp_function(&Process::SetInteractions, keep_alive<1, 2>()))
    .def_property("distributions", &SecondaryInjectionProcess::GetSecondaryInjectionDistributions, cpp_function(&SecondaryInjectionProcess::SetSecondaryInjectionDistributions, keep_alive<1, 2>()))
    ;

  // Injection

  class_<FailureLedger>(m, "FailureLedger")
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
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::shared_ptr<siren::utilities::SIREN_random>>(), keep_alive<1, 4>())
    .def(init<unsigned int, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PrimaryInjectionProcess>, std::vector<std::shared_ptr<SecondaryInjectionProcess>>, std::shared_ptr<siren::utilities::SIREN_random>>(), keep_alive<1, 4>(), keep_alive<1, 5>())
    .def("SetStoppingCondition",&Injector::SetStoppingCondition)
    .def("GetStoppingCondition",&Injector::GetStoppingCondition)
    .def("SetPrimaryProcess",&Injector::SetPrimaryProcess, keep_alive<1, 2>())
    .def("AddSecondaryProcess",&Injector::AddSecondaryProcess, keep_alive<1, 2>())
    .def("GetPrimaryProcess",&Injector::GetPrimaryProcess)
    .def("GetSecondaryProcesses",&Injector::GetSecondaryProcesses)
    .def("GetSecondaryProcessMap",&Injector::GetSecondaryProcessMap)
    .def("NewRecord",&Injector::NewRecord)
    .def("SetRandom",&Injector::SetRandom)
    .def("GetRandom",&Injector::GetRandom)
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
    .def("UnregisteredSecondaryCount",&Injector::UnregisteredSecondaryCount)
    .def("GetLastFailureMessage",&Injector::GetLastFailureMessage)
    .def("GetLastFailedTree",&Injector::GetLastFailedTree, pybind11::return_value_policy::reference_internal)
    .def("GetFailureLedger",&Injector::GetFailureLedger, pybind11::return_value_policy::reference_internal)
    .def("ResetInjectedEvents",overload_cast<unsigned int>(&Injector::ResetInjectedEvents))
    .def("ResetInjectedEvents",overload_cast<>(&Injector::ResetInjectedEvents))
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


  // Weighter classes

  class_<VertexWeightFactors>(m, "VertexWeightFactors")
    .def_readonly("injector_index", &VertexWeightFactors::injector_index)
    .def_readonly("depth", &VertexWeightFactors::depth)
    .def_readonly("vertex_pdg", &VertexWeightFactors::vertex_pdg)
    .def_readonly("generation", &VertexWeightFactors::generation)
    .def_readonly("physical", &VertexWeightFactors::physical)
    .def_readonly("interaction_prob", &VertexWeightFactors::interaction_prob)
    .def_readonly("position_prob", &VertexWeightFactors::position_prob)
    .def_readonly("channel_density_topology", &VertexWeightFactors::channel_density_topology)
    .def_readonly("channel_density_measure", &VertexWeightFactors::channel_density_measure)
    .def_readonly("channel_densities", &VertexWeightFactors::channel_densities)
    .def_readonly("cancelled", &VertexWeightFactors::cancelled)
    .def_readonly("flags", &VertexWeightFactors::flags);

  class_<EventWeightBreakdown>(m, "EventWeightBreakdown")
    .def_readonly("total", &EventWeightBreakdown::total)
    .def_readonly("vertices", &EventWeightBreakdown::vertices);

  class_<PrimaryProcessWeighter, std::shared_ptr<PrimaryProcessWeighter>>(m, "PrimaryProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<PrimaryInjectionProcess>, std::shared_ptr<siren::detector::DetectorModel>>())
    .def("InteractionProbability",&PrimaryProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&PrimaryProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",
         overload_cast<
             std::tuple<siren::math::Vector3D, siren::math::Vector3D> const &,
             siren::dataclasses::InteractionRecord const &>(
             &PrimaryProcessWeighter::PhysicalProbability, const_))
    .def("GenerationProbability",
         overload_cast<siren::dataclasses::InteractionTreeDatum const &>(
             &PrimaryProcessWeighter::GenerationProbability, const_))
    .def("EventWeight",&PrimaryProcessWeighter::EventWeight)
    ;

  class_<SecondaryProcessWeighter, std::shared_ptr<SecondaryProcessWeighter>>(m, "SecondaryProcessWeighter")
    .def(init<std::shared_ptr<PhysicalProcess>, std::shared_ptr<SecondaryInjectionProcess>, std::shared_ptr<siren::detector::DetectorModel>>())
    .def("InteractionProbability",&SecondaryProcessWeighter::InteractionProbability)
    .def("NormalizedPositionProbability",&SecondaryProcessWeighter::NormalizedPositionProbability)
    .def("PhysicalProbability",
         overload_cast<
             std::tuple<siren::math::Vector3D, siren::math::Vector3D> const &,
             siren::dataclasses::InteractionRecord const &>(
             &SecondaryProcessWeighter::PhysicalProbability, const_))
    .def("GenerationProbability",
         overload_cast<siren::dataclasses::InteractionTreeDatum const &>(
             &SecondaryProcessWeighter::GenerationProbability, const_))
    .def("EventWeight",&SecondaryProcessWeighter::EventWeight)
    ;

  class_<Weighter, std::shared_ptr<Weighter>>(m, "Weighter")
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>, std::vector<std::shared_ptr<PhysicalProcess>>>(), keep_alive<1, 2>(), keep_alive<1, 4>(), keep_alive<1, 5>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::shared_ptr<siren::detector::DetectorModel>, std::shared_ptr<PhysicalProcess>>(), keep_alive<1, 2>(), keep_alive<1, 4>())
    .def(init<std::vector<std::shared_ptr<Injector>>, std::string>(), keep_alive<1, 2>())
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
