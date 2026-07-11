
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
#include "../../public/SIREN/injection/GeometryVolume.h"
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

  // New topology/measure enums
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
    .value("KinematicallyForbidden", siren::utilities::FailureReason::KinematicallyForbidden)
    .value("UnregisteredSecondaryType", siren::utilities::FailureReason::UnregisteredSecondaryType)
    .value("PrimaryVertexFailure", siren::utilities::FailureReason::PrimaryVertexFailure)
    .value("TopLevelCatch", siren::utilities::FailureReason::TopLevelCatch);

  enum_<PhaseSpaceMeasure::Type>(m, "PhaseSpaceMeasureType")
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
      combined density g(x) = sum_i alpha_i g_i(x). Mixes physical and detector-
      directed channels for importance sampling: each channel proposes candidate
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
         "Per-channel mixture weights alpha_i; need not be pre-normalized.")
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
    .def("Accumulate", &MultiChannelPhaseSpace::Accumulate,
         arg("detector_model"), arg("record"), arg("weight"),
         arg("discount_fallback") = true, arg("recurse") = true,
         "Feed one weighted sample into the per-channel statistics used by "
         "UpdateWeights to retune the mixture weights.")
    .def("UpdateWeights", &MultiChannelPhaseSpace::UpdateWeights,
         arg("update_rule"), arg("damping"), arg("min_weight"),
         arg("recurse") = true, arg("failure_mode") = "throughput",
         "Retune channel weights from accumulated statistics using update_rule, "
         "damping toward the previous weights and floored at min_weight.")
    .def("ResetAccumulators", &MultiChannelPhaseSpace::ResetAccumulators,
         arg("recurse") = true,
         "Clear the statistics accumulated by Accumulate/AccumulateSelection.")
    .def("AccumulateSelection", &MultiChannelPhaseSpace::AccumulateSelection,
         arg("detector_model"), arg("record"), arg("failed"),
         "Feed one selection outcome (pass/fail) into the per-channel "
         "acceptance statistics used by UpdateWeights.")
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

  class_<Isotropic2BodyChannel, std::shared_ptr<Isotropic2BodyChannel>, PhaseSpaceChannel>(m, "Isotropic2BodyChannel",
      "Samples a two-body decay isotropically in the parent rest frame, with no "
      "detector direction bias. Topology Decay2Body, Measure SolidAngleRest. Use as "
      "the physical/fallback channel for a 2-body decay with no directing bias.")
    .def(init<int>(), arg("daughter_index") = 0)
    ;

  enum_<DetectorDirected2BodyChannel::Mode>(m, "DirectedMode")
    .value("Cone", DetectorDirected2BodyChannel::Mode::Cone,
           "Bias into a fixed cone around the detector direction; simple but "
           "leaves solid angle outside the cone unsampled.")
    .value("Volume", DetectorDirected2BodyChannel::Mode::Volume,
           "Bias toward the target volume's actual angular extent; adapts to "
           "detector geometry.");

  class_<DetectorDirected2BodyChannel, std::shared_ptr<DetectorDirected2BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected2BodyChannel",
      "Directs one daughter of a two-body decay toward a target volume. Topology "
      "Decay2Body, Measure SolidAngleRest. Use for 2-body decays where one "
      "daughter is the signal to point at the detector.")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int, DetectorDirected2BodyChannel::Mode, double>(),
         arg("target"), arg("daughter_index") = 0,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("volume") = -1.0,
         "target: geometry to bias toward.\n"
         "daughter_index: which of the two final-state particles is directed.\n"
         "mode: DirectedMode.Cone or DirectedMode.Volume.\n"
         "volume: fixed solid angle for Cone mode; -1 to derive it from geometry.")
    .def("SetVolume", &DetectorDirected2BodyChannel::SetVolume)
    .def("DirectingActive", &DetectorDirected2BodyChannel::DirectingActive, arg("record"))
    ;

  class_<DetectorDirectedAngularSectorChannel, std::shared_ptr<DetectorDirectedAngularSectorChannel>, PhaseSpaceChannel>(m, "DetectorDirectedAngularSectorChannel",
      "Directs a two-body daughter into one angular sector (u, phi bounds) of the "
      "target, for disjoint tiling of the target's angular extent across channels. "
      "Topology Decay2Body, Measure SolidAngleRest.")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, double, double, double, double, int>(),
         arg("target"), arg("u_lo"), arg("u_hi"), arg("phi_lo"), arg("phi_hi"),
         arg("daughter_index") = 0,
         "target: geometry the sector is defined against.\n"
         "u_lo, u_hi: bounds in cos(theta) (u = cos theta) for this sector.\n"
         "phi_lo, phi_hi: azimuthal bounds for this sector.\n"
         "daughter_index: which of the two final-state particles is directed.")
    .def("DirectingActive", &DetectorDirectedAngularSectorChannel::DirectingActive, arg("record"))
    ;

  // 1-D importance maps: one object provides BOTH the draw (Forward) and
  // its own normalized density (Density), so a model/channel routing both
  // through a shared instance cannot let sampling and density drift apart
  // (Contract C1).  Models consume these; they do not subclass in Python,
  // so no trampoline is needed.
  class_<Mapping1D, std::shared_ptr<Mapping1D>>(m, "Mapping1D",
      "Base interface for a 1-D importance map: one object provides both the draw "
      "(Forward, from a uniform variate) and its own normalized density (Density) "
      "over the same variable, so a model/channel routing both through a shared "
      "instance cannot let sampling and density drift apart.")
    .def("Forward", &Mapping1D::Forward, arg("r"))
    .def("Inverse", &Mapping1D::Inverse, arg("x"))
    .def("Density", &Mapping1D::Density, arg("x"))
    .def("Accumulate", &Mapping1D::Accumulate, arg("x"), arg("weight"))
    .def("Refine", &Mapping1D::Refine);

  class_<BreitWignerMapping, std::shared_ptr<BreitWignerMapping>, Mapping1D>(m, "BreitWignerMapping",
      "Breit-Wigner-shaped density in s over [s_min, s_max]; peaks at s = mass^2 "
      "with width set by width.")
    .def(init<double, double, double, double>(),
         arg("mass"), arg("width"), arg("s_min"), arg("s_max"));

  class_<PowerLawMapping, std::shared_ptr<PowerLawMapping>, Mapping1D>(m, "PowerLawMapping",
      "Power-law-shaped density (index nu, offset m2) over [s_min, s_max]; use for "
      "a heavy-tailed invariant-mass or Q^2 variable with no resonance.")
    .def(init<double, double, double, double>(),
         arg("nu"), arg("m2"), arg("s_min"), arg("s_max"));

  class_<TabulatedMapping, std::shared_ptr<TabulatedMapping>, Mapping1D>(m, "TabulatedMapping",
      "Density defined by a user-supplied cumulative table (s_nodes, cdf_nodes) "
      "over [s_min, s_max]; use when the target density has no closed form.")
    .def(init<std::vector<double>, std::vector<double>, double, double>(),
         arg("s_nodes"), arg("cdf_nodes"), arg("s_min"), arg("s_max"));

  class_<PropagatorMapping, std::shared_ptr<PropagatorMapping>, Mapping1D>(m, "PropagatorMapping",
      "1/(x^2 + m2)-shaped density for propagator-peaked variables like off-shell "
      "Q^2 or invariant mass.")
    .def(init<double, double, double>(),
         arg("m2"), arg("x_min"), arg("x_max"));

  class_<UniformMapping, std::shared_ptr<UniformMapping>, Mapping1D>(m, "UniformMapping",
      "Flat density over [s_min, s_max].")
    .def(init<double, double>(),
         arg("s_min"), arg("s_max"));

  class_<LogMapping, std::shared_ptr<LogMapping>, Mapping1D>(m, "LogMapping",
      "Log-uniform density over [x_min, x_max]; use for a variable spanning "
      "multiple orders of magnitude.")
    .def(init<double, double>(),
         arg("x_min"), arg("x_max"));

  class_<ExponentialMapping, std::shared_ptr<ExponentialMapping>, Mapping1D>(m, "ExponentialMapping",
      "Exponential density with mean tau over [x_min, x_max].")
    .def(init<double, double, double>(),
         arg("tau"), arg("x_min"), arg("x_max"));

  class_<GaussianMapping, std::shared_ptr<GaussianMapping>, Mapping1D>(m, "GaussianMapping",
      "Gaussian density (mean mu, width sigma) over [x_min, x_max].")
    .def(init<double, double, double, double>(),
         arg("mu"), arg("sigma"), arg("x_min"), arg("x_max"));

  class_<AdaptiveMapping, std::shared_ptr<AdaptiveMapping>, Mapping1D>(m, "AdaptiveMapping",
      "Piecewise-constant density over [x_min, x_max] in n_bins bins, refined "
      "toward accumulated samples via Refine; use when no analytic shape fits "
      "and the density should self-tune during generation.")
    .def(init<double, double, int, double, double>(),
         arg("x_min"), arg("x_max"), arg("n_bins") = 32,
         arg("damping") = 0.5, arg("floor_frac") = 1e-3)
    .def_readonly("p", &AdaptiveMapping::p)
    .def_readonly("n_bins", &AdaptiveMapping::n_bins);

  enum_<DetectorDirected3BodyChannel::InvariantMassMode>(m, "InvariantMassMode")
    .value("Uniform", DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
           "Flat invariant-mass sampling.")
    .value("BreitWigner", DetectorDirected3BodyChannel::InvariantMassMode::BreitWigner,
           "Breit-Wigner resonance shape (needs resonance_mass, resonance_width).")
    .value("PowerLaw", DetectorDirected3BodyChannel::InvariantMassMode::PowerLaw,
           "Power-law tail (needs power_law_nu, power_law_offset).")
    .value("Tabulated", DetectorDirected3BodyChannel::InvariantMassMode::Tabulated,
           "User-supplied cumulative table (mass_cdf_nodes/values).");

  enum_<DetectorDirected3BodyChannel::Factorization>(m, "ThreeBodyMode")
    .value("Direct", DetectorDirected3BodyChannel::Factorization::Direct,
           "Bias one daughter directly toward the detector; the invariant mass "
           "pairs the other two. Best when the biased daughter is the signal "
           "and the spectator pair has no resonance.")
    .value("Recursive", DetectorDirected3BodyChannel::Factorization::Recursive,
           "Bias a two-body sub-decay (a resonance pair) toward the detector, "
           "then decay the pair internally. Best when the pair has resonance "
           "structure (e.g. an off-shell mediator sampled with Breit-Wigner).");

  class_<DetectorDirected3BodyChannel, std::shared_ptr<DetectorDirected3BodyChannel>, PhaseSpaceChannel>(m, "DetectorDirected3BodyChannel",
      "Directs a three-body decay toward a target volume, factorized per "
      "ThreeBodyMode. Topology Decay3Body, Measure Recursive2Body. Use for "
      "3-body BSM decays where a signal daughter or resonance pair should "
      "point at the detector.")
    // Keyword ctor: factorization selects Direct vs Recursive explicitly, so
    // pybind disambiguates on the leading Factorization-typed argument rather
    // than on argument count. Constructs via the underlying Direct/Recursive
    // C++ ctors, so objects built here are byte-identical to ones built there.
    .def(init(
        [](DetectorDirected3BodyChannel::Factorization factorization,
           std::shared_ptr<siren::geometry::Geometry const> target,
           int directed_index,
           int spectator_index,
           int pair_first_index,
           int pair_second_index,
           DetectorDirected3BodyChannel::InvariantMassMode mass_mode,
           double resonance_mass,
           double resonance_width,
           double power_law_nu,
           double power_law_offset,
           DetectorDirected2BodyChannel::Mode mode,
           PhaseSpaceTopology topology,
           std::vector<double> mass_cdf_nodes,
           std::vector<double> mass_cdf_values,
           double volume) {
          if (factorization == DetectorDirected3BodyChannel::Factorization::Direct) {
            return std::make_shared<DetectorDirected3BodyChannel>(
                target, directed_index,
                mass_mode, resonance_mass, resonance_width,
                power_law_nu, power_law_offset,
                mode, topology, mass_cdf_nodes, mass_cdf_values, volume);
          }
          if (spectator_index < 0 || pair_first_index < 0 || pair_second_index < 0) {
            throw siren::utilities::ConfigurationError(
                "factorization=ThreeBodyMode.Recursive requires spectator_index, "
                "pair_first_index, and pair_second_index");
          }
          int directed_pair_index =
              (directed_index == pair_first_index || directed_index == pair_second_index)
                  ? directed_index : pair_first_index;
          return std::make_shared<DetectorDirected3BodyChannel>(
              target, spectator_index, pair_first_index, pair_second_index,
              directed_pair_index,
              mass_mode, resonance_mass, resonance_width,
              power_law_nu, power_law_offset,
              mode, topology, mass_cdf_nodes, mass_cdf_values, volume);
        }),
        arg("factorization"), arg("target"),
        arg("directed_index") = 0,
        arg("spectator_index") = -1,
        arg("pair_first_index") = -1,
        arg("pair_second_index") = -1,
        arg("mass_mode") = DetectorDirected3BodyChannel::InvariantMassMode::Uniform,
        arg("resonance_mass") = 0.0, arg("resonance_width") = 0.0,
        arg("power_law_nu") = 0.8, arg("power_law_offset") = 0.0,
        arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
        arg("topology") = PhaseSpaceTopology::Decay3Body,
        arg("mass_cdf_nodes") = std::vector<double>{},
        arg("mass_cdf_values") = std::vector<double>{},
        arg("volume") = -1.0,
        "Keyword constructor selecting the factorization explicitly via "
        "factorization=ThreeBodyMode.Direct or ThreeBodyMode.Recursive.\n"
        "factorization: ThreeBodyMode.Direct or ThreeBodyMode.Recursive.\n"
        "target: geometry to bias toward.\n"
        "directed_index: (Direct) which daughter is biased toward target; "
        "(Recursive) which member of the pair is biased, defaulting to "
        "pair_first_index if not one of the pair indices.\n"
        "spectator_index, pair_first_index, pair_second_index: (Recursive only) "
        "indices of the spectator daughter and the resonance-pair daughters.\n"
        "mass_mode: InvariantMassMode governing the paired/spectator invariant mass.\n"
        "resonance_mass, resonance_width: used when mass_mode is BreitWigner.\n"
        "power_law_nu, power_law_offset: used when mass_mode is PowerLaw.\n"
        "mode: DirectedMode.Cone or DirectedMode.Volume for the directed sub-step.\n"
        "topology: PhaseSpaceTopology tag carried by the channel.\n"
        "mass_cdf_nodes, mass_cdf_values: used when mass_mode is Tabulated.\n"
        "volume: exact target volume override (skips the viability estimate).")
    .def("SetVolume", &DetectorDirected3BodyChannel::SetVolume)
    .def("DirectingActive", &DetectorDirected3BodyChannel::DirectingActive, arg("record"))
    ;

  enum_<DetectorDirectedScatteringChannel::Variable>(m, "ScatteringVariable")
    .value("Q2", DetectorDirectedScatteringChannel::Variable::Q2,
           "Sample in momentum transfer Q^2.")
    .value("BjorkenY", DetectorDirectedScatteringChannel::Variable::BjorkenY,
           "Sample in inelasticity y.")
    .value("RecoilY", DetectorDirectedScatteringChannel::Variable::RecoilY,
           "Sample in recoil-energy fraction.");

  enum_<DetectorDirectedScatteringChannel::Q2Mode>(m, "ScatteringQ2Mode")
    .value("Geometry", DetectorDirectedScatteringChannel::Q2Mode::Geometry,
           "Q^2 range set by detector geometry (preserves current behavior).")
    .value("Propagator", DetectorDirectedScatteringChannel::Q2Mode::Propagator,
           "Q^2 sampled from the propagator peak (lower variance when the "
           "density is propagator-peaked).")
    .value("Tabulated", DetectorDirectedScatteringChannel::Q2Mode::Tabulated,
           "Q^2 from a user CDF table (q2_cdf_nodes/values).");

  class_<DetectorDirectedScatteringChannel, std::shared_ptr<DetectorDirectedScatteringChannel>, PhaseSpaceChannel>(m, "DetectorDirectedScatteringChannel",
      "Directs the recoil of a 2->2 scatter toward a target volume. Topology "
      "Scatter2to2, with an explicit-azimuth Q2 or y measure. Use for upscattering where the "
      "outgoing state should point at the detector.")
    .def(init<std::shared_ptr<siren::geometry::Geometry const>, int,
              DetectorDirectedScatteringChannel::Variable,
              DetectorDirected2BodyChannel::Mode,
              DetectorDirectedScatteringChannel::Q2Mode, double,
              std::vector<double>, std::vector<double>, double>(),
         arg("target"), arg("directed_index") = 0,
         arg("variable") = DetectorDirectedScatteringChannel::Variable::Q2,
         arg("mode") = DetectorDirected2BodyChannel::Mode::Volume,
         arg("q2_mode") = DetectorDirectedScatteringChannel::Q2Mode::Geometry,
         arg("mediator_mass") = 0.0,
         arg("q2_cdf_nodes") = std::vector<double>{},
         arg("q2_cdf_values") = std::vector<double>{},
         arg("volume") = -1.0,
         "target: geometry to bias the recoil toward.\n"
         "directed_index: which outgoing particle is directed.\n"
         "variable: ScatteringVariable sampled (Q2, BjorkenY, or RecoilY).\n"
         "mode: DirectedMode.Cone or DirectedMode.Volume for the directed sub-step.\n"
         "q2_mode: ScatteringQ2Mode governing how the Q^2 range/density is chosen.\n"
         "mediator_mass: propagator mass used when q2_mode is Propagator.\n"
         "q2_cdf_nodes, q2_cdf_values: used when q2_mode is Tabulated.\n"
         "volume: exact target volume override (skips the viability estimate).")
    .def("SetVolume", &DetectorDirectedScatteringChannel::SetVolume)
    ;

  // Physical channel adapters

  class_<PhysicalDecayChannel, std::shared_ptr<PhysicalDecayChannel>, PhaseSpaceChannel>(m, "PhysicalDecayChannel",
      "Samples the unbiased physical final state of a Decay (no detector "
      "direction). Serves as the fallback channel in a mixture so events that "
      "miss the target are still represented. Topology/Measure follow the "
      "underlying interaction.")
    .def(init<std::shared_ptr<siren::interactions::Decay>>())
    .def(init<std::shared_ptr<siren::interactions::Decay>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::Decay>,
              siren::dataclasses::InteractionSignature const &,
              PhaseSpaceConvention const &>(),
         arg("decay"), arg("signature"), arg("convention"))
    .def("GetDecay", &PhysicalDecayChannel::GetDecay)
    ;

  class_<PhysicalCrossSectionChannel, std::shared_ptr<PhysicalCrossSectionChannel>, PhaseSpaceChannel>(m, "PhysicalCrossSectionChannel",
      "Samples the unbiased physical final state of a CrossSection (no detector "
      "direction). Serves as the fallback channel in a mixture so events that "
      "miss the target are still represented. Topology/Measure follow the "
      "underlying interaction.")
    .def(init<std::shared_ptr<siren::interactions::CrossSection>>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>,
              siren::dataclasses::InteractionSignature const &>())
    .def(init<std::shared_ptr<siren::interactions::CrossSection>,
              siren::dataclasses::InteractionSignature const &,
              PhaseSpaceConvention const &>(),
         arg("cross_section"), arg("signature"), arg("convention"))
    .def("GetCrossSection", &PhysicalCrossSectionChannel::GetCrossSection)
    ;

  // Two-body kinematics utilities

  m.def("TwoBodyRestMomentum", &TwoBodyRestMomentum);
  m.def("TwoBodyRestEnergy", &TwoBodyRestEnergy);
  m.def("Kallen", &Kallen);

  // Exact analytic volume for the geometries the directed channels accept
  // in Volume mode; raises for composites, which need a caller-supplied
  // volume.
  m.def("geometry_volume", [](
      std::shared_ptr<siren::geometry::Geometry const> const & geometry) {
    if (!geometry) throw std::invalid_argument("geometry_volume requires a geometry");
    double volume = ExactGeometryVolume(*geometry);
    if (!std::isfinite(volume)) {
      throw std::invalid_argument(
          "exact volume is unavailable for this geometry; provide volume_fn");
    }
    return volume;
  });

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
    .def("GetLastFailureReason",&Injector::GetLastFailureReason)
    .def("GetLastFailedTree",&Injector::GetLastFailedTree, pybind11::return_value_policy::reference_internal)
    .def("ResetInjectedEvents",&Injector::ResetInjectedEvents)
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
