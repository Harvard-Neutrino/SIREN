
#include <set>
#include <limits>
#include <vector>
#include <string>

#include "../../public/SIREN/distributions/Distributions.h"
#include "../../public/SIREN/distributions/pyWeightableDistribution.h"
#include "../../public/SIREN/distributions/pyPrimaryInjectionDistribution.h"
#include "../../public/SIREN/distributions/pySecondaryInjectionDistribution.h"
#include "../../public/SIREN/distributions/primary/PrimaryExternalDistribution.h"
#include "../../public/SIREN/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "../../public/SIREN/distributions/primary/direction/pyPrimaryDirectionDistribution.h"
#include "../../public/SIREN/distributions/primary/direction/Cone.h"
#include "../../public/SIREN/distributions/primary/direction/FixedDirection.h"
#include "../../public/SIREN/distributions/primary/direction/IsotropicDirection.h"
#include "../../public/SIREN/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "../../public/SIREN/distributions/primary/energy/pyPrimaryEnergyDistribution.h"
#include "../../public/SIREN/distributions/primary/energy/Monoenergetic.h"
#include "../../public/SIREN/distributions/primary/energy/PowerLaw.h"
#include "../../public/SIREN/distributions/primary/energy/PiDARNuEDistribution.h"
#include "../../public/SIREN/distributions/primary/energy/TabulatedFluxDistribution.h"
#include "../../public/SIREN/distributions/primary/energy_direction/PrimaryEnergyDirectionDistribution.h"
#include "../../public/SIREN/distributions/primary/energy_direction/pyPrimaryEnergyDirectionDistribution.h"
#include "../../public/SIREN/distributions/primary/energy_direction/Tabulated2DFluxDistribution.h"
#include "../../public/SIREN/distributions/primary/helicity/PrimaryNeutrinoHelicityDistribution.h"
#include "../../public/SIREN/distributions/primary/mass/PrimaryMass.h"
#include "../../public/SIREN/distributions/primary/vertex/VertexPositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/pyVertexPositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/SphereVolumePositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/FixedTargetPositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/DecayRangeFunction.h"
#include "../../public/SIREN/distributions/primary/vertex/DecayRangePositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/DepthFunction.h"
#include "../../public/SIREN/distributions/primary/vertex/LeptonDepthFunction.h"
#include "../../public/SIREN/distributions/primary/vertex/PointSourcePositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/RangeFunction.h"
#include "../../public/SIREN/distributions/primary/vertex/RangePositionDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/PrimaryPhysicalVertexDistribution.h"
#include "../../public/SIREN/distributions/primary/vertex/PrimaryBoundedVertexDistribution.h"
#include "../../public/SIREN/distributions/primary/area/PrimaryAreaDistribution.h"
#include "../../public/SIREN/distributions/primary/area/pyPrimaryAreaDistribution.h"
#include "../../public/SIREN/distributions/primary/area/FixedTargetAreaDistribution.h"
#include "../../public/SIREN/distributions/secondary/vertex/pySecondaryVertexPositionDistribution.h"
#include "../../public/SIREN/distributions/secondary/vertex/SecondaryPhysicalVertexDistribution.h"
#include "../../public/SIREN/distributions/secondary/vertex/SecondaryBoundedVertexDistribution.h"
#include "../../public/SIREN/distributions/secondary/vertex/SecondaryDecayRangePositionDistribution.h"

#include "../../../utilities/public/SIREN/utilities/Random.h"
#include "../../../detector/public/SIREN/detector/DetectorModel.h"
#include "../../../interactions/public/SIREN/interactions/InteractionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(distributions,m) {
  using namespace siren::distributions;

  enum_<DistributionVariable>(m, "DistributionVariable")
    .value("PrimaryMass", DistributionVariable::PrimaryMass)
    .value("PrimaryEnergy", DistributionVariable::PrimaryEnergy)
    .value("PrimaryDirection", DistributionVariable::PrimaryDirection)
    .value("PrimaryHelicity", DistributionVariable::PrimaryHelicity)
    .value("PrimaryArea", DistributionVariable::PrimaryArea)
    .value("InitialPosition", DistributionVariable::InitialPosition)
    .value("InteractionVertex", DistributionVariable::InteractionVertex)
    .value("InteractionParameters", DistributionVariable::InteractionParameters);

  class_<PhysicallyNormalizedDistribution, std::shared_ptr<PhysicallyNormalizedDistribution>>(m, "PhysicallyNormalizedDistribution")
    .def(init<>())
    .def(init<double>())
    .def_property("normalization",&PhysicallyNormalizedDistribution::GetNormalization,&PhysicallyNormalizedDistribution::SetNormalization)
    .def("IsNormalizationSet",&PhysicallyNormalizedDistribution::IsNormalizationSet);

  class_<WeightableDistribution, std::shared_ptr<WeightableDistribution>, pyWeightableDistribution>(m, "WeightableDistribution",
    "Base class for weight-only distributions. Subclass in python and\n"
    "implement GenerationProbability(detector_model, interactions, record);\n"
    "optionally override Name, DensityVariables, equal, and less.")
    .def(init<>())
    .def("__eq__", [](const WeightableDistribution &self, const WeightableDistribution &other){ return self == other; })
    .def("GenerationProbability",&WeightableDistribution::GenerationProbability)
    .def("PhysicalDensity",&WeightableDistribution::PhysicalDensity,
         "Density contributed on the physical side of the weight ratio;\n"
         "defaults to GenerationProbability. Differs only for distributions\n"
         "whose physical role is not their sampling role (an external table\n"
         "with importance weights).")
    .def("PhysicalDensityDiffers",&WeightableDistribution::PhysicalDensityDiffers)
    .def("DensityVariables",&WeightableDistribution::DensityVariables)
    .def("Name",&WeightableDistribution::Name)
    .def("AreEquivalent",&WeightableDistribution::AreEquivalent)
    TrampolinePickleMethods(pyWeightableDistribution);

  class_<NormalizationConstant, std::shared_ptr<NormalizationConstant>, WeightableDistribution, PhysicallyNormalizedDistribution>(m, "NormalizationConstant")
    .def(init<double>())
    .def("GenerationProbability",&NormalizationConstant::GenerationProbability)
    .def("Name",&NormalizationConstant::Name);

  class_<PrimaryInjectionDistribution, std::shared_ptr<PrimaryInjectionDistribution>, pyPrimaryInjectionDistribution, WeightableDistribution>(m, "PrimaryInjectionDistribution",
    "Base class for primary injection distributions. Subclass in python and\n"
    "implement Sample(rand, detector_model, interactions, record) and\n"
    "GenerationProbability(detector_model, interactions, record).")
    .def(init<>())
    .def("Sample",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::PrimaryDistributionRecord &>(&PrimaryInjectionDistribution::Sample, const_))
    .def("SetVariables",&PrimaryInjectionDistribution::SetVariables)
    .def("RequiredVariables",&PrimaryInjectionDistribution::RequiredVariables)
    .def("clone",&PrimaryInjectionDistribution::clone)
    TrampolinePickleMethods(pyPrimaryInjectionDistribution);

  // NOTE: PrimaryExternalDistribution is defined after VertexPositionDistribution
  // because it inherits from VertexPositionDistribution.

  // Direction distributions

  class_<PrimaryDirectionDistribution, std::shared_ptr<PrimaryDirectionDistribution>, pyPrimaryDirectionDistribution, PrimaryInjectionDistribution>(m, "PrimaryDirectionDistribution",
    "Base class for primary direction distributions. Subclass in python and\n"
    "implement SampleDirection(rand, detector_model, interactions, record)\n"
    "returning a Vector3D, plus GenerationProbability.")
    .def(init<>())
    .def("Sample",&PrimaryDirectionDistribution::Sample)
    .def("DensityVariables",&PrimaryDirectionDistribution::DensityVariables)
    .def("GenerationProbability",&PrimaryDirectionDistribution::GenerationProbability)
    TrampolinePickleMethods(pyPrimaryDirectionDistribution);

  class_<Cone, std::shared_ptr<Cone>, PrimaryDirectionDistribution>(m, "Cone")
    .def(init<siren::math::Vector3D, double>())
    .def(init([](std::array<double, 3> const & dir, double opening_angle) {
        return Cone(siren::math::Vector3D(dir), opening_angle);
    }), arg("direction"), arg("opening_angle"));

  class_<IsotropicDirection, std::shared_ptr<IsotropicDirection>, PrimaryDirectionDistribution>(m, "IsotropicDirection")
    .def(init<>());

  class_<FixedDirection, std::shared_ptr<FixedDirection>, PrimaryDirectionDistribution>(m, "FixedDirection")
    .def(init<siren::math::Vector3D>())
    .def(init([](std::array<double, 3> const & dir) {
        return FixedDirection(siren::math::Vector3D(dir));
    }), arg("direction"));

  // Energy distributions

  class_<PrimaryEnergyDistribution, std::shared_ptr<PrimaryEnergyDistribution>, pyPrimaryEnergyDistribution, PrimaryInjectionDistribution, PhysicallyNormalizedDistribution>(m, "PrimaryEnergyDistribution",
    "Base class for primary energy distributions. Subclass in python and\n"
    "implement SampleEnergy(rand, detector_model, interactions, record)\n"
    "returning the sampled energy, plus GenerationProbability.")
    .def(init<>())
    .def("Sample",&PrimaryEnergyDistribution::Sample)
    TrampolinePickleMethods(pyPrimaryEnergyDistribution);

  class_<Monoenergetic, std::shared_ptr<Monoenergetic>, PrimaryEnergyDistribution>(m, "Monoenergetic")
    .def(init<double>())
    .def("pdf",&Monoenergetic::pdf)
    .def("SampleEnergy",&Monoenergetic::SampleEnergy)
    .def("GenerationProbability",&Monoenergetic::GenerationProbability)
    .def("Name",&Monoenergetic::Name);

  class_<PowerLaw, std::shared_ptr<PowerLaw>, PrimaryEnergyDistribution>(m, "PowerLaw")
    .def(init<double, double, double>())
    .def("pdf",&PowerLaw::pdf)
    .def("SampleEnergy",&PowerLaw::SampleEnergy)
    .def("GenerationProbability",&PowerLaw::GenerationProbability)
    .def("SetNormalizationAtEnergy",&PowerLaw::SetNormalizationAtEnergy)
    .def("Name",&PowerLaw::Name);

  class_<TabulatedFluxDistribution, std::shared_ptr<TabulatedFluxDistribution>, PrimaryEnergyDistribution>(m, "TabulatedFluxDistribution")
    .def(init<std::string, bool, bool>(), arg("fluxTableFilename"), arg("has_physical_normalization")=false, arg("romberg")=true)
    .def(init<double, double, std::string, bool, bool>(), arg("energyMin"), arg("energyMax"), arg("fluxTableFilename"), arg("has_physical_normalization")=false, arg("romberg")=true)
    .def(init<std::vector<double>, std::vector<double>, bool, bool>(), arg("energies"), arg("flux"), arg("has_physical_normalization")=false, arg("romberg")=true)
    .def(init<double, double, std::vector<double>, std::vector<double>, bool, bool>(), arg("energyMin"), arg("energyMax"), arg("energies"), arg("flux"), arg("has_physical_normalization")=false, arg("romberg")=true)
    .def("InverseCDF",&TabulatedFluxDistribution::InverseCDF)
    .def("SampleEnergy",&TabulatedFluxDistribution::SampleEnergy)
    .def("GenerationProbability",&TabulatedFluxDistribution::GenerationProbability)
    .def("SetEnergyBounds",&TabulatedFluxDistribution::SetEnergyBounds)
    .def("Name",&TabulatedFluxDistribution::Name)
    .def("GetIntegral",&TabulatedFluxDistribution::GetIntegral)
    .def("SamplePDF",&TabulatedFluxDistribution::SamplePDF)
    .def("SampleUnnormedPDF",&TabulatedFluxDistribution::SampleUnnormedPDF)
    .def("ComputeCDF",&TabulatedFluxDistribution::ComputeCDF)
    .def("GetCDF",&TabulatedFluxDistribution::GetCDF)
    .def("GetCDFEnergyNodes",&TabulatedFluxDistribution::GetCDFEnergyNodes)
    .def("GetEnergyNodes",&TabulatedFluxDistribution::GetEnergyNodes);

  class_<PiDARNuEDistribution, std::shared_ptr<PiDARNuEDistribution>, PrimaryEnergyDistribution>(m, "PiDARNuEDistribution")
    .def(init<>())
    .def("pdf",&PiDARNuEDistribution::pdf)
    .def("cdf",&PiDARNuEDistribution::cdf)
    .def("inv_cdf",&PiDARNuEDistribution::inv_cdf)
    .def("SampleEnergy",&PiDARNuEDistribution::SampleEnergy)
    .def("GenerationProbability",&PiDARNuEDistribution::GenerationProbability)
    .def("Name",&PiDARNuEDistribution::Name);

  // Energy Direction distributions

  class_<PrimaryEnergyDirectionDistribution, std::shared_ptr<PrimaryEnergyDirectionDistribution>, pyPrimaryEnergyDirectionDistribution, PrimaryInjectionDistribution, PhysicallyNormalizedDistribution>(m, "PrimaryEnergyDirectionDistribution",
    "Base class for joint energy-direction distributions. Subclass in python\n"
    "and implement SampleEnergyAndDirection(rand, detector_model,\n"
    "interactions, record) returning (energy, Vector3D), plus\n"
    "GenerationProbability.")
    .def(init<>())
    .def("Sample",&PrimaryEnergyDirectionDistribution::Sample)
    .def("SampleEnergyAndDirection",&PrimaryEnergyDirectionDistribution::SampleEnergyAndDirection)
    TrampolinePickleMethods(pyPrimaryEnergyDirectionDistribution);

  class_<Tabulated2DFluxDistribution, std::shared_ptr<Tabulated2DFluxDistribution>, PrimaryEnergyDirectionDistribution>(m, "Tabulated2DFluxDistribution")
    .def(init<std::string, bool>(), arg("fluxTableFilename"), arg("has_physical_normalization")=false)
    .def(init<double, double, std::string, bool>(), arg("energyMin"), arg("energyMax"), arg("fluxTableFilename"), arg("has_physical_normalization")=false)
    .def(init<double, double, double, double, std::string, bool>(), arg("energyMin"), arg("energyMax"), arg("cosZenithMin"), arg("cosZenithMax"), arg("fluxTableFilename"), arg("has_physical_normalization")=false)
    .def(init<std::vector<double>, std::vector<double>, std::vector<double>, bool>(), arg("energies"), arg("cosZeniths"), arg("flux"), arg("has_physical_normalization")=false)
    .def(init<double, double, std::vector<double>, std::vector<double>, std::vector<double>, bool>(), arg("energyMin"), arg("energyMax"), arg("energies"), arg("cosZeniths"), arg("flux"), arg("has_physical_normalization")=false)
    .def(init<double, double, double, double, std::vector<double>, std::vector<double>, std::vector<double>, bool>(), arg("energyMin"), arg("energyMax"), arg("cosZenithMin"), arg("cosZenithMax"), arg("energies"), arg("cosZeniths"), arg("flux"), arg("has_physical_normalization")=false)
    .def("SampleEnergyAndDirection",&Tabulated2DFluxDistribution::SampleEnergyAndDirection)
    .def("GenerationProbability",&Tabulated2DFluxDistribution::GenerationProbability)
    .def("SetEnergyBounds",&Tabulated2DFluxDistribution::SetEnergyBounds)
    .def("SetCosZenithBounds",&Tabulated2DFluxDistribution::SetCosZenithBounds)
    .def("Name",&Tabulated2DFluxDistribution::Name)
    .def("GetIntegral",&Tabulated2DFluxDistribution::GetIntegral)
    .def("SamplePDF",&Tabulated2DFluxDistribution::SamplePDF)
    .def("SampleUnnormedPDF",&Tabulated2DFluxDistribution::SampleUnnormedPDF)
    .def("GetEnergyNodes",&Tabulated2DFluxDistribution::GetEnergyNodes)
    .def("GetCosZenithNodes",&Tabulated2DFluxDistribution::GetCosZenithNodes);

  // Helicity distributions

  class_<PrimaryNeutrinoHelicityDistribution, std::shared_ptr<PrimaryNeutrinoHelicityDistribution>, PrimaryInjectionDistribution>(m, "PrimaryNeutrinoHelicityDistribution")
    .def(init<>())
    .def("Sample",&PrimaryNeutrinoHelicityDistribution::Sample)
    .def("GenerationProbability",&PrimaryNeutrinoHelicityDistribution::GenerationProbability)
    .def("DensityVariables",&PrimaryNeutrinoHelicityDistribution::DensityVariables)
    .def("Name",&PrimaryNeutrinoHelicityDistribution::Name);

  // Type distributions

  class_<PrimaryMass, std::shared_ptr<PrimaryMass>, PrimaryInjectionDistribution>(m, "PrimaryMass")
    .def(init<double>())
    .def("GetPrimaryMass",&PrimaryMass::GetPrimaryMass)
    .def("Sample",&PrimaryMass::Sample)
    .def("GenerationProbability",&PrimaryMass::GenerationProbability)
    .def("DensityVariables",&PrimaryMass::DensityVariables)
    .def("Name",&PrimaryMass::Name);

  // Vertex distributions

  class_<VertexPositionDistribution, std::shared_ptr<VertexPositionDistribution>, pyVertexPositionDistribution, PrimaryInjectionDistribution>(m, "VertexPositionDistribution",
    "Base class for vertex position distributions. Subclass in python and\n"
    "implement SamplePosition(rand, detector_model, interactions, record)\n"
    "returning (initial_position, interaction_vertex) as Vector3D objects,\n"
    "InjectionBounds(detector_model, interactions, record) returning the\n"
    "(start, end) bounds used for interaction probabilities, and\n"
    "GenerationProbability.")
    .def(init<>())
    .def("DensityVariables",&VertexPositionDistribution::DensityVariables)
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&VertexPositionDistribution::InjectionBounds, const_))
    .def("AreEquivalent",&VertexPositionDistribution::AreEquivalent)
    TrampolinePickleMethods(pyVertexPositionDistribution);

  // First, some range and depth functions

  class_<RangeFunction, std::shared_ptr<RangeFunction>>(m, "RangeFunction");
    //.def(init<>());
    //.def((siren::dataclasses::InteractionSignature const &, double));

  class_<DecayRangeFunction, std::shared_ptr<DecayRangeFunction>, RangeFunction>(m, "DecayRangeFunction")
    .def(init<double, double, double, double>())
    //.def((siren::dataclasses::InteractionSignature const &, double))
    .def("DecayLength",overload_cast<siren::dataclasses::ParticleType const &, double>(&DecayRangeFunction::DecayLength, const_))
    .def("DecayLength",overload_cast<double, double, double>(&DecayRangeFunction::DecayLength))
    .def("Range",&DecayRangeFunction::Range)
    .def("Multiplier",&DecayRangeFunction::Multiplier)
    .def("ParticleMass",&DecayRangeFunction::ParticleMass)
    .def("DecayWidth",&DecayRangeFunction::DecayWidth)
    .def("MaxDistance",&DecayRangeFunction::MaxDistance);

  class_<DepthFunction, std::shared_ptr<DepthFunction>>(m, "DepthFunction");
    //.def(init<>());
    //.def((siren::dataclasses::InteractionSignature const &, double));

  class_<LeptonDepthFunction, std::shared_ptr<LeptonDepthFunction>, DepthFunction>(m, "LeptonDepthFunction")
    .def(init<>())
    .def("__call__", &LeptonDepthFunction::operator())
    .def("SetMuParams",&LeptonDepthFunction::SetMuParams)
    .def("SetTauParams",&LeptonDepthFunction::SetTauParams)
    .def("SetScale",&LeptonDepthFunction::SetScale)
    .def("SetMaxDepth",&LeptonDepthFunction::SetMaxDepth)
    .def("GetMuAlpha",&LeptonDepthFunction::GetMuAlpha)
    .def("GetMuBeta",&LeptonDepthFunction::GetMuBeta)
    .def("GetTauAlpha",&LeptonDepthFunction::GetTauAlpha)
    .def("GetTauBeta",&LeptonDepthFunction::GetTauBeta)
    .def("GetScale",&LeptonDepthFunction::GetScale)
    .def("GetMaxDepth",&LeptonDepthFunction::GetMaxDepth)
    .def("GetLeptonDepthFunctionReturnValue",&LeptonDepthFunction::GetLeptonDepthFunctionReturnValue);
    //.def((siren::dataclasses::InteractionSignature const &, double));

  // VertexPositionDistribution subclasses

  class_<CylinderVolumePositionDistribution, std::shared_ptr<CylinderVolumePositionDistribution>, VertexPositionDistribution>(m, "CylinderVolumePositionDistribution")
    .def(init<siren::geometry::Cylinder>())
    .def("GenerationProbability",&CylinderVolumePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&CylinderVolumePositionDistribution::InjectionBounds)
    .def("Name",&CylinderVolumePositionDistribution::Name);

  class_<SphereVolumePositionDistribution, std::shared_ptr<SphereVolumePositionDistribution>, VertexPositionDistribution>(m, "SphereVolumePositionDistribution")
    .def(init<siren::geometry::Sphere>())
    .def("GenerationProbability",&SphereVolumePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&SphereVolumePositionDistribution::InjectionBounds)
    .def("Name",&SphereVolumePositionDistribution::Name);

  class_<ColumnDepthPositionDistribution, std::shared_ptr<ColumnDepthPositionDistribution>, VertexPositionDistribution>(m, "ColumnDepthPositionDistribution")
    .def(init<double, double, std::shared_ptr<DepthFunction>>())
    .def("GenerationProbability",&ColumnDepthPositionDistribution::GenerationProbability)
    .def("InjectionBounds",&ColumnDepthPositionDistribution::InjectionBounds)
    .def("Name",&ColumnDepthPositionDistribution::Name)
    .def("GetSamplePosition",&ColumnDepthPositionDistribution::GetSamplePosition);

  class_<DecayRangePositionDistribution, std::shared_ptr<DecayRangePositionDistribution>, VertexPositionDistribution>(m, "DecayRangePositionDistribution")
    .def(init<>())
    .def(init<double, double, std::shared_ptr<DecayRangeFunction>>())
    .def("GenerationProbability",&DecayRangePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&DecayRangePositionDistribution::InjectionBounds)
    .def("Name",&DecayRangePositionDistribution::Name);

  class_<PointSourcePositionDistribution, std::shared_ptr<PointSourcePositionDistribution>, VertexPositionDistribution>(m, "PointSourcePositionDistribution")
    .def(init<>())
    .def(init<siren::math::Vector3D, double>())
    .def(init([](std::array<double, 3> const & origin, double max_dist) {
        return PointSourcePositionDistribution(siren::math::Vector3D(origin), max_dist);
    }), arg("origin"), arg("max_distance"))
    .def("GenerationProbability",&PointSourcePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&PointSourcePositionDistribution::InjectionBounds)
    .def("Name",&PointSourcePositionDistribution::Name);

  class_<RangePositionDistribution, std::shared_ptr<RangePositionDistribution>, VertexPositionDistribution>(m, "RangePositionDistribution")
    .def(init<>())
    .def(init<double, double, std::shared_ptr<RangeFunction>, std::set<siren::dataclasses::ParticleType>>())
    .def("GenerationProbability",&RangePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&RangePositionDistribution::InjectionBounds)
    .def("Name",&RangePositionDistribution::Name);

  class_<PrimaryPhysicalVertexDistribution, std::shared_ptr<PrimaryPhysicalVertexDistribution>, VertexPositionDistribution>(m, "PrimaryPhysicalVertexDistribution")
    .def(init<>())
    .def("SamplePosition",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::PrimaryDistributionRecord &>(&PrimaryPhysicalVertexDistribution::SamplePosition, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&PrimaryPhysicalVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&PrimaryPhysicalVertexDistribution::InjectionBounds, const_))
    .def("Name",&PrimaryPhysicalVertexDistribution::Name);

  class_<PrimaryBoundedVertexDistribution, std::shared_ptr<PrimaryBoundedVertexDistribution>, VertexPositionDistribution>(m, "PrimaryBoundedVertexDistribution")
    .def(init<>())
    .def(init<double>())
    .def(init<std::shared_ptr<siren::geometry::Geometry>>())
    .def(init<std::shared_ptr<siren::geometry::Geometry>, double>())
    .def("SamplePosition",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::PrimaryDistributionRecord &>(&PrimaryBoundedVertexDistribution::SamplePosition, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&PrimaryBoundedVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&PrimaryBoundedVertexDistribution::InjectionBounds, const_))
    .def("Name",&PrimaryBoundedVertexDistribution::Name);

  class_<FixedTargetPositionDistribution, std::shared_ptr<FixedTargetPositionDistribution>, VertexPositionDistribution>(m, "FixedTargetPositionDistribution")
    .def(init<>())
    .def(init<siren::geometry::Cylinder, double>())
    .def(init<siren::geometry::Cylinder, std::shared_ptr<siren::geometry::Geometry>>())
    .def(init<siren::geometry::Cylinder, std::shared_ptr<siren::geometry::Geometry>, double>())
    .def("GenerationProbability",&FixedTargetPositionDistribution::GenerationProbability)
    .def("InjectionBounds",&FixedTargetPositionDistribution::InjectionBounds)
    .def("Name",&FixedTargetPositionDistribution::Name);

  // External distribution (inherits VertexPositionDistribution, must come after it)
  class_<PrimaryExternalDistribution, std::shared_ptr<PrimaryExternalDistribution>, VertexPositionDistribution, PhysicallyNormalizedDistribution>(m,"PrimaryExternalDistribution",
    "Primary distribution driven by an external CSV file.\n\n"
    "The first line is a comma-separated header naming each column and\n"
    "each subsequent line is one candidate event. Recognized columns:\n"
    "  x0, y0, z0  initial position of the primary\n"
    "  x, y, z     interaction vertex (also used as the initial position\n"
    "              when x0, y0, z0 are absent)\n"
    "  px, py, pz  primary momentum\n"
    "  E           primary energy; rows with E below emin are dropped\n"
    "  m           primary mass\n"
    "  t0          primary initial time, in SIREN time units where one\n"
    "              second is 1e9, so t0 is expressed in nanoseconds; it is\n"
    "              propagated to the vertex by time of flight\n"
    "  weight      physical weight of the row (primaries per unit exposure,\n"
    "              e.g. nimpwt/POT for dk2nu tables). Never affects how rows\n"
    "              are sampled: the supplied list is the intended injection\n"
    "              ensemble, drawn uniformly unless explicit sampling_weights\n"
    "              bias the selection. The weights encode the return to the\n"
    "              physical ensemble: PhysicalDensity reports the physical\n"
    "              row density on the physical side of the weight ratio and\n"
    "              the weight total is the distribution's physical\n"
    "              normalization, so absolute rates are correct with one\n"
    "              instance shared between the injection and physical sides.\n"
    "Any other column is stored as a named interaction parameter. Each of\n"
    "the x0/y0/z0, x/y/z, and px/py/pz groups must be given in full or\n"
    "omitted entirely.")
    .def(init<std::string>())
    .def(init<std::string, double>())
    .def(init<std::vector<std::string>, std::vector<std::vector<double>>>())
    .def(init<std::vector<std::string>, std::vector<std::vector<double>>, double>())
    .def(init<std::vector<std::string>, std::vector<std::vector<double>>, std::vector<double>>(),
         arg("keys"), arg("data"), arg("sampling_weights"))
    .def(init<std::vector<std::string>, std::vector<std::vector<double>>, std::vector<double>, double>(),
         arg("keys"), arg("data"), arg("sampling_weights"), arg("emin"))
    .def("Sample",&PrimaryExternalDistribution::Sample)
    .def("GetPhysicalNumEvents",&PrimaryExternalDistribution::GetPhysicalNumEvents)
    .def("DensityVariables",&PrimaryExternalDistribution::DensityVariables)
    .def("GenerationProbability",&PrimaryExternalDistribution::GenerationProbability)
    .def("Name",&PrimaryExternalDistribution::Name);

  class_<PrimaryAreaDistribution, std::shared_ptr<PrimaryAreaDistribution>, pyPrimaryAreaDistribution, PrimaryInjectionDistribution>(m, "PrimaryAreaDistribution",
    "Base class for primary area distributions. Subclass in python and\n"
    "implement SamplePointOfClosestApproach(rand, detector_model,\n"
    "interactions, record) returning a Vector3D, plus GenerationProbability.")
    .def(init<>())
    .def("DensityVariables",&PrimaryAreaDistribution::DensityVariables)
    .def("AreEquivalent",&PrimaryAreaDistribution::AreEquivalent)
    TrampolinePickleMethods(pyPrimaryAreaDistribution);

  class_<FixedTargetAreaDistribution, std::shared_ptr<FixedTargetAreaDistribution>, PrimaryAreaDistribution>(m, "FixedTargetAreaDistribution")
    .def(init<>())
    .def(init<siren::geometry::Cylinder>())
    .def("GenerationProbability",&FixedTargetAreaDistribution::GenerationProbability)
    .def("Name",&FixedTargetAreaDistribution::Name);

  class_<SecondaryInjectionDistribution, std::shared_ptr<SecondaryInjectionDistribution>, pySecondaryInjectionDistribution, WeightableDistribution>(m, "SecondaryInjectionDistribution",
    "Base class for secondary injection distributions. Subclass in python\n"
    "and implement Sample(rand, detector_model, interactions, record) and\n"
    "GenerationProbability(detector_model, interactions, record).")
    .def(init<>())
    .def("Sample",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::SecondaryDistributionRecord &>(&SecondaryInjectionDistribution::Sample, const_))
    .def("clone",&SecondaryInjectionDistribution::clone)
    TrampolinePickleMethods(pySecondaryInjectionDistribution);

  class_<SecondaryVertexPositionDistribution, std::shared_ptr<SecondaryVertexPositionDistribution>, pySecondaryVertexPositionDistribution, SecondaryInjectionDistribution>(m, "SecondaryVertexPositionDistribution",
    "Base class for secondary vertex position distributions. Subclass in\n"
    "python and implement SampleVertex(rand, detector_model, interactions,\n"
    "record), InjectionBounds(detector_model, interactions, record), and\n"
    "GenerationProbability.")
    .def(init<>())
    .def("DensityVariables",&SecondaryVertexPositionDistribution::DensityVariables)
    .def("AreEquivalent",&SecondaryVertexPositionDistribution::AreEquivalent)
    .def("Sample",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::SecondaryDistributionRecord &>(&SecondaryVertexPositionDistribution::Sample, const_))
    .def("SampleVertex",&SecondaryVertexPositionDistribution::SampleVertex)
    .def("InjectionBounds",&SecondaryVertexPositionDistribution::InjectionBounds)
    TrampolinePickleMethods(pySecondaryVertexPositionDistribution);

  class_<SecondaryPhysicalVertexDistribution, std::shared_ptr<SecondaryPhysicalVertexDistribution>, SecondaryVertexPositionDistribution>(m, "SecondaryPhysicalVertexDistribution")
    .def(init<>())
    .def("SampleVertex",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::SecondaryDistributionRecord &>(&SecondaryPhysicalVertexDistribution::SampleVertex, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryPhysicalVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryPhysicalVertexDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryPhysicalVertexDistribution::Name);

  class_<SecondaryBoundedVertexDistribution, std::shared_ptr<SecondaryBoundedVertexDistribution>, SecondaryVertexPositionDistribution>(m, "SecondaryBoundedVertexDistribution")
    .def(init<>())
    .def(init<double>())
    .def(init<std::shared_ptr<siren::geometry::Geometry>>())
    .def(init<std::shared_ptr<siren::geometry::Geometry>, double>())
    .def("SampleVertex",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::SecondaryDistributionRecord &>(&SecondaryBoundedVertexDistribution::SampleVertex, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryBoundedVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryBoundedVertexDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryBoundedVertexDistribution::Name);

  class_<SecondaryDecayRangePositionDistribution, std::shared_ptr<SecondaryDecayRangePositionDistribution>, SecondaryVertexPositionDistribution>(m, "SecondaryDecayRangePositionDistribution",
    "Bias a secondary interaction vertex by the probability that a collinear\n"
    "proxy daughter subsequently interacts or decays inside a fiducial volume.\n"
    "The current and daughter legs are evaluated in detector interaction depth,\n"
    "including arbitrary material density profiles and decay lengths.")
    .def(init<std::shared_ptr<siren::geometry::Geometry>,
              std::shared_ptr<siren::interactions::InteractionCollection>,
              double, double, double>(),
         arg("fiducial_volume"), arg("daughter_interactions"),
         arg("daughter_mass"), arg("daughter_energy_fraction") = 1.0,
         arg("max_length") = std::numeric_limits<double>::infinity(),
         keep_alive<1, 2>(), keep_alive<1, 3>())
    .def("SampleVertex",overload_cast<std::shared_ptr<siren::utilities::SIREN_random>, std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::SecondaryDistributionRecord &>(&SecondaryDecayRangePositionDistribution::SampleVertex, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryDecayRangePositionDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<siren::detector::DetectorModel const>, std::shared_ptr<siren::interactions::InteractionCollection const>, siren::dataclasses::InteractionRecord const &>(&SecondaryDecayRangePositionDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryDecayRangePositionDistribution::Name);
}
