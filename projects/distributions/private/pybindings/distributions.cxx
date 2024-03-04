
#include <vector>
#include <string>

#include "../../public/LeptonInjector/distributions/Distributions.h"
#include "../../public/LeptonInjector/distributions/primary/direction/PrimaryDirectionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/direction/Cone.h"
#include "../../public/LeptonInjector/distributions/primary/direction/FixedDirection.h"
#include "../../public/LeptonInjector/distributions/primary/direction/IsotropicDirection.h"
#include "../../public/LeptonInjector/distributions/primary/energy/Monoenergetic.h"
#include "../../public/LeptonInjector/distributions/primary/energy/PowerLaw.h"
#include "../../public/LeptonInjector/distributions/primary/energy/TabulatedFluxDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/helicity/PrimaryNeutrinoHelicityDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/mass/PrimaryMass.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/VertexPositionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/ColumnDepthPositionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/CylinderVolumePositionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/DecayRangeFunction.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/DecayRangePositionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/DepthFunction.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/LeptonDepthFunction.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/PointSourcePositionDistribution.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/RangeFunction.h"
#include "../../public/LeptonInjector/distributions/primary/vertex/RangePositionDistribution.h"
#include "../../public/LeptonInjector/distributions/secondary/vertex/SecondaryPhysicalVertexDistribution.h"
#include "../../public/LeptonInjector/distributions/secondary/vertex/SecondaryBoundedVertexDistribution.h"

#include "../../../utilities/public/LeptonInjector/utilities/Random.h"
#include "../../../detector/public/LeptonInjector/detector/DetectorModel.h"
#include "../../../interactions/public/LeptonInjector/interactions/InteractionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

PYBIND11_DECLARE_HOLDER_TYPE(T__,std::shared_ptr<T__>)

using namespace pybind11;

PYBIND11_MODULE(distributions,m) {
  using namespace LI::distributions;

  class_<PhysicallyNormalizedDistribution, std::shared_ptr<PhysicallyNormalizedDistribution>>(m, "PhysicallyNormalizedDistribution")
    .def(init<>())
    .def(init<double>())
    .def_property("normalization",&PhysicallyNormalizedDistribution::GetNormalization,&PhysicallyNormalizedDistribution::SetNormalization)
    .def("IsNormalizationSet",&PhysicallyNormalizedDistribution::IsNormalizationSet);

  class_<WeightableDistribution, std::shared_ptr<WeightableDistribution>>(m, "WeightableDistribution")
    .def("DensityVariables",&WeightableDistribution::DensityVariables)
    .def("AreEquivalent",&WeightableDistribution::AreEquivalent);

  class_<NormalizationConstant, std::shared_ptr<NormalizationConstant>, WeightableDistribution, PhysicallyNormalizedDistribution>(m, "NormalizationConstant")
    .def(init<double>())
    .def("GenerationProbability",&NormalizationConstant::GenerationProbability)
    .def("Name",&NormalizationConstant::Name);

  class_<PrimaryInjectionDistribution, std::shared_ptr<PrimaryInjectionDistribution>, WeightableDistribution>(m, "PrimaryInjectionDistribution")
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::PrimaryDistributionRecord &>(&PrimaryInjectionDistribution::Sample, const_))
    ;

  // Direciton distributions

  class_<PrimaryDirectionDistribution, std::shared_ptr<PrimaryDirectionDistribution>, PrimaryInjectionDistribution>(m, "PrimaryDirectionDistribution")
    .def("Sample",&PrimaryDirectionDistribution::Sample)
    .def("DensityVariables",&PrimaryDirectionDistribution::DensityVariables)
    .def("GenerationProbability",&PrimaryDirectionDistribution::GenerationProbability);

  class_<Cone, std::shared_ptr<Cone>, PrimaryDirectionDistribution>(m, "Cone")
    .def(init<LI::math::Vector3D, double>());
    //.def("GenerationProbability",&Cone::GenerationProbability);

  class_<IsotropicDirection, std::shared_ptr<IsotropicDirection>, PrimaryDirectionDistribution>(m, "IsotropicDirection")
    .def(init<>());
    //.def("GenerationProbability",&IsotropicDirection::GenerationProbability);

  class_<FixedDirection, std::shared_ptr<FixedDirection>, PrimaryDirectionDistribution>(m, "FixedDirection")
    .def(init<LI::math::Vector3D>());
    //.def("GenerationProbability",&FixedDirection::GenerationProbability);

  // Energy distributions

  class_<PrimaryEnergyDistribution, std::shared_ptr<PrimaryEnergyDistribution>, PrimaryInjectionDistribution, PhysicallyNormalizedDistribution>(m, "PrimaryEnergyDistribution")
    .def("Sample",&PrimaryEnergyDistribution::Sample);

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
    .def(init<std::string, bool>())
    .def(init<double, double, std::string, bool>())
    .def(init<std::vector<double>, std::vector<double>, bool>())
    .def(init<double, double, std::vector<double>, std::vector<double>, bool>())
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

  class_<VertexPositionDistribution, std::shared_ptr<VertexPositionDistribution>, PrimaryInjectionDistribution>(m, "VertexPositionDistribution")
    .def("DensityVariables",&VertexPositionDistribution::DensityVariables)
    //.def("InjectionBounds",&VertexPositionDistribution::InjectionBounds)
    .def("AreEquivalent",&VertexPositionDistribution::AreEquivalent);

  // First, some range and depth functions

  class_<RangeFunction, std::shared_ptr<RangeFunction>>(m, "RangeFunction");
    //.def(init<>());
    //.def((LI::dataclasses::InteractionSignature const &, double));

  class_<DecayRangeFunction, std::shared_ptr<DecayRangeFunction>, RangeFunction>(m, "DecayRangeFunction")
    .def(init<double, double, double, double>())
    //.def((LI::dataclasses::InteractionSignature const &, double))
    .def("DecayLength",overload_cast<LI::dataclasses::ParticleType const &, double>(&DecayRangeFunction::DecayLength, const_))
    .def("DecayLength",overload_cast<double, double, double>(&DecayRangeFunction::DecayLength))
    .def("Range",&DecayRangeFunction::Range)
    .def("Multiplier",&DecayRangeFunction::Multiplier)
    .def("ParticleMass",&DecayRangeFunction::ParticleMass)
    .def("DecayWidth",&DecayRangeFunction::DecayWidth)
    .def("MaxDistance",&DecayRangeFunction::MaxDistance);

  class_<DepthFunction, std::shared_ptr<DepthFunction>>(m, "DepthFunction");
    //.def(init<>());
    //.def((LI::dataclasses::InteractionSignature const &, double));

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
    //.def((LI::dataclasses::InteractionSignature const &, double));

  // VertexPositionDistribution subclasses

  class_<CylinderVolumePositionDistribution, std::shared_ptr<CylinderVolumePositionDistribution>, VertexPositionDistribution>(m, "CylinderVolumePositionDistribution")
    .def(init<LI::geometry::Cylinder>())
    .def("GenerationProbability",&CylinderVolumePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&CylinderVolumePositionDistribution::InjectionBounds)
    .def("Name",&CylinderVolumePositionDistribution::Name);

  class_<ColumnDepthPositionDistribution, std::shared_ptr<ColumnDepthPositionDistribution>, VertexPositionDistribution>(m, "ColumnDepthPositionDistribution")
    .def(init<double, double, std::shared_ptr<DepthFunction>, std::set<LI::dataclasses::ParticleType>>())
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
    .def(init<LI::math::Vector3D, double, std::set<LI::dataclasses::ParticleType>>())
    .def("GenerationProbability",&PointSourcePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&PointSourcePositionDistribution::InjectionBounds)
    .def("Name",&PointSourcePositionDistribution::Name);

  class_<RangePositionDistribution, std::shared_ptr<RangePositionDistribution>, VertexPositionDistribution>(m, "RangePositionDistribution")
    .def(init<>())
    .def(init<double, double, std::shared_ptr<RangeFunction>, std::set<LI::dataclasses::ParticleType>>())
    .def("GenerationProbability",&RangePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&RangePositionDistribution::InjectionBounds)
    .def("Name",&RangePositionDistribution::Name);


  class_<SecondaryInjectionDistribution, std::shared_ptr<SecondaryInjectionDistribution>, WeightableDistribution>(m, "SecondaryInjectionDistribution")
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::SecondaryDistributionRecord &>(&SecondaryInjectionDistribution::Sample, const_));

  class_<SecondaryVertexPositionDistribution, std::shared_ptr<SecondaryVertexPositionDistribution>, SecondaryInjectionDistribution>(m, "SecondaryVertexPositionDistribution")
    .def("DensityVariables",&SecondaryVertexPositionDistribution::DensityVariables)
    .def("AreEquivalent",&SecondaryVertexPositionDistribution::AreEquivalent)
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::SecondaryDistributionRecord &>(&SecondaryVertexPositionDistribution::Sample, const_));

  class_<SecondaryPhysicalVertexDistribution, std::shared_ptr<SecondaryPhysicalVertexDistribution>, SecondaryVertexPositionDistribution>(m, "SecondaryPhysicalVertexDistribution")
    .def(init<>())
    .def("SampleVertex",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::SecondaryDistributionRecord &>(&SecondaryPhysicalVertexDistribution::SampleVertex, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryPhysicalVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryPhysicalVertexDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryPhysicalVertexDistribution::Name);

  class_<SecondaryBoundedVertexDistribution, std::shared_ptr<SecondaryBoundedVertexDistribution>, SecondaryVertexPositionDistribution>(m, "SecondaryBoundedVertexDistribution")
    .def(init<>())
    .def(init<double>())
    .def(init<std::shared_ptr<LI::geometry::Geometry>>())
    .def(init<std::shared_ptr<LI::geometry::Geometry>, double>())
    .def("SampleVertex",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::SecondaryDistributionRecord &>(&SecondaryBoundedVertexDistribution::SampleVertex, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryBoundedVertexDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<LI::detector::DetectorModel const>, std::shared_ptr<LI::interactions::InteractionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryBoundedVertexDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryBoundedVertexDistribution::Name);
}


