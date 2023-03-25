
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
#include "../../public/LeptonInjector/distributions/primary/type/PrimaryInjector.h"
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
#include "../../public/LeptonInjector/distributions/primary/vertex/SecondaryPositionDistribution.h"

#include "../../../utilities/public/LeptonInjector/utilities/Random.h"
#include "../../../detector/public/LeptonInjector/detector/EarthModel.h"
#include "../../../crosssections/public/LeptonInjector/crosssections/CrossSectionCollection.h"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

using namespace pybind11;

/*PYBIND11_MODULE(Utilities,m) {
  using namespace LI::utilities;

  class_<LI_random>(m, "LI_random")
    .def(init<>())
    .def(init<unsigned int>())
    .def("Uniform",&LI_random::Uniform)
    .def("set_seed",&LI_random::set_seed);
}*/

PYBIND11_MODULE(Distributions,m) {
  using namespace LI::distributions;

  class_<PhysicallyNormalizedDistribution>(m, "PhysicallyNormalizedDistribution")
    .def(init<>())
    .def(init<double>())
    .def_property("normalization",&PhysicallyNormalizedDistribution::GetNormalization,&PhysicallyNormalizedDistribution::SetNormalization)
    .def("IsNormalizationSet",&PhysicallyNormalizedDistribution::IsNormalizationSet);

  class_<WeightableDistribution>(m, "WeightableDistribution")
    .def("DensityVariables",&WeightableDistribution::DensityVariables)
    .def("AreEquivalent",&WeightableDistribution::AreEquivalent);

  class_<NormalizationConstant, WeightableDistribution, PhysicallyNormalizedDistribution>(m, "NormalizationConstant")
    .def(init<double>())
    .def("GenerationProbability",&NormalizationConstant::GenerationProbability)
    .def("Name",&NormalizationConstant::Name);

  class_<InjectionDistribution, WeightableDistribution>(m, "InjectionDistribution")
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord &>(&InjectionDistribution::Sample, const_))
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionTreeDatum &>(&InjectionDistribution::Sample, const_))
    .def("IsPositionDistribution",&InjectionDistribution::IsPositionDistribution);

  // Direciton distributions
  
  class_<PrimaryDirectionDistribution, InjectionDistribution>(m, "PrimaryDirectionDistribution")
    .def("Sample",&PrimaryDirectionDistribution::Sample)
    .def("DensityVariables",&PrimaryDirectionDistribution::DensityVariables)
    .def("GenerationProbability",&PrimaryDirectionDistribution::GenerationProbability);

  class_<Cone, PrimaryDirectionDistribution>(m, "Cone")
    .def(init<LI::math::Vector3D, double>());
    //.def("GenerationProbability",&Cone::GenerationProbability);
  
  class_<IsotropicDirection, PrimaryDirectionDistribution>(m, "IsotropicDirection")
    .def(init<>());
    //.def("GenerationProbability",&IsotropicDirection::GenerationProbability);
  
  class_<FixedDirection, PrimaryDirectionDistribution>(m, "FixedDirection")
    .def(init<LI::math::Vector3D>());
    //.def("GenerationProbability",&FixedDirection::GenerationProbability);

  // Energy distributions

  class_<PrimaryEnergyDistribution, InjectionDistribution>(m, "PrimaryEnergyDistribution")
    .def("Sample",&PrimaryEnergyDistribution::Sample);

  class_<Monoenergetic, PrimaryEnergyDistribution>(m, "Monoenergetic")
    .def(init<double>())
    .def("pdf",&Monoenergetic::pdf)
    .def("SampleEnergy",&Monoenergetic::SampleEnergy)
    .def("GenerationProbability",&Monoenergetic::GenerationProbability)
    .def("Name",&Monoenergetic::Name);
  
  class_<PowerLaw, PrimaryEnergyDistribution>(m, "PowerLaw")
    .def(init<double, double, double>())
    .def("pdf",&PowerLaw::pdf)
    .def("SampleEnergy",&PowerLaw::SampleEnergy)
    .def("GenerationProbability",&PowerLaw::GenerationProbability)
    .def("SetNormalizationAtEnergy",&PowerLaw::SetNormalizationAtEnergy)
    .def("Name",&PowerLaw::Name);

  class_<TabulatedFluxDistribution, PrimaryEnergyDistribution>(m, "TabulatedFluxDistribution")
    .def(init<std::string, bool>())
    .def(init<double, double, std::string, bool>())
    .def("SampleEnergy",&TabulatedFluxDistribution::SampleEnergy)
    .def("GenerationProbability",&TabulatedFluxDistribution::GenerationProbability)
    .def("SetEnergyBounds",&TabulatedFluxDistribution::SetEnergyBounds)
    .def("Name",&TabulatedFluxDistribution::Name);

  // Helicity distributions
  
  class_<PrimaryNeutrinoHelicityDistribution, InjectionDistribution>(m, "PrimaryNeutrinoHelicityDistribution")
    .def(init<>())
    .def("Sample",&PrimaryNeutrinoHelicityDistribution::Sample)
    .def("GenerationProbability",&PrimaryNeutrinoHelicityDistribution::GenerationProbability)
    .def("DensityVariables",&PrimaryNeutrinoHelicityDistribution::DensityVariables)
    .def("Name",&PrimaryNeutrinoHelicityDistribution::Name);

  // Type distributions
  
  class_<PrimaryInjector, InjectionDistribution>(m, "PrimaryInjector")
    .def(init<LI::dataclasses::Particle::ParticleType, double>())
    .def("PrimaryMass",&PrimaryInjector::PrimaryMass)
    .def("Sample",&PrimaryInjector::Sample)
    .def("GenerationProbability",&PrimaryInjector::GenerationProbability)
    .def("DensityVariables",&PrimaryInjector::DensityVariables)
    .def("Name",&PrimaryInjector::Name);

  // Vertex distributions

  class_<VertexPositionDistribution, InjectionDistribution>(m, "VertexPositionDistribution")
    .def("IsPositionDistribution",&VertexPositionDistribution::IsPositionDistribution)
    .def("DensityVariables",&VertexPositionDistribution::DensityVariables)
    //.def("InjectionBounds",&VertexPositionDistribution::InjectionBounds)
    .def("AreEquivalent",&VertexPositionDistribution::AreEquivalent);

  // First, some range and depth functions

  class_<RangeFunction>(m, "RangeFunction");
    //.def(init<>());
    //.def((LI::dataclasses::InteractionSignature const &, double));
  
  class_<DecayRangeFunction, RangeFunction>(m, "DecayRangeFunction")
    .def(init<double, double, double, double>())
    //.def((LI::dataclasses::InteractionSignature const &, double))
    .def("DecayLength",overload_cast<LI::dataclasses::InteractionSignature const &, double>(&DecayRangeFunction::DecayLength, const_))
    .def("DecayLength",overload_cast<double, double, double>(&DecayRangeFunction::DecayLength))
    .def("Range",&DecayRangeFunction::Range)
    .def("Multiplier",&DecayRangeFunction::Multiplier)
    .def("ParticleMass",&DecayRangeFunction::ParticleMass)
    .def("DecayWidth",&DecayRangeFunction::DecayWidth)
    .def("MaxDistance",&DecayRangeFunction::MaxDistance);
  
  class_<DepthFunction>(m, "DepthFunction");
    //.def(init<>());
    //.def((LI::dataclasses::InteractionSignature const &, double));

  class_<LeptonDepthFunction, DepthFunction>(m, "LeptonDepthFunction")
    .def(init<>())
    .def("SetMuParams",&LeptonDepthFunction::SetMuParams)
    .def("SetTauParams",&LeptonDepthFunction::SetTauParams)
    .def("SetScale",&LeptonDepthFunction::SetScale)
    .def("SetMaxDepth",&LeptonDepthFunction::SetMaxDepth)
    .def("GetMuAlpha",&LeptonDepthFunction::GetMuAlpha)
    .def("GetMuBeta",&LeptonDepthFunction::GetMuBeta)
    .def("GetTauAlpha",&LeptonDepthFunction::GetTauAlpha)
    .def("GetTauBeta",&LeptonDepthFunction::GetTauBeta)
    .def("GetScale",&LeptonDepthFunction::GetScale)
    .def("GetMaxDepth",&LeptonDepthFunction::GetMaxDepth);
    //.def((LI::dataclasses::InteractionSignature const &, double));

  // VertexPositionDistribution subclasses

  class_<CylinderVolumePositionDistribution, VertexPositionDistribution>(m, "CylinderVolumePositionDistribution")
    .def(init<LI::geometry::Cylinder>())
    .def("GenerationProbability",&CylinderVolumePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&CylinderVolumePositionDistribution::InjectionBounds)
    .def("Name",&CylinderVolumePositionDistribution::Name);
  
  class_<ColumnDepthPositionDistribution, VertexPositionDistribution>(m, "ColumnDepthPositionDistribution")
    .def(init<double, double, std::shared_ptr<DepthFunction>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def("GenerationProbability",&ColumnDepthPositionDistribution::GenerationProbability)
    .def("InjectionBounds",&ColumnDepthPositionDistribution::InjectionBounds)
    .def("Name",&ColumnDepthPositionDistribution::Name);
  
  class_<DecayRangePositionDistribution, VertexPositionDistribution>(m, "DecayRangePositionDistribution")
    .def(init<>())
    .def(init<double, double, std::shared_ptr<DecayRangeFunction>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def("GenerationProbability",&DecayRangePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&DecayRangePositionDistribution::InjectionBounds)
    .def("Name",&DecayRangePositionDistribution::Name);
  
  class_<PointSourcePositionDistribution, VertexPositionDistribution>(m, "PointSourcePositionDistribution")
    .def(init<>())
    .def(init<LI::math::Vector3D, double, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def("GenerationProbability",&PointSourcePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&PointSourcePositionDistribution::InjectionBounds)
    .def("Name",&PointSourcePositionDistribution::Name);
  
  class_<RangePositionDistribution, VertexPositionDistribution>(m, "RangePositionDistribution")
    .def(init<>())
    .def(init<double, double, std::shared_ptr<RangeFunction>, std::set<LI::dataclasses::Particle::ParticleType>>())
    .def("GenerationProbability",&RangePositionDistribution::GenerationProbability)
    .def("InjectionBounds",&RangePositionDistribution::InjectionBounds)
    .def("Name",&RangePositionDistribution::Name);
  
  class_<SecondaryPositionDistribution, VertexPositionDistribution>(m, "SecondaryPositionDistribution")
    .def(init<>())
    .def(init<double>())
    .def(init<double, std::shared_ptr<LI::geometry::Geometry>>())
    .def(init<std::shared_ptr<LI::geometry::Geometry>>())
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord &>(&SecondaryPositionDistribution::Sample, const_))
    .def("Sample",overload_cast<std::shared_ptr<LI::utilities::LI_random>, std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionTreeDatum &>(&SecondaryPositionDistribution::Sample, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryPositionDistribution::GenerationProbability, const_))
    .def("GenerationProbability",overload_cast<std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionTreeDatum const &>(&SecondaryPositionDistribution::GenerationProbability, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionRecord const &>(&SecondaryPositionDistribution::InjectionBounds, const_))
    .def("InjectionBounds",overload_cast<std::shared_ptr<LI::detector::EarthModel const>, std::shared_ptr<LI::crosssections::CrossSectionCollection const>, LI::dataclasses::InteractionTreeDatum const &>(&SecondaryPositionDistribution::InjectionBounds, const_))
    .def("Name",&SecondaryPositionDistribution::Name);

}


