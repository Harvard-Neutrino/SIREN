#pragma once
#ifndef LI_RadialAxisExponentialDensityDistribution_H
#define LI_RadialAxisExponentialDensityDistribution_H
#include <memory>
#include <string>
#include <exception>
#include <functional>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/detector/Axis1D.h"
#include "LeptonInjector/detector/RadialAxis1D.h"
#include "LeptonInjector/detector/DensityDistribution.h"
#include "LeptonInjector/detector/DensityDistribution1D.h"
#include "LeptonInjector/detector/ExponentialDistribution1D.h"

namespace LI {
namespace detector {

typedef DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D> RadialAxisExponentialDensityDistribution;

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::RadialAxisExponentialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(LI::detector::RadialAxisExponentialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::detector::DensityDistribution, LI::detector::RadialAxisExponentialDensityDistribution);

#endif // LI_RadialAxisExponentialDensityDistribution.h

