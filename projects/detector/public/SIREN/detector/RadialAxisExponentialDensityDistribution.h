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

#include "SIREN/math/Vector3D.h"

#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/RadialAxis1D.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/DensityDistribution1D.h"
#include "SIREN/detector/ExponentialDistribution1D.h"

namespace SI {
namespace detector {

typedef DensityDistribution1D<RadialAxis1D,ExponentialDistribution1D> RadialAxisExponentialDensityDistribution;

} // namespace detector
} // namespace SI

CEREAL_CLASS_VERSION(SI::detector::RadialAxisExponentialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(SI::detector::RadialAxisExponentialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::detector::DensityDistribution, SI::detector::RadialAxisExponentialDensityDistribution);

#endif // LI_RadialAxisExponentialDensityDistribution.h

