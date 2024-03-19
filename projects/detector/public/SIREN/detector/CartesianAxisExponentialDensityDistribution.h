#pragma once
#ifndef SIREN_CartesianAxisDensityDistribution_H
#define SIREN_CartesianAxisDensityDistribution_H
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/CartesianAxis1D.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/DensityDistribution1D.h"
#include "SIREN/detector/ExponentialDistribution1D.h"

namespace siren {
namespace detector {

typedef DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> CartesianAxisExponentialDensityDistribution;

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::CartesianAxisExponentialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(siren::detector::CartesianAxisExponentialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::DensityDistribution, siren::detector::CartesianAxisExponentialDensityDistribution);

#endif // SIREN_CartesianAxisExponentialDensityDistribution.h
