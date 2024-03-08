#pragma once
#ifndef LI_CartesianAxisDensityDistribution_H
#define LI_CartesianAxisDensityDistribution_H
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/CartesianAxis1D.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/DensityDistribution1D.h"
#include "SIREN/detector/PolynomialDistribution1D.h"

namespace SI {
namespace detector {

typedef DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> CartesianAxisPolynomialDensityDistribution;

} // namespace detector
} // namespace SI

CEREAL_CLASS_VERSION(SI::detector::CartesianAxisPolynomialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(SI::detector::CartesianAxisPolynomialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(SI::detector::DensityDistribution, SI::detector::CartesianAxisPolynomialDensityDistribution);

#endif // LI_CartesianAxisPolynomialDensityDistribution.h
