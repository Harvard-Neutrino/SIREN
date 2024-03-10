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
#include "SIREN/detector/PolynomialDistribution1D.h"

namespace siren {
namespace detector {

typedef DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> CartesianAxisPolynomialDensityDistribution;

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::CartesianAxisPolynomialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(siren::detector::CartesianAxisPolynomialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::DensityDistribution, siren::detector::CartesianAxisPolynomialDensityDistribution);

#endif // SIREN_CartesianAxisPolynomialDensityDistribution.h
