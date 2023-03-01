#pragma once
#ifndef LI_CartesianAxisDensityDistribution_H
#define LI_CartesianAxisDensityDistribution_H
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "LeptonInjector/detector/Axis1D.h"
#include "LeptonInjector/detector/CartesianAxis1D.h"
#include "LeptonInjector/detector/DensityDistribution.h"
#include "LeptonInjector/detector/DensityDistribution1D.h"
#include "LeptonInjector/detector/PolynomialDistribution1D.h"

namespace LI {
namespace detector {

typedef DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> CartesianAxisPolynomialDensityDistribution;

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::CartesianAxisPolynomialDensityDistribution, 0);
CEREAL_REGISTER_TYPE(LI::detector::CartesianAxisPolynomialDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::detector::DensityDistribution, LI::detector::CartesianAxisPolynomialDensityDistribution);

#endif // LI_CartesianAxisPolynomialDensityDistribution.h
