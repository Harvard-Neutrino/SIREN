#include "SIREN/detector/RadialAxisPolynomialDensityDistribution.h"

// Explicit template instantiation
template class siren::detector::DensityDistribution1D<
    siren::detector::RadialAxis1D,
    siren::detector::PolynomialDistribution1D>;

CEREAL_REGISTER_DYNAMIC_INIT(siren_RadialAxisPolynomialDensityDistribution);

