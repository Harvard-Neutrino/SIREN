#include "SIREN/detector/ConstantDensityDistribution.h"

// Explicit template instantiation
template class siren::detector::DensityDistribution1D<
    siren::detector::CartesianAxis1D,
    siren::detector::ConstantDistribution1D>;

CEREAL_REGISTER_DYNAMIC_INIT(siren_ConstantDensityDistribution);

