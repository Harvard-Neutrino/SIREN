#include "SIREN/detector/DensityDistribution.h"

namespace siren {
namespace detector {

DensityDistribution::DensityDistribution() {}

DensityDistribution::DensityDistribution(const DensityDistribution& density_distr) {}

bool DensityDistribution::operator==(const DensityDistribution& dens_distr) const
{
    if (!this->compare(dens_distr) )
        return false;
    return true;
}


bool DensityDistribution::operator!=(const DensityDistribution& dens_distr) const {
    return !(*this == dens_distr);
}

} // namespace siren
} // namespace detector

CEREAL_REGISTER_DYNAMIC_INIT(siren_DensityDistribution);

