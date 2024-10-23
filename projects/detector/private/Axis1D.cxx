#include "SIREN/detector/Axis1D.h"

#include "SIREN/math/Vector3D.h"

namespace siren {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%        Axis        %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Axis1D::Axis1D() {}

Axis1D::Axis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0) : fAxis_(fAxis), fp0_(fp0) {}

Axis1D::Axis1D(const Axis1D& axis) : fAxis_(axis.fAxis_), fp0_(axis.fp0_) {}

bool Axis1D::operator==(const Axis1D& axis) const {
    if (!this->compare(axis) )
        return false;
    return true;
}

bool Axis1D::operator!=(const Axis1D& axis) const {
    return !(*this == axis);
}

} // namespace siren
} // namespace detector

CEREAL_REGISTER_DYNAMIC_INIT(siren_Axis1D);

