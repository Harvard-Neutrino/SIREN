#include "SIREN/detector/CartesianAxis1D.h"

#include "SIREN/math/Vector3D.h"

namespace SI {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%      Cartesian     %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CartesianAxis1D::CartesianAxis1D() : Axis1D() {
    fAxis_.SetCartesianCoordinates(1, 0, 0);
    fp0_.SetCartesianCoordinates(0, 0, 0);
}

CartesianAxis1D::CartesianAxis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0) : Axis1D(fAxis, fp0) {}

bool CartesianAxis1D::compare(const Axis1D& ax) const {
    const CartesianAxis1D* c_ax = dynamic_cast<const CartesianAxis1D*>(&ax);
    if(!c_ax)
        return false;
    if(fp0_ != c_ax->fp0_)
        return false;
    if(fAxis_ != c_ax->fAxis_)
        return false;
    return true;
}

double CartesianAxis1D::GetX(const math::Vector3D& xi) const {
    return fAxis_ * (xi - fp0_);
}

double CartesianAxis1D::GetdX(const math::Vector3D& xi, const math::Vector3D& direction) const {
    (void)xi;

    return fAxis_ * direction;
}

} // namespace SI
} // namespace detector
