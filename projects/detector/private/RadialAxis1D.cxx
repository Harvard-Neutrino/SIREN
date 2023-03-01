#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/detector/RadialAxis1D.h"

namespace LI {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%       Radial       %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RadialAxis1D::RadialAxis1D() : Axis1D() {
    fp0_.SetCartesianCoordinates(0, 0, 0);
    fAxis_.SetSphericalCoordinates(1, 0, 0);
}

bool RadialAxis1D::compare(const Axis1D& ax) const {
    const RadialAxis1D* r_ax = dynamic_cast<const RadialAxis1D*>(&ax);
    if(!r_ax)
        return false;
    if(fp0_ != r_ax->fp0_)
        return false;
    return true;
}

RadialAxis1D::RadialAxis1D(const math::Vector3D& fAxis, const math::Vector3D& fp0) : Axis1D(fAxis, fp0) {}

RadialAxis1D::RadialAxis1D(const math::Vector3D& fp0) : Axis1D(math::Vector3D(), fp0) {}

double RadialAxis1D::GetX(const math::Vector3D& xi) const {
    return (xi - fp0_).magnitude();
}

double RadialAxis1D::GetdX(const math::Vector3D& xi, const math::Vector3D& direction) const {
    math::Vector3D aux{xi - fp0_};
    aux.normalize();

    return aux * direction;
}

} // namespace LI
} // namespace detector
