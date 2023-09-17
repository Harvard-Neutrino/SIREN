#include "LeptonInjector/detector/Distribution1D.h"

#include "LeptonInjector/math/Vector3D.h"

namespace LI {
namespace detector {

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%    Distribution    %%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool Distribution1D::operator==(const Distribution1D& dist) const {
    if (!this->compare(dist) )
        return false;
    return true;
}

bool Distribution1D::operator!=(const Distribution1D& dist) const {
    return !(*this == dist);
}

} // namespace LI
} // namespace detector
