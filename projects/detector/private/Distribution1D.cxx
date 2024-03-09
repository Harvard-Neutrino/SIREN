#include "SIREN/detector/Distribution1D.h"

namespace siren {
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

} // namespace siren
} // namespace detector
