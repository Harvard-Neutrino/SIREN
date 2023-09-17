#include "LeptonInjector/crosssections/CrossSection.h"

#include <vector>

namespace LI {
namespace crosssections {

CrossSection::CrossSection() {}

bool CrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace crosssections
} // namespace LI

