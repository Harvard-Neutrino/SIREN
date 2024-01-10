#include "LeptonInjector/interactions/CrossSection.h"

namespace LI {
namespace interactions {

CrossSection::CrossSection() {}

bool CrossSection::operator==(CrossSection const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace interactions
} // namespace LI

