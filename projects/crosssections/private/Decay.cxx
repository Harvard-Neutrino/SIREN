#include "LeptonInjector/crosssections/Decay.h"

namespace LI {
namespace crosssections {

Decay::Decay() {}

bool Decay::operator==(Decay const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace crosssections
} // namespace LI

