#include "SIREN/interactions/Hadronization.h"
#include "SIREN/dataclasses/InteractionRecord.h"

namespace siren {
namespace interactions {

Hadronization::Hadronization() {}

bool Hadronization::operator==(Hadronization const & other) const {
    if(this == &other)
        return true;
    else
        return this->equal(other);
}

} // namespace interactions
} // namespace siren