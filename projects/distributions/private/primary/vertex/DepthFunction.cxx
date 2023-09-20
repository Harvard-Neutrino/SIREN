#include "LeptonInjector/distributions/primary/vertex/DepthFunction.h"

#include <typeinfo>   // for type_info
#include <typeindex>  // for type_index

namespace LI { namespace dataclasses { struct InteractionSignature; } }

namespace LI {
namespace distributions {

//---------------
// class DepthFunction
//---------------
DepthFunction::DepthFunction() {}

double DepthFunction::operator()(LI::dataclasses::InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool DepthFunction::operator==(DepthFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool DepthFunction::operator<(DepthFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

} // namespace distributions
} // namespace LI
