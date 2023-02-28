#include "LeptonInjector/dataclasses/InteractionSignature.h"
#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/vertex/DepthFunction.h"

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
