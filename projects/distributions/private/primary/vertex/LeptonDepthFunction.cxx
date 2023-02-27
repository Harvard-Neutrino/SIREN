#include <tuple>
#include "LeptonInjector/distributions/primary/vertex/DepthFunction.h"
#include "LeptonInjector/distributions/primary/vertex/LeptonDepthFunction.h"

namespace LI {
namespace distributions {

//---------------
// class LeptonDepthFunction : DepthFunction
//---------------
LeptonDepthFunction::LeptonDepthFunction() {}

void LeptonDepthFunction::SetMuParams(double mu_alpha, double mu_beta) {
    this->mu_alpha = mu_alpha;
    this->mu_beta = mu_beta;
}

void LeptonDepthFunction::SetTauParams(double tau_alpha, double tau_beta) {
    this->tau_alpha = tau_alpha;
    this->tau_beta = tau_beta;

}

void LeptonDepthFunction::SetScale(double scale) {
    this->scale = scale;
}

void LeptonDepthFunction::SetMaxDepth(double max_depth) {
    this->max_depth = max_depth;
}

double LeptonDepthFunction::GetMuAlpha() const {
    return mu_alpha;
}

double LeptonDepthFunction::GetMuBeta() const {
    return mu_beta;
}

double LeptonDepthFunction::GetTauAlpha() const {
    return tau_alpha;
}

double LeptonDepthFunction::GetTauBeta() const {
    return tau_beta;
}

double LeptonDepthFunction::GetScale() const {
    return scale;
}

double LeptonDepthFunction::GetMaxDepth() const {
    return max_depth;
}

double LeptonDepthFunction::operator()(LI::crosssections::InteractionSignature const & signature, double energy) const {
    double range = log(1.0 + energy * mu_beta / mu_alpha) / mu_beta;
    if(tau_primaries.count(signature.primary_type) > 0)
        range += log(1.0 + energy * tau_beta / tau_alpha) / tau_beta;
    return std::min(range, max_depth);
}

bool LeptonDepthFunction::equal(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        ==
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

bool LeptonDepthFunction::less(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        <
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

} // namespace distributions
} // namespace LI
