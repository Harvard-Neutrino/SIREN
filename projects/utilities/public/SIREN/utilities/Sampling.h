#pragma once
#ifndef SIREN_Sampling_H
#define SIREN_Sampling_H

#include <cmath>
#include <memory>
#include <utility>
#include <vector>
#include <cassert>

#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren {
namespace utilities {

/**
* @brief Performs Independent Metropolis-Hastings sampling of a differential
* distribution.
*
* Each call to func_proposal must return a pair (sample, proposal_density),
* where sample is independent of the current chain state.  The density is
* used to form the correct independent-MH acceptance ratio:
*     alpha = min(1, [pi(x') * q(x)] / [pi(x) * q(x')])
* Density only needs to be correct up to a constant that cancels in the ratio,
* so a uniform proposal can return any fixed positive value (e.g. 1.0).
*
* @param func_proposal callable returning std::pair<std::vector<double>, double>
*                       (sample, proposal density)
* @param func_likelihood callable returning the target density of a sample
* @param random random number generator
* @param burnin number of burn-in iterations before returning a sample
*/
template<typename FuncTypeProposal, typename FuncTypeLikelihood>
std::vector<double> MetropolisHasting_Sample(const FuncTypeProposal& func_proposal, const FuncTypeLikelihood& func_likelihood, std::shared_ptr<siren::utilities::SIREN_random> random, const size_t burnin = 40) {

    auto initial = func_proposal();
    std::vector<double> vars = initial.first;
    double q = initial.second;
    double llh = func_likelihood(vars);
    // Importance weight pi/q; ratios of these give the independent-MH alpha.
    double weight = (q > 0) ? llh / q : 0.0;

    for (size_t j = 0; j <= burnin; j++) {
        auto proposed = func_proposal();
        std::vector<double> const & test_vars = proposed.first;
        double test_q = proposed.second;
        double test_llh = func_likelihood(test_vars);
        double test_weight = (test_q > 0) ? test_llh / test_q : 0.0;

        bool accept;
        if (weight <= 0) {
            accept = true;
        } else {
            accept = (test_weight >= weight) || (random->Uniform(0,1) < test_weight / weight);
        }
        if(accept) {
            vars = test_vars;
            llh = test_llh;
            weight = test_weight;
        }
    }
    return vars;
}

}
}

#endif // SIREN_Sampling_H