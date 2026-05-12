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
* @brief Independent Metropolis-Hastings sampler.
*
* @param func_proposal returns pair(sample, proposal_density); density up to a constant
* @param func_likelihood target density at a sample point
* @param random RNG
* @param burnin iterations before returning
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