#pragma once
#ifndef SIREN_Sampling_H
#define SIREN_Sampling_H

#include <cmath>
#include <memory>
#include <vector>
#include <cassert>

#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren {
namespace utilities {

/**
* @brief Performs a Metropolis-Hastings sampling of a differential distribution
*
* @param func_proposal callable returning a proposed sample (vector<double>)
* @param func_likelihood callable returning the likelihood of a sample
* @param random random number generator
* @param burnin number of burn-in iterations before returning a sample
*/
template<typename FuncTypeProposal, typename FuncTypeLikelihood>
std::vector<double> MetropolisHasting_Sample(const FuncTypeProposal& func_proposal, const FuncTypeLikelihood& func_likelihood, std::shared_ptr<siren::utilities::SIREN_random> random, const size_t burnin = 40) {

    std::vector<double> vars, test_vars;
    bool accept;

    // sample an initial point
    vars = func_proposal();

    double llh, test_llh;
    llh = func_likelihood(vars);

    for (size_t j = 0; j <= burnin; j++) {
        test_vars = func_proposal();
        test_llh = func_likelihood(test_vars);
        if (llh <= 0) {
            accept = true;
        } else {
            accept = (test_llh >= llh) || (random->Uniform(0,1) < test_llh / llh);
        }
        if(accept) {
            vars = test_vars;
            llh = test_llh;
        }
    }
    return vars;
}

}
}

#endif // SIREN_Sampling_H