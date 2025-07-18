#pragma once
#ifndef SIREN_Sampling_H
#define SIREN_Sampling_H

#include <cmath>
#include <vector>
#include <cassert>

#include "SIREN/utilities/Random.h"               // for SIREN_random

namespace siren {
namespace utilities {

/**
* @brief Performs a Metropolic Hasting sampling of a differential distribution
*
*
* @param func the function to integrate
* @param tol the (absolute) tolerance on the error of the integral
*/
template<typename FuncTypeProposal, typename FuncTypeLikelihood>
std::vector<double> MetropolisHasting_Sample(const FuncTypeProposal& func_proposal, const FuncTypeLikelihood& func_likelihood, std::shared_ptr<siren::utilities::SIREN_random> random, const size_t burnin = 40) {

    std::vector<double> vars, test_vars;

    // sample an initial point
    vars = func_proposal();

    double llh, test_llh;
    llh = func_likelihood(vars);

    for (size_t j = 0; j <= burnin; j++) {
        test_vars = func_proposal;
        test_llh = func_likelihood(test_vars);
        double odds = test_llh / llh;
        accept = (llh==0 || odds>1. || random->Uniform(0,1) < odds);
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