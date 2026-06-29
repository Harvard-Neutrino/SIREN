#pragma once
#ifndef SIREN_CharmDecayTestHelpers_H
#define SIREN_CharmDecayTestHelpers_H

// Shared numeric oracle + record helpers for the CharmMesonDecay closure tests.
// Both decay test files include this so the V-A angle-average cross-check and
// the q^2 reconstruction stay in one place.

#include <cmath>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "SIREN/dataclasses/InteractionRecord.h"

namespace siren {
namespace interactions {
namespace charm_decay_test {

// q^2 = (p_D - p_K)^2 from a finalized record.
inline double reconstruct_q2(const siren::dataclasses::InteractionRecord & rec) {
    double qE  = rec.primary_momentum[0] - rec.secondary_momenta[0][0];
    double qpx = rec.primary_momentum[1] - rec.secondary_momenta[0][1];
    double qpy = rec.primary_momentum[2] - rec.secondary_momenta[0][2];
    double qpz = rec.primary_momentum[3] - rec.secondary_momenta[0][3];
    return qE * qE - qpx * qpx - qpy * qpy - qpz * qpz;
}

// Independent numeric quadrature of the clamped V-A weight, evaluated with the
// same rk::P4 boosts SampleFinalState uses. Cross-checks the closed-form
// charm_decay::VAWeightAngleAverage used by the weighting code.
inline double numericVAWeightAngleAverage(double mD, double mK, double ml, double m23) {
    double mnu = 0.0;
    double p1Abs = 0.5 * std::sqrt((mD - mK - m23) * (mD + mK + m23)
                                 * (mD + mK - m23) * (mD - mK + m23)) / mD;
    double p23Abs = 0.5 * std::sqrt((m23 - ml - mnu) * (m23 + ml + mnu)
                                  * (m23 + ml - mnu) * (m23 - ml + mnu)) / m23;
    if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
    double wtMEmax = std::min(std::pow(mD, 4) / 16.0,
                              mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));
    rk::P4 p4m23(geom3::Vector3(0.0, 0.0, -p1Abs), m23);
    rk::Boost boost = p4m23.labBoost();
    rk::P4 p4K(geom3::Vector3(0.0, 0.0, p1Abs), mK);
    const int N = 20000;            // even -> composite Simpson
    double h = 2.0 / N, sum = 0.0;
    for (int i = 0; i <= N; ++i) {
        double c = -1.0 + i * h;
        double s = std::sqrt(std::max(0.0, 1.0 - c * c));
        geom3::Vector3 dir(s, 0.0, c);
        rk::P4 p4l = rk::P4(p23Abs * dir, ml).boost(boost);
        rk::P4 p4nu = rk::P4(-p23Abs * dir, mnu).boost(boost);
        double w = mD * p4l.e() * p4nu.dot(p4K);
        if (w < 0.0) w = 0.0;
        if (w > wtMEmax) w = wtMEmax;
        double wgt = (i == 0 || i == N) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += wgt * w;
    }
    return 0.5 * (sum * h / 3.0);
}

} // namespace charm_decay_test
} // namespace interactions
} // namespace siren

#endif // SIREN_CharmDecayTestHelpers_H
