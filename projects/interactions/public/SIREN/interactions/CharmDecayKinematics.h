#pragma once
#ifndef SIREN_CharmDecayKinematics_H
#define SIREN_CharmDecayKinematics_H

// Shared closure kinematics for the charm-meson semileptonic decay samplers.
//
// CharmMesonDecay and CharmMesonDecay3Body are independent classes (no shared
// base), but their SampleFinalState <-> FinalStateProbability closure relies on
// the same q^2 density. These inline free helpers are the single source of that
// math so both classes stay byte-for-byte consistent.
//
// Closure rationale (authoritative copy): FinalStateProbability must be the
// normalized q^2 density that SampleFinalState produces, because the Weighter
// consumes it as the physical final-state density. The sampler draws
// m23 = sqrt(q^2) flat with accept-reject on the phase-space weight
// p1Abs*p23Abs, then (for D+/D0) accept-rejects on the V-A weight
// wtME = mD*E_l*(p_nu . p_K). The accepted m23 density is therefore proportional
// to p1Abs(m23)*p23Abs(m23)*<max(0,wtME)>_angle, and the q^2 density carries the
// Jacobian dm23/dq^2 = 1/(2 m23). SampledQ2Density reproduces exactly this; the
// per-event numerator and the normalization integral both call it, so
// Sample == Density by construction. The form-factor DifferentialDecayWidth is a
// separate physical quantity the sampler never used and is NOT the closure
// density.

#include <cmath>
#include <map>
#include <functional>

#include "SIREN/dataclasses/Particle.h"
#include "SIREN/math/Kinematics.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Integration.h"

namespace siren {
namespace interactions {
namespace charm_decay {

inline double particleMass(siren::dataclasses::ParticleType particle) {
    switch(particle){
			case siren::dataclasses::ParticleType::D0:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::D0Bar:
				return( siren::utilities::Constants::D0Mass);
			case siren::dataclasses::ParticleType::DPlus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::DMinus:
				return( siren::utilities::Constants::DPlusMass);
			case siren::dataclasses::ParticleType::K0:
				return( siren::utilities::Constants::K0Mass);
			case siren::dataclasses::ParticleType::K0Bar:
				return( siren::utilities::Constants::K0Mass);
			case siren::dataclasses::ParticleType::KPlus:
				return( siren::utilities::Constants::KPlusMass);
			case siren::dataclasses::ParticleType::KMinus:
				return( siren::utilities::Constants::KMinusMass);
      case siren::dataclasses::ParticleType::EPlus:
        return( siren::utilities::Constants::electronMass );
      case siren::dataclasses::ParticleType::EMinus:
        return( siren::utilities::Constants::electronMass );
      case siren::dataclasses::ParticleType::MuPlus:
        return( siren::utilities::Constants::muonMass );
      case siren::dataclasses::ParticleType::MuMinus:
        return( siren::utilities::Constants::muonMass );
      case siren::dataclasses::ParticleType::TauPlus:
        return( siren::utilities::Constants::tauMass );
      case siren::dataclasses::ParticleType::TauMinus:
        return( siren::utilities::Constants::tauMass );
      case siren::dataclasses::ParticleType::DsPlus:
        return( siren::utilities::Constants::DsPlusMass );
      case siren::dataclasses::ParticleType::DsMinus:
        return( siren::utilities::Constants::DsMinusMass );
      default:
        return(0.0);
    }
}

// Single source of truth for the K*(892) mass shared by the sampler and the
// FinalStateProbability normalizer.
inline double KStarMass() {
  return siren::utilities::Constants::KPrimePlusMass;
}

// Integral over [lo, hi] of max(0, a2*c^2 + a1*c + a0). Closed form with a
// linear/constant fallback for a2 ~ 0 and numerically stable quadratic roots.
// Used to evaluate the angle-averaged V-A weight analytically.
inline double integratePositivePartQuadratic(double a2, double a1, double a0,
                                      double lo, double hi) {
    if (hi <= lo) return 0.0;
    auto F = [&](double c) { return a2 * c * c * c / 3.0 + a1 * c * c / 2.0 + a0 * c; };
    auto seg = [&](double x0, double x1) -> double {
        if (x1 <= x0) return 0.0;
        return F(x1) - F(x0);
    };
    double scale = std::abs(a1) + std::abs(a0) + 1.0;
    if (std::abs(a2) < 1e-12 * scale) {
        // Linear q(c) = a1*c + a0 (or constant if a1 ~ 0).
        if (std::abs(a1) < 1e-300) return (a0 > 0.0) ? seg(lo, hi) : 0.0;
        double root = -a0 / a1;
        if (a1 > 0.0) return seg(std::max(lo, root), hi);   // q > 0 for c > root
        return seg(lo, std::min(hi, root));                 // q > 0 for c < root
    }
    double disc = a1 * a1 - 4.0 * a2 * a0;
    if (disc <= 0.0) {
        // No real roots: q keeps the sign of a2 everywhere.
        return (a2 > 0.0) ? seg(lo, hi) : 0.0;
    }
    double sq = std::sqrt(disc);
    double qq = -0.5 * (a1 + ((a1 >= 0.0) ? sq : -sq));     // stable roots
    double r1 = qq / a2;
    double r2 = a0 / qq;
    double rlo = std::min(r1, r2), rhi = std::max(r1, r2);
    if (a2 > 0.0) {
        // q > 0 outside [rlo, rhi].
        return seg(lo, std::min(hi, rlo)) + seg(std::max(lo, rhi), hi);
    }
    // q > 0 inside (rlo, rhi).
    return seg(std::max(lo, rlo), std::min(hi, rhi));
}

// Analytic angle-average of the sampler's ACCEPTED V-A weight,
//   (1/2) * integral_{-1}^{1} clamp(wtME(c), 0, wtMEmax) dc,
// with c = cos(theta_lepton) in the (l,nu) rest frame. After the boost to the D
// rest frame wtME = pref*(a0 + a1 c + a2 c^2) is an exact quadratic, so
// clamp(q,0,M) = max(0,q) - max(0,q-M) gives two closed-form positive-part
// integrals (no quadrature error). The unit tests cross-check this against a
// numeric quadrature of the identical clamped weight.
inline double VAWeightAngleAverage(double mD, double mK, double ml, double m23) {
    double mnu = 0.0;
    double p1Abs = siren::math::TwoBodyRestMomentum(mD, mK, m23);
    double p23Abs = siren::math::TwoBodyRestMomentum(m23, ml, mnu);
    if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
    double E23 = std::sqrt(p1Abs * p1Abs + m23 * m23);
    double bz = -p1Abs / E23;              // (l,nu)-system velocity along -kaon axis
    double gamma = E23 / m23;
    double Elrest = std::sqrt(p23Abs * p23Abs + ml * ml);
    double EK = std::sqrt(p1Abs * p1Abs + mK * mK);
    // E_l = gamma (C + D c);  (p_nu . p_K) = gamma p23Abs (A + B c).
    double C = Elrest;
    double D = bz * p23Abs;
    double A = EK - p1Abs * bz;
    double B = p1Abs - EK * bz;
    double pref = mD * gamma * gamma * p23Abs;   // > 0
    if (pref <= 0.0) return 0.0;
    // wtME = pref * q(c), q(c) = (C + D c)(A + B c) = a0 + a1 c + a2 c^2.
    double a0 = C * A;
    double a1 = C * B + D * A;
    double a2 = D * B;
    // Same matrix-element ceiling as SampleFinalState (meMode == 22).
    double wtMEmax = std::min(std::pow(mD, 4) / 16.0,
                              mD * (mD - mK - ml) * (mD - mK - mnu) * (mD - ml - mnu));
    double qmax = wtMEmax / pref;
    double integral = integratePositivePartQuadratic(a2, a1, a0, -1.0, 1.0)
                    - integratePositivePartQuadratic(a2, a1, a0 - qmax, -1.0, 1.0);
    return 0.5 * pref * integral;          // average over c (measure width 2)
}

// Unnormalized sampler density in q^2 for one hadron-mass component. apply_va
// toggles the V-A angle-average (D+/D0) vs pure phase space (Ds).
inline double SampledQ2Density(double mD, double mK, double ml, double q2, bool apply_va) {
  double m23 = std::sqrt(q2);
  double m23Max = mD - mK;
  if (m23 <= ml || m23 >= m23Max) return 0.0;
  double mnu = 0.0;
  double p1Abs = siren::math::TwoBodyRestMomentum(mD, mK, m23);
  double p23Abs = siren::math::TwoBodyRestMomentum(m23, ml, mnu);
  if (p1Abs <= 0.0 || p23Abs <= 0.0) return 0.0;
  // Phase-space weight WITH the sampler's rejection ceiling: where
  // p1Abs*p23Abs exceeds wtPSmax the sampler saturates (always accepts), so the
  // accepted density is flat at wtPSmax. Reproduce the clip for exact closure.
  double m23Min = ml + mnu;
  double p1Max = siren::math::TwoBodyRestMomentum(mD, mK, m23Min);
  double p23Max = siren::math::TwoBodyRestMomentum(m23Max, ml, mnu);
  double wtPSmax = 0.5 * p1Max * p23Max;
  double wtPS = p1Abs * p23Abs;
  if (wtPS > wtPSmax) wtPS = wtPSmax;
  double me = apply_va ? VAWeightAngleAverage(mD, mK, ml, m23) : 1.0;
  // Jacobian dm23/dq^2 = 1/(2 m23) converts the m23 density to a q^2 density.
  return wtPS * me / (2.0 * m23);
}

// Integral of SampledQ2Density over the allowed q^2 range, normalizing each
// mixture component to a proper pdf. Cached in norm_cache per (mD, mK, ml,
// apply_va) so the per-event FinalStateProbability call does not re-integrate.
inline double SampledQ2Normalization(double mD, double mK, double ml, bool apply_va,
                                     std::map<long, double> & norm_cache) {
  double key = ((mD * 1e3 + mK) * 1e3 + ml) * 2.0 + (apply_va ? 1.0 : 0.0);
  long ikey = (long)std::llround(key * 1e3);
  auto it = norm_cache.find(ikey);
  if (it != norm_cache.end()) return it->second;
  double q2min = ml * ml;
  double q2max = (mD - mK) * (mD - mK);
  std::function<double(double)> integrand = [&](double q2) -> double {
    return SampledQ2Density(mD, mK, ml, q2, apply_va);
  };
  double n = siren::utilities::rombergIntegrate(integrand, q2min, q2max);
  norm_cache[ikey] = n;
  return n;
}

} // namespace charm_decay
} // namespace interactions
} // namespace siren

#endif // SIREN_CharmDecayKinematics_H
