#pragma once
#ifndef SIREN_InvariantMassMapping_H
#define SIREN_InvariantMassMapping_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "SIREN/math/InterpolationUtils.h"

namespace siren {
namespace injection {

// Abstract 1-D importance map.  A single object BOTH draws a value
// (Forward) and reports its own normalized density (Density) at an
// arbitrary value, so the sampler and the reported density cannot drift
// apart.  The concrete maps below are plain value types; they derive from
// this base only so that a physics model or channel can hold a shared
// mapping polymorphically, and so a future adaptive map has a place to
// hook batch updates.
class Mapping1D {
public:
    virtual ~Mapping1D() = default;
    // Map a uniform deviate r in [0,1] to the variable.
    virtual double Forward(double r) const = 0;
    // Map the variable back to its uniform deviate in [0,1].
    virtual double Inverse(double x) const = 0;
    // Normalized proposal density at the variable value x.
    virtual double Density(double x) const = 0;
    // Optional adaptive hooks: fixed maps ignore them; an adaptive map
    // accumulates sampled (x, weight) pairs and rebins on Refine().
    virtual void Accumulate(double /*x*/, double /*weight*/) {}
    virtual void Refine() {}
};

// Maps a uniform random number r in [0,1] to an invariant mass squared
// value s in [s_min, s_max], concentrating samples near a Breit-Wigner
// resonance at (M, Gamma).
//
// Forward:  s = M*Gamma*tan(u1 + (u2-u1)*r) + M^2
//   where u1 = atan((s_min - M^2) / (M*Gamma))
//         u2 = atan((s_max - M^2) / (M*Gamma))
//
// Density: g(s) = M*Gamma / ((u2-u1) * [(s-M^2)^2 + M^2*Gamma^2])
//
// This is the standard MadGraph Breit-Wigner mapping.
struct BreitWignerMapping : public Mapping1D {
    double M;         // resonance mass
    double Gamma;     // resonance width
    double s_min;     // minimum invariant mass squared
    double s_max;     // maximum invariant mass squared
    double u_min;     // atan((s_min - M^2) / (M*Gamma))
    double u_max;     // atan((s_max - M^2) / (M*Gamma))
    double norm;      // u_max - u_min

    BreitWignerMapping(double M_, double Gamma_, double s_min_, double s_max_)
        : M(M_), Gamma(Gamma_), s_min(s_min_), s_max(s_max_)
    {
        double MG = M * Gamma;
        double M2 = M * M;
        u_min = std::atan((s_min - M2) / MG);
        u_max = std::atan((s_max - M2) / MG);
        norm = u_max - u_min;
    }

    // Map r in [0,1] to s in [s_min, s_max]
    double Forward(double r) const {
        return M * Gamma * std::tan(u_min + norm * r) + M * M;
    }

    // Map s back to r in [0,1]
    double Inverse(double s) const {
        return (std::atan((s - M * M) / (M * Gamma)) - u_min) / norm;
    }

    // Density g(s) = probability density at invariant mass squared s
    double Density(double s) const {
        if (!std::isfinite(s) || s < s_min || s > s_max) return 0.0;
        double ds = s - M * M;
        double MG = M * Gamma;
        return MG / (norm * (ds * ds + MG * MG));
    }
};

// Maps r in [0,1] to s in [s_min, s_max] with a power-law distribution.
// Useful for non-resonant propagators.
//
// Forward:  s = [r * (s_max - m^2)^(1-nu) + (1-r) * (s_min - m^2)^(1-nu)]^(1/(1-nu)) + m^2
//
// For nu=1 (logarithmic):
//   s = (s_min - m^2)^(1-r) * (s_max - m^2)^r + m^2
//
// m^2 is an offset (typically 0 or a threshold mass squared).
struct PowerLawMapping : public Mapping1D {
    double nu;        // power law exponent (typically 0.8)
    double m2;        // mass offset squared
    double s_min;
    double s_max;
    double a_min;     // (s_min - m2)^(1-nu)
    double a_max;     // (s_max - m2)^(1-nu)
    double one_minus_nu;
    double inv_one_minus_nu;
    bool logarithmic;
    double log_min;   // log(s_min - m2), for nu == 1
    double log_range; // log((s_max - m2) / (s_min - m2)), for nu == 1

    PowerLawMapping(double nu_, double m2_, double s_min_, double s_max_)
        : nu(nu_), m2(m2_), s_min(s_min_), s_max(s_max_)
    {
        if (!std::isfinite(nu) || !std::isfinite(m2) ||
            !std::isfinite(s_min) || !std::isfinite(s_max) ||
            !(s_max > s_min) || s_min < m2) {
            throw std::runtime_error(
                "PowerLawMapping requires finite parameters with "
                "m2 <= s_min < s_max");
        }

        one_minus_nu = 1.0 - nu;
        logarithmic = (nu == 1.0);
        if (logarithmic) {
            double x_min = s_min - m2;
            double x_max = s_max - m2;
            if (!(x_min > 0.0)) {
                throw std::runtime_error(
                    "PowerLawMapping with nu == 1 requires m2 < s_min");
            }
            log_min = std::log(x_min);
            log_range = std::log(x_max / x_min);
            if (!(log_range > 0.0) || !std::isfinite(log_range)) {
                throw std::runtime_error(
                    "PowerLawMapping logarithmic range is degenerate");
            }
            inv_one_minus_nu = 0.0;
            a_min = 0.0;
            a_max = 0.0;
            return;
        }

        // At x_min = 0, x^-nu is normalizable only for nu < 1.
        if (s_min == m2 && one_minus_nu <= 0.0) {
            throw std::runtime_error(
                "PowerLawMapping requires m2 < s_min when nu >= 1");
        }

        inv_one_minus_nu = 1.0 / one_minus_nu;
        a_min = std::pow(s_min - m2, one_minus_nu);
        a_max = std::pow(s_max - m2, one_minus_nu);
        log_min = 0.0;
        log_range = 0.0;
        if (!std::isfinite(a_min) || !std::isfinite(a_max) ||
            a_max == a_min) {
            throw std::runtime_error(
                "PowerLawMapping has a non-finite or degenerate normalization");
        }
    }

    double Forward(double r) const {
        double s;
        if (logarithmic) {
            s = std::exp(log_min + r * log_range) + m2;
        } else {
            s = std::pow(
                r * a_max + (1.0 - r) * a_min,
                inv_one_minus_nu) + m2;
        }
        return std::clamp(s, s_min, s_max);
    }

    double Inverse(double s) const {
        if (s <= s_min) return 0.0;
        if (s >= s_max) return 1.0;
        if (logarithmic) {
            return (std::log(s - m2) - log_min) / log_range;
        }
        double a = std::pow(s - m2, one_minus_nu);
        return (a - a_min) / (a_max - a_min);
    }

    double Density(double s) const {
        if (!std::isfinite(s) || s < s_min || s > s_max) return 0.0;
        double ds = s - m2;
        if (logarithmic) return 1.0 / (ds * log_range);
        return std::pow(ds, -nu) * one_minus_nu / (a_max - a_min);
    }
};

// Maps r in [0,1] to s in [s_min, s_max] by inverting a tabulated CDF.
//
// The caller supplies monotone-increasing s_nodes and the corresponding
// cumulative distribution cdf_nodes (cdf_nodes[0] need not be 0 nor
// cdf_nodes[back] be 1 -- the table is renormalized internally to the
// overlap with [s_min, s_max]).  Between nodes the CDF is treated as
// piecewise-linear, so the density is piecewise-constant per bin.  This
// is the general, near-exact importance map: build the table once from a
// model's own differential rate / cross section and the proposal matches
// the physical marginal for any masses (and is coupling-independent,
// since the normalization cancels).
//
// Forward:  invert the (clipped, renormalized) CDF for s.
// Density:  per-bin slope d(cdf)/ds, renormalized over [s_min, s_max].
struct TabulatedMappingTable {
    std::vector<double> s;        // node positions (ascending)
    std::vector<double> cdf;      // cumulative weight at each node (ascending)

    TabulatedMappingTable(
        std::vector<double> s_nodes,
        std::vector<double> cdf_nodes)
        : s(std::move(s_nodes)), cdf(std::move(cdf_nodes))
    {
        if (s.size() < 2 || s.size() != cdf.size()) {
            throw std::runtime_error(
                "TabulatedMapping requires matching s/cdf arrays of length >= 2");
        }
        for (std::size_t i = 0; i < s.size(); ++i) {
            if (!std::isfinite(s[i]) ||
                (i > 0 && !(s[i] > s[i - 1]))) {
                throw std::runtime_error(
                    "TabulatedMapping s nodes must be finite and strictly increasing");
            }
            if (!std::isfinite(cdf[i]) ||
                (i > 0 && cdf[i] < cdf[i - 1])) {
                throw std::runtime_error(
                    "TabulatedMapping CDF nodes must be finite and nondecreasing");
            }
        }
    }

    double InterpCdf(double sv) const {
        if (sv <= s.front()) return cdf.front();
        if (sv >= s.back()) return cdf.back();
        auto [lo, hi] = siren::math::InterpolationBracket(s, sv);
        double ds = s[hi] - s[lo];
        double frac = ds > 0.0 ? (sv - s[lo]) / ds : 0.0;
        return cdf[lo] + frac * (cdf[hi] - cdf[lo]);
    }

    double SlopeAt(double sv) const {
        if (sv < s.front() || sv > s.back()) return 0.0;
        auto [lo, hi] = siren::math::InterpolationBracket(s, sv);
        double ds = s[hi] - s[lo];
        if (ds <= 0.0) return 0.0;
        return (cdf[hi] - cdf[lo]) / ds;
    }
};

struct TabulatedMapping : public Mapping1D {
    std::shared_ptr<TabulatedMappingTable const> table;
    double s_min;
    double s_max;
    double c_lo;                  // interpolated cumulative at s_min
    double c_hi;                  // interpolated cumulative at s_max
    double norm;                  // c_hi - c_lo (cumulative mass in range)

    TabulatedMapping(std::vector<double> s_nodes,
                     std::vector<double> cdf_nodes,
                     double s_min_, double s_max_)
        : TabulatedMapping(
            std::make_shared<TabulatedMappingTable>(
                std::move(s_nodes), std::move(cdf_nodes)),
            s_min_, s_max_)
    {}

    TabulatedMapping(
        std::shared_ptr<TabulatedMappingTable const> table_,
        double s_min_, double s_max_)
        : table(std::move(table_)), s_min(s_min_), s_max(s_max_)
    {
        if (!table) throw std::runtime_error("TabulatedMapping requires a table");
        c_lo = InterpCdf(s_min);
        c_hi = InterpCdf(s_max);
        norm = c_hi - c_lo;
    }

    // A per-event kinematic window may lie entirely outside the immutable
    // table, or overlap only a zero-weight plateau. That is an inactive
    // proposal for this event, not a malformed table: Density must remain
    // evaluable so another channel in the mixture can supply the point.
    bool HasSupport() const {
        return std::isfinite(norm) && norm > 0.0;
    }

    // Linear-interpolate the stored cumulative at an arbitrary s (clamped).
    double InterpCdf(double sv) const {
        return table->InterpCdf(sv);
    }

    // Local density d(cdf)/ds at sv (bin slope), unnormalized.
    double SlopeAt(double sv) const {
        return table->SlopeAt(sv);
    }

    double Forward(double r) const {
        if (!HasSupport()) {
            throw std::runtime_error(
                "TabulatedMapping cannot sample a window with no cumulative support");
        }
        double target = c_lo + r * norm;             // cumulative we want
        // Find the bin [s[lo], s[hi]] whose cumulative brackets target.
        if (target <= table->cdf.front()) return s_min;
        if (target >= table->cdf.back()) return s_max;
        auto [lo, hi] = siren::math::InterpolationBracket(
            table->cdf, target);
        double dc = table->cdf[hi] - table->cdf[lo];
        double frac = dc > 0.0 ? (target - table->cdf[lo]) / dc : 0.0;
        double sv = table->s[lo] + frac * (table->s[hi] - table->s[lo]);
        // Guard against rounding outside the channel's allowed range.
        if (sv < s_min) sv = s_min;
        if (sv > s_max) sv = s_max;
        return sv;
    }

    double Inverse(double sv) const {
        if (!HasSupport()) return 0.0;
        if (sv <= s_min) return 0.0;
        if (sv >= s_max) return 1.0;
        return (InterpCdf(sv) - c_lo) / norm;
    }

    // Normalized proposal density in s over [s_min, s_max].
    double Density(double sv) const {
        if (!HasSupport()) return 0.0;
        if (sv < s_min || sv > s_max) return 0.0;
        return SlopeAt(sv) / norm;
    }
};

// Maps r in [0,1] to x in [x_min, x_max] following a t-channel
// (spacelike) propagator density g(x) ~ 1 / (x + m2)^2.
//
// This is the natural importance map for the momentum transfer of a
// t-channel exchange: x = Q^2 and m2 = m_mediator^2.  The dominant
// shape of dsigma/dQ^2 ~ 1/(Q^2 + m_V^2)^2 (the propagator squared) is
// captured exactly; the residual matrix-element numerator and nuclear
// form factor are smooth and slowly varying, so sampling from this map
// leaves only mild residual weight variance.
//
// Density:  g(x) = [1 / (x + m2)^2] / C,
//   with C = 1/(x_min + m2) - 1/(x_max + m2)   (so g integrates to 1).
//
// Forward (inverse CDF):
//   x = 1 / ( 1/(x_min + m2) - r * C ) - m2
//
// Inverse:
//   r = ( 1/(x_min + m2) - 1/(x + m2) ) / C
//
// The closed form requires the pole at x = -m2 to lie outside the
// closed window: (x_min + m2) and (x_max + m2) must share a sign.  The
// density is non-normalizable across the pole, so a window containing
// it is rejected at construction; callers whose kinematic range spans
// the pole must restrict the window to one side first.  Windows
// entirely below the pole (both offsets negative) are valid: the
// normalization C = 1/(x_min + m2) - 1/(x_max + m2) is positive on
// either side.
struct PropagatorMapping : public Mapping1D {
    double m2;        // mediator mass squared (propagator offset)
    double x_min;
    double x_max;
    double inv_lo;    // 1 / (x_min + m2)
    double inv_hi;    // 1 / (x_max + m2)
    double C;         // inv_lo - inv_hi (normalization; > 0 for x_max > x_min)

    PropagatorMapping(double m2_, double x_min_, double x_max_)
        : m2(m2_), x_min(x_min_), x_max(x_max_)
    {
        double d_lo = x_min + m2;
        double d_hi = x_max + m2;
        bool same_side = (d_lo > 0.0 && d_hi > 0.0) ||
                         (d_lo < 0.0 && d_hi < 0.0);
        if (!std::isfinite(d_lo) || !std::isfinite(d_hi) || !same_side) {
            throw std::invalid_argument(
                "PropagatorMapping window may not contain the pole at "
                "x = -m2; restrict the window to one side of the pole");
        }
        inv_lo = 1.0 / d_lo;
        inv_hi = 1.0 / d_hi;
        C = inv_lo - inv_hi;
    }

    double Forward(double r) const {
        // Guard the degenerate (zero-width) range.
        if (!(C > 0.0)) return x_min;
        double denom = inv_lo - r * C;
        // Above the pole denom stays positive up to subtraction
        // rounding at r -> 1; below the pole it is negative throughout
        // and the closed form is exact.
        if (inv_lo > 0.0 && denom <= 0.0) return x_max;
        double x = 1.0 / denom - m2;
        return std::clamp(x, x_min, x_max);
    }

    double Inverse(double x) const {
        if (!(C > 0.0)) return 0.0;
        return (inv_lo - 1.0 / (x + m2)) / C;
    }

    double Density(double x) const {
        if (!std::isfinite(x) || x < x_min || x > x_max) return 0.0;
        if (!(C > 0.0)) return 0.0;
        double d = x + m2;
        return 1.0 / (d * d * C);
    }
};

// Uniform mapping: s = s_min + r * (s_max - s_min)
struct UniformMapping : public Mapping1D {
    double s_min;
    double s_max;
    double range;

    UniformMapping(double s_min_, double s_max_)
        : s_min(s_min_), s_max(s_max_), range(s_max_ - s_min_) {}

    double Forward(double r) const { return s_min + r * range; }
    double Inverse(double s) const { return (s - s_min) / range; }
    double Density(double s) const {
        if (!std::isfinite(s) || s < s_min || s > s_max) return 0.0;
        return 1.0 / range;
    }
};

// Log (scale-free) mapping: g(x) ~ 1/x on [x_min, x_max], x_min > 0.
// This is the offset-free equivalent of PowerLawMapping's nu == 1 branch,
// exposed separately for callers that want an explicitly logarithmic map.
// It is the natural map for a massless/scale-free propagator (pure 1/Q^2)
// or any variable swept over decades.
//   CDF:     F(x) = ln(x/x_min) / L,   L = ln(x_max/x_min)
//   Forward: x = x_min * exp(r * L)
//   Density: g(x) = 1 / (x * L)
struct LogMapping : public Mapping1D {
    double x_min;
    double x_max;
    double L;          // ln(x_max / x_min)

    LogMapping(double x_min_, double x_max_)
        : x_min(x_min_), x_max(x_max_)
    {
        if (!(x_min > 0.0) || !(x_max > x_min)) {
            throw std::runtime_error("LogMapping requires 0 < x_min < x_max");
        }
        L = std::log(x_max / x_min);
    }

    double Forward(double r) const { return x_min * std::exp(r * L); }
    double Inverse(double x) const { return std::log(x / x_min) / L; }
    double Density(double x) const {
        if (!std::isfinite(x) || x < x_min || x > x_max) return 0.0;
        return 1.0 / (x * L);
    }
};

// Exponential mapping: g(x) ~ exp(-x/tau) on [x_min, x_max], tau != 0.
// The natural importance map for a decay length / proper time (the
// exponential decay law) and for any falling (tau > 0) or rising (tau < 0)
// exponential variable.
//   a = exp(-x_min/tau), b = exp(-x_max/tau), norm = tau*(a-b) (> 0 either sign)
//   Forward: x = -tau * ln(a - r*(a-b))
//   Density: g(x) = exp(-x/tau) / norm
struct ExponentialMapping : public Mapping1D {
    double tau;
    double x_min;
    double x_max;
    double a;          // exp(-x_min/tau)
    double b;          // exp(-x_max/tau)
    double norm;       // tau*(a-b) = integral of exp(-x/tau) over the range

    ExponentialMapping(double tau_, double x_min_, double x_max_)
        : tau(tau_), x_min(x_min_), x_max(x_max_)
    {
        if (tau == 0.0 || !(x_max > x_min)) {
            throw std::runtime_error(
                "ExponentialMapping requires tau != 0 and x_max > x_min");
        }
        a = std::exp(-x_min / tau);
        b = std::exp(-x_max / tau);
        norm = tau * (a - b);
        if (!(norm > 0.0)) {
            throw std::runtime_error("ExponentialMapping has degenerate range");
        }
    }

    double Forward(double r) const { return -tau * std::log(a - r * (a - b)); }
    double Inverse(double x) const { return (a - std::exp(-x / tau)) / (a - b); }
    double Density(double x) const {
        if (!std::isfinite(x) || x < x_min || x > x_max) return 0.0;
        return std::exp(-x / tau) / norm;
    }
};

namespace detail {

// Standard normal CDF, Phi(z) = 0.5*erfc(-z/sqrt(2)).
inline double NormalCdf(double z) {
    return 0.5 * std::erfc(-z * 0.70710678118654752440);  // 1/sqrt(2)
}

// Inverse standard normal CDF (probit).  Acklam's rational approximation,
// refined with one Halley step using erfc, reaches ~machine precision -- so a
// Gaussian map's Forward exactly inverts the CDF its Density reports (C1).
inline double NormalQuantile(double p) {
    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return std::numeric_limits<double>::infinity();
    static const double a[6] = {
        -3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00};
    static const double b[5] = {
        -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[6] = {
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[4] = {
         7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
         3.754408661907416e+00};
    const double p_low = 0.02425, p_high = 1.0 - p_low;
    double x;
    if (p < p_low) {
        double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (p <= p_high) {
        double q = p - 0.5, r = q * q;
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    // One Halley refinement: drives the error to ~machine precision.
    double e = NormalCdf(x) - p;
    double u = e * 2.50662827463100050242 * std::exp(0.5 * x * x);  // sqrt(2*pi)
    x -= u / (1.0 + 0.5 * x * u);
    return x;
}

} // namespace detail

// Gaussian mapping: g(x) ~ exp(-(x-mu)^2/(2 sigma^2)) truncated to
// [x_min, x_max] (sigma > 0).  Useful for a smooth peaked variable (beam
// energy spread, resolution smearing).  Density is the exact truncated-normal
// pdf; Forward uses the probit refined to ~machine precision, so it inverts
// exactly the CDF Density reports (C1).
struct GaussianMapping : public Mapping1D {
    double mu;
    double sigma;
    double x_min;
    double x_max;
    double p_lo;       // Phi((x_min - mu)/sigma)
    double p_hi;       // Phi((x_max - mu)/sigma)
    double norm;       // p_hi - p_lo

    GaussianMapping(double mu_, double sigma_, double x_min_, double x_max_)
        : mu(mu_), sigma(sigma_), x_min(x_min_), x_max(x_max_)
    {
        if (!(sigma > 0.0) || !(x_max > x_min)) {
            throw std::runtime_error(
                "GaussianMapping requires sigma > 0 and x_max > x_min");
        }
        p_lo = detail::NormalCdf((x_min - mu) / sigma);
        p_hi = detail::NormalCdf((x_max - mu) / sigma);
        norm = p_hi - p_lo;
        if (!(norm > 0.0)) {
            throw std::runtime_error(
                "GaussianMapping has zero probability mass over [x_min, x_max]");
        }
    }

    double Forward(double r) const {
        double x = mu + sigma * detail::NormalQuantile(p_lo + r * norm);
        if (x < x_min) x = x_min;
        if (x > x_max) x = x_max;
        return x;
    }

    double Inverse(double x) const {
        return (detail::NormalCdf((x - mu) / sigma) - p_lo) / norm;
    }

    double Density(double x) const {
        if (!std::isfinite(x) || x < x_min || x > x_max) return 0.0;
        double z = (x - mu) / sigma;
        // 1/sqrt(2*pi) = 0.39894228040143267794
        return (0.39894228040143267794 / sigma) * std::exp(-0.5 * z * z) / norm;
    }
};

// Adaptive (VEGAS-style) mapping: a piecewise-constant proposal on a fixed
// grid of n_bins equal-width bins over [x_min, x_max], whose per-bin
// probabilities self-tune from sampled weights via the Accumulate/Refine
// hooks.  It starts uniform (so it is a no-op until refined) and converges
// toward any 1-D shape with no hand-tuned parameters -- the foundation for
// adaptive intra-channel sampling.
//
// Contract: feed Accumulate(x, w) with the importance weight w = f(x)/g(x)
// for points x drawn from THIS map; then sum over a batch of bin i estimates
// the bin integral of f (the g-sampling cancels), so Refine sets the new
// per-bin probability proportional to that, damping-blended with the current
// one and floored so no bin dies.  The result drives the proposal toward f.
struct AdaptiveMapping : public Mapping1D {
    double x_min;
    double x_max;
    int n_bins;
    double damping;            // blend new/old per-bin probability in [0,1]
    double floor_frac;         // min bin probability = floor_frac / n_bins
    double width;              // (x_max - x_min) / n_bins
    std::vector<double> p;     // per-bin probability (sums to 1)
    std::vector<double> cum;   // cumulative probability at bin edges (n_bins+1)
    std::vector<double> acc;   // accumulated weight per bin (reset on Refine)
    std::vector<long> cnt;     // sample count per bin (reset on Refine)

    AdaptiveMapping(double x_min_, double x_max_, int n_bins_ = 32,
                    double damping_ = 0.5, double floor_frac_ = 1e-3)
        : x_min(x_min_), x_max(x_max_), n_bins(n_bins_),
          damping(damping_), floor_frac(floor_frac_)
    {
        if (n_bins < 1 || !(x_max > x_min)) {
            throw std::runtime_error(
                "AdaptiveMapping requires n_bins >= 1 and x_max > x_min");
        }
        width = (x_max - x_min) / n_bins;
        p.assign(n_bins, 1.0 / n_bins);
        acc.assign(n_bins, 0.0);
        cnt.assign(n_bins, 0);
        cum.assign(n_bins + 1, 0.0);
        RebuildCum();
    }

    void RebuildCum() {
        cum[0] = 0.0;
        for (int i = 0; i < n_bins; ++i) cum[i + 1] = cum[i] + p[i];
        double tot = cum[n_bins];
        if (tot > 0.0) for (double & c : cum) c /= tot;
    }

    int BinOf(double x) const {
        if (x <= x_min) return 0;
        if (x >= x_max) return n_bins - 1;
        int i = static_cast<int>((x - x_min) / width);
        if (i < 0) i = 0;
        if (i >= n_bins) i = n_bins - 1;
        return i;
    }

    double Forward(double r) const {
        if (r <= 0.0) return x_min;
        if (r >= 1.0) return x_max;
        int lo = 0, hi = n_bins;  // search the n_bins+1 edge cumulatives
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (cum[mid] <= r) lo = mid; else hi = mid;
        }
        double pi = cum[lo + 1] - cum[lo];
        double frac = pi > 0.0 ? (r - cum[lo]) / pi : 0.0;
        return x_min + (lo + frac) * width;
    }

    double Inverse(double x) const {
        if (x <= x_min) return 0.0;
        if (x >= x_max) return 1.0;
        int i = BinOf(x);
        double frac = (x - (x_min + i * width)) / width;
        return cum[i] + frac * (cum[i + 1] - cum[i]);
    }

    double Density(double x) const {
        if (x < x_min || x > x_max) return 0.0;
        return p[BinOf(x)] / width;
    }

    void Accumulate(double x, double weight) override {
        if (x < x_min || x > x_max || !std::isfinite(weight)) return;
        int i = BinOf(x);
        acc[i] += weight;
        cnt[i] += 1;
    }

    void Refine() override {
        double tot = 0.0;
        for (double v : acc) tot += v;
        if (tot > 0.0) {
            std::vector<double> np(n_bins);
            double min_p = floor_frac / n_bins;
            double s = 0.0;
            for (int i = 0; i < n_bins; ++i) {
                np[i] = std::max(acc[i] / tot, min_p);
                s += np[i];
            }
            for (int i = 0; i < n_bins; ++i) np[i] /= s;
            double ps = 0.0;
            for (int i = 0; i < n_bins; ++i) {
                p[i] = damping * np[i] + (1.0 - damping) * p[i];
                ps += p[i];
            }
            for (double & v : p) v /= ps;
            RebuildCum();
        }
        std::fill(acc.begin(), acc.end(), 0.0);
        std::fill(cnt.begin(), cnt.end(), 0);
    }
};

} // namespace injection
} // namespace siren

#endif // SIREN_InvariantMassMapping_H
