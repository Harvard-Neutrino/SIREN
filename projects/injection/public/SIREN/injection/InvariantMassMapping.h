#pragma once
#ifndef SIREN_InvariantMassMapping_H
#define SIREN_InvariantMassMapping_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace siren {
namespace injection {

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
struct BreitWignerMapping {
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
struct PowerLawMapping {
    double nu;        // power law exponent (typically 0.8)
    double m2;        // mass offset squared
    double s_min;
    double s_max;
    double a_min;     // (s_min - m2)^(1-nu)
    double a_max;     // (s_max - m2)^(1-nu)
    double one_minus_nu;
    double inv_one_minus_nu;

    PowerLawMapping(double nu_, double m2_, double s_min_, double s_max_)
        : nu(nu_), m2(m2_), s_min(s_min_), s_max(s_max_)
    {
        one_minus_nu = 1.0 - nu;
        inv_one_minus_nu = 1.0 / one_minus_nu;
        a_min = std::pow(s_min - m2, one_minus_nu);
        a_max = std::pow(s_max - m2, one_minus_nu);
    }

    double Forward(double r) const {
        return std::pow(r * a_max + (1.0 - r) * a_min, inv_one_minus_nu) + m2;
    }

    double Inverse(double s) const {
        double a = std::pow(s - m2, one_minus_nu);
        return (a - a_min) / (a_max - a_min);
    }

    double Density(double s) const {
        double ds = s - m2;
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
struct TabulatedMapping {
    std::vector<double> s;        // node positions (ascending)
    std::vector<double> cdf;      // cumulative weight at each node (ascending)
    double s_min;
    double s_max;
    double c_lo;                  // interpolated cumulative at s_min
    double c_hi;                  // interpolated cumulative at s_max
    double norm;                  // c_hi - c_lo (cumulative mass in range)

    TabulatedMapping(std::vector<double> s_nodes,
                     std::vector<double> cdf_nodes,
                     double s_min_, double s_max_)
        : s(std::move(s_nodes)), cdf(std::move(cdf_nodes)),
          s_min(s_min_), s_max(s_max_)
    {
        if (s.size() < 2 || s.size() != cdf.size()) {
            throw std::runtime_error(
                "TabulatedMapping requires matching s/cdf arrays of length >= 2");
        }
        c_lo = InterpCdf(s_min);
        c_hi = InterpCdf(s_max);
        norm = c_hi - c_lo;
        if (!(norm > 0.0)) {
            throw std::runtime_error(
                "TabulatedMapping has non-positive cumulative mass over [s_min, s_max]");
        }
    }

    // Linear-interpolate the stored cumulative at an arbitrary s (clamped).
    double InterpCdf(double sv) const {
        if (sv <= s.front()) return cdf.front();
        if (sv >= s.back()) return cdf.back();
        std::size_t lo = 0, hi = s.size() - 1;
        while (hi - lo > 1) {
            std::size_t mid = (lo + hi) / 2;
            if (s[mid] <= sv) lo = mid; else hi = mid;
        }
        double ds = s[hi] - s[lo];
        double frac = ds > 0.0 ? (sv - s[lo]) / ds : 0.0;
        return cdf[lo] + frac * (cdf[hi] - cdf[lo]);
    }

    // Local density d(cdf)/ds at sv (bin slope), unnormalized.
    double SlopeAt(double sv) const {
        if (sv < s.front() || sv > s.back()) return 0.0;
        std::size_t lo = 0, hi = s.size() - 1;
        while (hi - lo > 1) {
            std::size_t mid = (lo + hi) / 2;
            if (s[mid] <= sv) lo = mid; else hi = mid;
        }
        double ds = s[hi] - s[lo];
        if (ds <= 0.0) return 0.0;
        return (cdf[hi] - cdf[lo]) / ds;
    }

    double Forward(double r) const {
        double target = c_lo + r * norm;             // cumulative we want
        // Find the bin [s[lo], s[hi]] whose cumulative brackets target.
        if (target <= cdf.front()) return s_min;
        if (target >= cdf.back()) return s_max;
        std::size_t lo = 0, hi = cdf.size() - 1;
        while (hi - lo > 1) {
            std::size_t mid = (lo + hi) / 2;
            if (cdf[mid] <= target) lo = mid; else hi = mid;
        }
        double dc = cdf[hi] - cdf[lo];
        double frac = dc > 0.0 ? (target - cdf[lo]) / dc : 0.0;
        double sv = s[lo] + frac * (s[hi] - s[lo]);
        // Guard against rounding outside the channel's allowed range.
        if (sv < s_min) sv = s_min;
        if (sv > s_max) sv = s_max;
        return sv;
    }

    double Inverse(double sv) const {
        if (sv <= s_min) return 0.0;
        if (sv >= s_max) return 1.0;
        return (InterpCdf(sv) - c_lo) / norm;
    }

    // Normalized proposal density in s over [s_min, s_max].
    double Density(double sv) const {
        if (sv < s_min || sv > s_max) return 0.0;
        return SlopeAt(sv) / norm;
    }
};

// Uniform mapping: s = s_min + r * (s_max - s_min)
struct UniformMapping {
    double s_min;
    double s_max;
    double range;

    UniformMapping(double s_min_, double s_max_)
        : s_min(s_min_), s_max(s_max_), range(s_max_ - s_min_) {}

    double Forward(double r) const { return s_min + r * range; }
    double Inverse(double s) const { return (s - s_min) / range; }
    double Density(double) const { return 1.0 / range; }
};

} // namespace injection
} // namespace siren

#endif // SIREN_InvariantMassMapping_H
