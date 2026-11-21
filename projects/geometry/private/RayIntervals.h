// Interval-set algebra along a ray parameter, shared by the solid-path hit
// generation of shapes assembled from quadric sign conditions (Polycone,
// Sphere). The in-solid ray-parameter set of such a shape is built here
// from closed quadratic and half-line conditions with exact set operations;
// emitting each connected piece of the final set as an (enter, exit) pair
// makes the hit list parity-consistent by construction, with no per-surface
// tolerance windows that can disagree at surface seams (the failure mode
// fixed for Cylinder by the same interval method, done inline there).

#pragma once
#ifndef SIREN_RayIntervals_H
#define SIREN_RayIntervals_H

#include <cmath>
#include <limits>
#include <cassert>
#include <algorithm>

namespace siren {
namespace geometry {
namespace ray_intervals {

// A set of disjoint closed intervals [lo, hi] with lo < hi, sorted by lo.
// Endpoints may be +-infinity. Every operation below preserves the
// invariant. CAP is far above what the algebra in this codebase can
// produce (each constraint contributes at most two pieces and only a
// handful of constraints are combined per solid); if it were ever hit,
// whole pieces would be dropped, which loses a sliver of volume but can
// never break enter/exit pairing.
struct IntervalSet {
    static constexpr int CAP = 8;
    double lo[CAP];
    double hi[CAP];
    int n = 0;

    // Append a piece; the caller is responsible for keeping the sorted,
    // disjoint invariant. Empty and inverted pieces are dropped (the
    // comparison also rejects NaN endpoints).
    void Add(double a, double b) {
        if(!(a < b)) return;
        assert(n < CAP);
        if(n >= CAP) return;
        lo[n] = a;
        hi[n] = b;
        ++n;
    }

    void Shift(double dt) {
        for(int i = 0; i < n; ++i) {
            lo[i] += dt;
            hi[i] += dt;
        }
    }

    static IntervalSet Full() {
        IntervalSet s;
        double const inf = std::numeric_limits<double>::infinity();
        s.Add(-inf, inf);
        return s;
    }
};

// Intersection of two sets: two-pointer sweep over the sorted piece lists.
inline IntervalSet Intersect(IntervalSet const & A, IntervalSet const & B) {
    IntervalSet out;
    int i = 0, j = 0;
    while(i < A.n && j < B.n) {
        out.Add(std::max(A.lo[i], B.lo[j]), std::min(A.hi[i], B.hi[j]));
        if(A.hi[i] < B.hi[j]) ++i;
        else ++j;
    }
    return out;
}

// Union of two sets: merge-sweep; pieces that overlap or touch exactly
// (shared endpoint, bit-identical) coalesce into one.
inline IntervalSet Union(IntervalSet const & A, IntervalSet const & B) {
    IntervalSet out;
    int i = 0, j = 0;
    double cur_lo = 0, cur_hi = 0;
    bool open = false;
    while(i < A.n || j < B.n) {
        double a, b;
        if(j >= B.n || (i < A.n && A.lo[i] <= B.lo[j])) {
            a = A.lo[i]; b = A.hi[i]; ++i;
        } else {
            a = B.lo[j]; b = B.hi[j]; ++j;
        }
        if(!open) {
            cur_lo = a; cur_hi = b; open = true;
        } else if(a <= cur_hi) {
            cur_hi = std::max(cur_hi, b);
        } else {
            out.Add(cur_lo, cur_hi);
            cur_lo = a; cur_hi = b;
        }
    }
    if(open) out.Add(cur_lo, cur_hi);
    return out;
}

// A minus B. Remnant pieces keep the boundary parameters of B (closed
// remnants), so a hit landing on a subtracted surface stays on it.
inline IntervalSet Subtract(IntervalSet const & A, IntervalSet const & B) {
    IntervalSet out;
    for(int i = 0; i < A.n; ++i) {
        double cursor = A.lo[i];
        double end = A.hi[i];
        for(int j = 0; j < B.n && cursor < end; ++j) {
            if(B.hi[j] <= cursor) continue;
            if(B.lo[j] >= end) break;
            out.Add(cursor, std::min(end, B.lo[j]));
            cursor = std::max(cursor, B.hi[j]);
        }
        out.Add(cursor, end);
    }
    return out;
}

// The closed set {t : A*t^2 + 2*B_half*t + C <= 0}. Zero-measure solutions
// (tangencies) come out empty, matching the tangent-line-is-a-miss policy
// of the interval method. Roots use the numerically stable split (one root
// from the large-magnitude numerator, the other from C over it), so a tiny
// leading coefficient degrades gracefully toward the linear case instead
// of amplifying cancellation; an overflowed root becomes an infinite
// endpoint, which is the correct limit.
inline IntervalSet QuadraticLEQ(double A, double B_half, double C) {
    IntervalSet out;
    double const inf = std::numeric_limits<double>::infinity();
    if(A == 0) {
        if(B_half > 0) out.Add(-inf, -C / (2.0 * B_half));
        else if(B_half < 0) out.Add(-C / (2.0 * B_half), inf);
        else if(C <= 0) out.Add(-inf, inf);
        return out;
    }
    double det = B_half * B_half - A * C;
    if(det <= 0) {
        // No sign change: the parabola stays on one side (touching zero at
        // most at a point). Below the axis everywhere iff it opens down.
        if(A < 0) out.Add(-inf, inf);
        return out;
    }
    double sq = std::sqrt(det);
    double q = -(B_half + std::copysign(sq, B_half));
    double t1 = q / A;
    double t2 = C / q;
    if(t1 > t2) std::swap(t1, t2);
    if(A > 0) {
        out.Add(t1, t2);
    } else {
        out.Add(-inf, t1);
        out.Add(t2, inf);
    }
    return out;
}

inline IntervalSet QuadraticGEQ(double A, double B_half, double C) {
    return QuadraticLEQ(-A, -B_half, -C);
}

// The closed set {t : pz + t*dz >= 0}. The complement's boundary parameter
// is bit-identical (negate both arguments), so unions of the two half-lines
// coalesce exactly at the plane crossing.
inline IntervalSet HalfLineGEQ(double pz, double dz) {
    IntervalSet out;
    double const inf = std::numeric_limits<double>::infinity();
    if(dz > 0) out.Add(-pz / dz, inf);
    else if(dz < 0) out.Add(-inf, -pz / dz);
    else if(pz >= 0) out.Add(-inf, inf);
    return out;
}

} // namespace ray_intervals
} // namespace geometry
} // namespace siren

#endif // SIREN_RayIntervals_H
