#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"

namespace siren { namespace injection { namespace detail {

// A conservative product of a polar interval union and one azimuth interval
// in the parent rest frame. The boost preserves azimuth. The forward polar
// map is monotone on each side of u=-v*/beta, so bisection finds all preimages,
// including the massive-daughter backward branch without inverse-root loss.
struct RestFrameEnvelope {
    using Vec = siren::math::Vector3D;
    Vec axis, first, second;
    std::vector<std::pair<double, double>> intervals;
    double half_phi = M_PI;
    double length = 2.0;
    bool active = false;

    double Area() const { return 2.0 * half_phi * length; }
    bool Contains(double u, double phi) const {
        if (!active) return true;
        if (std::abs(phi) > half_phi) return false;
        for (auto const & interval : intervals)
            if (u >= interval.first && u <= interval.second) return true;
        return false;
    }
};

inline RestFrameEnvelope BuildRestFrameEnvelope(
    double parent_energy, RestFrameEnvelope::Vec const & parent_p,
    double parent_mass, double rest_energy, double rest_momentum,
    RestFrameEnvelope::Vec const & position,
    siren::geometry::Geometry const & target)
{
    using Vec = RestFrameEnvelope::Vec;
    RestFrameEnvelope out;
    auto bounds = target.GetWorldBoundingBox();
    Vec displacement = (bounds.min_corner + bounds.max_corner)*0.5-position;
    double distance = displacement.magnitude();
    double radius = (bounds.max_corner-bounds.min_corner).magnitude()*0.5;
    double pp = parent_p.magnitude();
    out.axis = pp > 0 ? parent_p/pp : (distance > 0 ? displacement/distance : Vec(0,0,1));
    Vec to_center = distance > 0 ? displacement/distance : out.axis;
    double cosine = std::clamp(siren::math::scalar_product(out.axis, to_center), -1.0, 1.0);
    out.first = to_center-out.axis*cosine;
    if (out.first.magnitude() < 1e-12) {
        Vec reference = std::abs(out.axis.GetX()) < .9 ? Vec(1,0,0) : Vec(0,1,0);
        out.first = reference-out.axis*siren::math::scalar_product(reference,out.axis);
    }
    out.first.normalize();
    out.second = siren::math::vector_product(out.axis,out.first);
    out.second.normalize();
    out.intervals = {{-1.0,1.0}};
    if (!(distance > radius) || !(radius > 0) || !(rest_momentum > 0)) return out;
    double alpha = std::acos(cosine), delta = std::asin(std::min(1.0,radius/distance));
    double cmin = std::cos(std::min(M_PI,alpha+delta));
    double cmax = std::cos(std::max(0.0,alpha-delta));
    double beta = pp/parent_energy, gamma = parent_energy/parent_mass;
    auto lab_cosine = [&](double u) {
        double z = gamma*(rest_momentum*u+beta*rest_energy);
        double transverse = rest_momentum*std::sqrt(std::max(0.0,1-u*u));
        double norm = std::hypot(z,transverse);
        return norm > 0 ? z/norm : 0.0;
    };
    std::vector<double> breaks{-1.0};
    if (beta*rest_energy > rest_momentum)
        breaks.push_back(-rest_momentum/(beta*rest_energy));
    breaks.push_back(1.0);
    out.intervals.clear();
    for (std::size_t i=1; i<breaks.size(); ++i) {
        double a=breaks[i-1], b=breaks[i], fa=lab_cosine(a), fb=lab_cosine(b);
        double low=std::max(cmin,std::min(fa,fb));
        double high=std::min(cmax,std::max(fa,fb));
        if (!(high > low)) continue;
        auto invert = [&](double value) {
            if (value == fa) return a;
            if (value == fb) return b;
            double left=a, right=b;
            for (int step=0; step<60; ++step) {
                double middle=left+(right-left)*.5;
                if ((lab_cosine(middle)<value) == (fb>fa)) left=middle;
                else right=middle;
            }
            return left+(right-left)*.5;
        };
        double u0=invert(low), u1=invert(high);
        if (u0>u1) std::swap(u0,u1);
        // Enlarge by roundoff, never remove a finite angular strip.
        double pad=64*std::numeric_limits<double>::epsilon();
        u0=std::max(-1.0,u0-pad); u1=std::min(1.0,u1+pad);
        if (u1>u0) out.intervals.emplace_back(u0,u1);
    }
    std::sort(out.intervals.begin(),out.intervals.end());
    if (out.intervals.size()==2 && out.intervals[0].second>=out.intervals[1].first) {
        out.intervals[0].second=out.intervals[1].second;
        out.intervals.resize(1);
    }
    out.length=0;
    for (auto const & interval : out.intervals) out.length+=interval.second-interval.first;
    // A cone containing either boost-axis pole spans every azimuth.
    if (alpha>delta && M_PI-alpha>delta) {
        out.half_phi=std::asin(std::clamp(std::sin(delta)/std::sin(alpha),0.0,1.0));
        out.half_phi=std::min(M_PI,out.half_phi+64*std::numeric_limits<double>::epsilon());
    }
    // A preimage narrower than representable boost-angle resolution cannot
    // support stable inversion/density evaluation. Use full physical support
    // rather than pruning it or emitting a nearly singular numerical box.
    if (out.length <= 1024*std::numeric_limits<double>::epsilon()
        || !(out.Area()>0) || !std::isfinite(1/out.Area())) {
        out.intervals={{-1.0,1.0}}; out.length=2; out.half_phi=M_PI;
        return out;
    }
    out.active=out.Area()<4*M_PI;
    return out;
}

}}} // namespace siren::injection::detail
