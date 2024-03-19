#pragma once
#ifndef SIREN_DensityDistribution1D_H
#define SIREN_DensityDistribution1D_H
#include <memory>
#include <string>
#include <exception>
#include <functional>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Polynomial.h"

#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/Distribution1D.h"
#include "SIREN/detector/DensityDistribution.h"

#include "SIREN/utilities/Integration.h"

namespace siren {
namespace detector {

template <typename AxisT, typename DistributionT, class E = typename std::enable_if<std::is_base_of<Axis1D, AxisT>::value && std::is_base_of<Distribution1D, DistributionT>::value>::type>
class DensityDistribution1D
    : public DensityDistribution {
    using T = DensityDistribution1D<AxisT,DistributionT>;
private:
    AxisT axis;
    DistributionT dist;
public:
    DensityDistribution1D() : axis(), dist() {};
    DensityDistribution1D(const AxisT& axis, const DistributionT& dist)
        : axis(axis), dist(dist) {};
    DensityDistribution1D(const DensityDistribution1D& other)
        : axis(other.axis), dist(other.dist) {};

    bool compare(const DensityDistribution& d) const override {
        const DensityDistribution1D* d_1d = dynamic_cast<const DensityDistribution1D*>(&d);
        if(!d_1d)
            return false;
        if(axis != d_1d->axis or dist != d_1d->dist)
            return false;
        return true;
    };

    DensityDistribution* clone() const override { return new T(*this); };
    std::shared_ptr<DensityDistribution> create() const override {
        return std::shared_ptr<DensityDistribution>(new T(*this));
    };

    double Derivative(const math::Vector3D& xi,
                      const math::Vector3D& direction) const override {
        return dist.Derivative(axis.GetX(xi))*axis.GetdX(xi, direction);
    };

    double AntiDerivative(const math::Vector3D& xi,
                          const math::Vector3D& direction) const override {
        math::Vector3D cap = xi - (xi*direction)*direction;
        return Integral(cap, direction, xi*direction);
    };

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& direction,
                    double distance) const override {
        std::function<double(double)> f = [&](double x)->double {
            return Evaluate(xi+x*direction);
        };
        return siren::utilities::rombergIntegrate(f, 0, distance, 1e-6);
    }

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& xj) const override {
        math::Vector3D direction = xj-xi;
        double distance = direction.magnitude();
        direction.normalize();
        return Integral(xi, direction, distance);
    };

    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double integral,
                           double max_distance) const override {
        std::function<double(double)> F = [&](double x)->double {
            return Integral(xi, direction, x) - integral;
        };

        std::function<double(double)> dF = [&](double x)->double {
            return Evaluate(xi+direction*x);
        };

        double res;
        try {
            double init = max_distance / 2.0;
            if(std::isinf(init)) {
                init = dF(0.0);
            }
            res = siren::math::NewtonRaphson(F, dF, 0, max_distance, init);
        } catch(siren::math::MathException& e) {
            res = -1;
        }
        return res;
    };
    
    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double constant,
                           double integral,
                           double max_distance) const override {
        std::function<double(double)> F = [&](double x)->double {
            return Integral(xi, direction, x) + constant*x - integral;
        };

        std::function<double(double)> dF = [&](double x)->double {
            return Evaluate(xi+direction*x) + constant;
        };

        double res;
        try {
            double init = max_distance / 2.0;
            if(std::isinf(init)) {
                init = dF(0.0);
            }
            res = siren::math::NewtonRaphson(F, dF, 0, max_distance, init);
        } catch(siren::math::MathException& e) {
            res = -1;
        }
        return res;
    };

    double Evaluate(const math::Vector3D& xi) const override {
        return dist.Evaluate(axis.GetX(xi));
    };

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Axis", axis));
            archive(::cereal::make_nvp("Distribution", dist));
            archive(cereal::virtual_base_class<DensityDistribution>(this));
        } else {
            throw std::runtime_error("DensityDistribution1D only supports version <= 0");
        }
    };
};

} // namespace detector
} // namespace siren

#endif // SIREN_DensityDistribution1D.h
