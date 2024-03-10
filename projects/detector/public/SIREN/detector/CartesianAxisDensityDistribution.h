#pragma once
#ifndef SIREN_CartesianAxisDensityDistribution_H
#define SIREN_CartesianAxisDensityDistribution_H
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

#include "SIREN/detector/Axis1D.h"
#include "SIREN/detector/CartesianAxis1D.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/DensityDistribution1D.h"

namespace siren {
namespace detector {

template<typename DistributionT>
class DensityDistribution1D<CartesianAxis1D, DistributionT, typename std::enable_if<std::is_base_of<Distribution1D, DistributionT>::value && !std::is_same<ConstantDistribution1D,DistributionT>::value>::type> : public DensityDistribution {
    using AxisT = CartesianAxis1D;
    using T = DensityDistribution1D<CartesianAxis1D, DistributionT>;
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
        const T* d_1d = dynamic_cast<const T*>(&d);
        if(!d_1d)
            return false;
        if(axis != d_1d->axis or dist != d_1d->dist)
            return false;
        return true;
    };

    DensityDistribution* clone() const override { return new T(*this); };
    std::shared_ptr<DensityDistribution> create() override {
        return std::shared_ptr<const DensityDistribution>(new T(*this));
    };

    double Derivative(const math::Vector3D& xi,
                      const math::Vector3D& direction) const override {
        return dist.Derivative(axis.GetX(xi))*axis.GetdX(xi, direction);
    };

    double AntiDerivative(const math::Vector3D& xi,
                          const math::Vector3D& direction) const override {
        double a = 0;
        double b = axis.GetX(xi);
        double dxdt = axis.GetdX(xi, direction);
        return (dist.AntiDerivative(b) - dist.AntiDerivative(a)) / dxdt;
    };

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& direction,
                    double distance) const override {
        double a = axis.GetX(xi);
        math::Vector3D xj = xi + direction*distance;
        double b = axis.GetX(xj);
        double dxdt = axis.GetdX(xi, direction);
        return (dist.AntiDerivative(b) - dist.AntiDerivative(a)) / dxdt;
    }

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& xj) const override {
        double a = axis.GetX(xi);
        double b = axis.GetX(xj);
        math::Vector3D direction = xj-xi;
        direction.normalize();
        double dxdt = axis.GetdX(xi, direction);
        return (dist.AntiDerivative(b) - dist.AntiDerivative(a)) / dxdt;
    };

    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double integral,
                           double max_distance) const override {
        double a = axis.GetX(xi);
        double b = axis.GetX(xi + direction*max_distance);
        double dxdt = axis.GetdX(xi, direction);

        double dist_integral = integral * dxdt;

        double Ia = dist.AntiDerivative(a);
        std::function<double(double)> F = [&](double x)->double {
            return (dist.AntiDerivative(x) - Ia) - dist_integral;
        };

        std::function<double(double)> dF = [&](double x)->double {
            return dist.Evaluate(x);
        };

        try {
            double b_res = siren::math::NewtonRaphson(F, dF, a, b, (a+b)/2.0);
            return (b_res - a)/dxdt;
        } catch(siren::math::MathException& e) {
            return -1;
        }
    };

    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double constant,
                           double integral,
                           double max_distance) const override {
        double a = axis.GetX(xi);
        double b = axis.GetX(xi + direction*max_distance);
        double dxdt = axis.GetdX(xi, direction);

        double dist_integral = integral * dxdt;

        double Ia = dist.AntiDerivative(a);
        std::function<double(double)> F = [&](double x)->double {
            return (dist.AntiDerivative(x) - Ia) + constant*x - dist_integral;
        };

        std::function<double(double)> dF = [&](double x)->double {
            return dist.Evaluate(x) + constant;
        };

        try {
            double b_res = siren::math::NewtonRaphson(F, dF, a, b, (a+b)/2.0);
            return (b_res - a)/dxdt;
        } catch(siren::math::MathException& e) {
            return -1;
        }
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

#endif // SIREN_CartesianAxisDensityDistribution.h
