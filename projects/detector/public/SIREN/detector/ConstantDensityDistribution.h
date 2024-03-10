#pragma once
#ifndef SIREN_ConstantDensityDistribution_H
#define SIREN_ConstantDensityDistribution_H
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
#include "SIREN/detector/ConstantDistribution1D.h"

namespace siren {
namespace detector {

template<typename AxisT>
class DensityDistribution1D<AxisT, ConstantDistribution1D, typename std::enable_if<std::is_base_of<Axis1D, AxisT>::value>::type> : public DensityDistribution {
    using DistributionT = ConstantDistribution1D;
    using T = DensityDistribution1D<AxisT, ConstantDistribution1D>;
   private:
    AxisT axis;
    ConstantDistribution1D dist;
   public:
    DensityDistribution1D() : axis(), dist() {};
    DensityDistribution1D(const AxisT& axis, const ConstantDistribution1D& dist)
        : axis(axis), dist(dist) {};
    DensityDistribution1D(const AxisT& axis, double val)
        : axis(axis), dist(val) {};
    DensityDistribution1D(double val)
        : axis(), dist(val) {};
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
    std::shared_ptr<DensityDistribution> create() const override {
        return std::shared_ptr<DensityDistribution>(new T(*this));
    };

    double Derivative(const math::Vector3D& xi,
                      const math::Vector3D& direction) const override {
        return 0;
    };

    double AntiDerivative(const math::Vector3D& xi,
                          const math::Vector3D& direction) const override {
        //return (xi*axis.GetAxis())/(direction*axis.GetAxis()) * dist.Evaluate;
        return (xi*direction)*dist.Evaluate(0);
    };

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& direction,
                    double distance) const override {
        (void)direction;
        return distance*dist.Evaluate(0);
    }

    double Integral(const math::Vector3D& xi,
                    const math::Vector3D& xj) const override {
        double distance = (xj-xi).magnitude();
        return distance*dist.Evaluate(0);
    };

    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double integral,
                           double max_distance) const override {
        (void)direction;
        double distance = integral / dist.Evaluate(0);
        if(distance > max_distance) {
            distance = -1;
        }
        return distance;
    };
    
    double InverseIntegral(const math::Vector3D& xi,
                           const math::Vector3D& direction,
                           double constant,
                           double integral,
                           double max_distance) const override {
        (void)direction;
        double distance = integral / (dist.Evaluate(0) + constant);
        if(distance > max_distance) {
            distance = -1;
        }
        return distance;
    };

    double Evaluate(const math::Vector3D& xi) const override {
        (void)xi;
        return dist.Evaluate(0);
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

// Define the specialization that we plan to use
typedef DensityDistribution1D<CartesianAxis1D, ConstantDistribution1D> ConstantDensityDistribution;

// Technically we could define this with the radial axis as well, but there is not point since it will have the same functionality
// typedef DensityDistribution1D<siren::detector::RadialAxis1D ConstantDensityDistribution;

} // namespace detector
} // namespace siren

CEREAL_CLASS_VERSION(siren::detector::ConstantDensityDistribution, 0);
CEREAL_REGISTER_TYPE(siren::detector::ConstantDensityDistribution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::detector::DensityDistribution, siren::detector::ConstantDensityDistribution);

#endif // SIREN_ConstantDensityDistribution.h
