#pragma once
#ifndef LI_Coordinates_H
#define LI_Coordinates_H

#include "LeptonInjector/math/Vector3D.h"           // for Vector3D
#include <NamedType/named_type.hpp>

namespace fluent {

template<typename O>
struct AltArithmetic {

template <typename T>
struct MultiplicableBy : crtp<T, MultiplicableBy>
{
    FLUENT_NODISCARD constexpr T operator*(O const& other) const
    {
        return T(this->underlying().get() * other);
    }
    FLUENT_CONSTEXPR17 T& operator*=(O const& other)
    {
        this->underlying().get() *= other;
        return this->underlying();
    }
};

template <typename T>
struct DivisibleBy : crtp<T, DivisibleBy>
{
    FLUENT_NODISCARD constexpr T operator/(O const& other) const
    {
        return T(this->underlying().get() / other);
    }
    FLUENT_CONSTEXPR17 T& operator/=(O const& other)
    {
        this->underlying().get() /= other;
        return this->underlying();
    }
};

template <typename T>
struct MultiplicableTo : crtp<T, MultiplicableTo>
{
    FLUENT_NODISCARD constexpr O operator*(T const& other) const
    {
        return O(this->underlying().get() * other.get());
    }
};

};

} // namespace fluent


namespace LI {
namespace detector {

using GeometryPosition = fluent::NamedType<LI::math::Vector3D, struct GeometryPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::AltArithmetic<double>::MultiplicableBy, fluent::AltArithmetic<double>::DivisibleBy, fluent::AltArithmetic<double>::MultiplicableTo>;
using GeometryDirection = fluent::NamedType<LI::math::Vector3D, struct GeometryDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::AltArithmetic<double>::MultiplicableBy, fluent::AltArithmetic<double>::DivisibleBy, fluent::AltArithmetic<double>::MultiplicableTo>;
using DetectorPosition = fluent::NamedType<LI::math::Vector3D, struct DetectorPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::AltArithmetic<double>::MultiplicableBy, fluent::AltArithmetic<double>::DivisibleBy, fluent::AltArithmetic<double>::MultiplicableTo>;
using DetectorDirection = fluent::NamedType<LI::math::Vector3D, struct DetectorDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::AltArithmetic<double>::MultiplicableBy, fluent::AltArithmetic<double>::DivisibleBy, fluent::AltArithmetic<double>::MultiplicableTo>;

} // namespace detector
} // namespace LI

#endif // LI_Coordinates_H
