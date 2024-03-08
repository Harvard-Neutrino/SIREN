#pragma once
#ifndef LI_Coordinates_H
#define LI_Coordinates_H

#include <utility>
#include "SIREN/math/Vector3D.h"           // for Vector3D
#include <NamedType/named_type.hpp>

namespace fluent {

template<typename O>
struct A {

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

};

} // namespace fluent


namespace SI {
namespace detector {

using GeometryPosition = fluent::NamedType<SI::math::Vector3D, struct GeometryPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::A<double>::MultiplicableBy, fluent::A<double>::DivisibleBy>;
using GeometryDirection = fluent::NamedType<SI::math::Vector3D, struct GeometryDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::A<double>::MultiplicableBy, fluent::A<double>::DivisibleBy>;
using DetectorPosition = fluent::NamedType<SI::math::Vector3D, struct DetectorPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::A<double>::MultiplicableBy, fluent::A<double>::DivisibleBy>;
using DetectorDirection = fluent::NamedType<SI::math::Vector3D, struct DetectorDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable, fluent::A<double>::MultiplicableBy, fluent::A<double>::DivisibleBy>;

} // namespace detector
} // namespace SI

#endif // LI_Coordinates_H
