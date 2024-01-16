#pragma once
#ifndef LI_Coordinates_H
#define LI_Coordinates_H

#include "LeptonInjector/math/Vector3D.h"           // for Vector3D
#include <NamedType/named_type.hpp>


namespace LI {
namespace detector {

using GeometryPosition = fluent::NamedType<LI::math::Vector3D, struct GeometryPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable>;
using GeometryDirection = fluent::NamedType<LI::math::Vector3D, struct GeometryDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable>;
using DetectorPosition = fluent::NamedType<LI::math::Vector3D, struct DetectorPositionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable>;
using DetectorDirection = fluent::NamedType<LI::math::Vector3D, struct DetectorDirectionTag, fluent::Callable, fluent::Comparable, fluent::BinaryAddable, fluent::Subtractable>;

} // namespace detector
} // namespace LI

#endif // LI_Coordinates_H
