#pragma once
#ifndef SIREN_BooleanGeometry_H
#define SIREN_BooleanGeometry_H

#include <memory>
#include <vector>
#include <cstdint>
#include <utility>
#include <iostream>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/memory.hpp>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Geometry.h"

namespace siren {
namespace geometry {

enum class BooleanOperation { UNION, SUBTRACTION, INTERSECTION };

class BooleanGeometry : public Geometry {
public:
    BooleanGeometry();
    BooleanGeometry(BooleanOperation op, std::shared_ptr<const Geometry> left, std::shared_ptr<const Geometry> right);
    BooleanGeometry(Placement const &, BooleanOperation op, std::shared_ptr<const Geometry> left, std::shared_ptr<const Geometry> right);
    BooleanGeometry(const BooleanGeometry&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Operation", op_));
            archive(::cereal::make_nvp("Left", left_));
            archive(::cereal::make_nvp("Right", right_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("BooleanGeometry only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>(new BooleanGeometry(*this)); }
    void swap(Geometry&) override;

    virtual ~BooleanGeometry() {}

    // Operators
    BooleanGeometry& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    BooleanOperation GetOperation() const { return op_; }
    std::shared_ptr<const Geometry> GetLeft() const { return left_; }
    std::shared_ptr<const Geometry> GetRight() const { return right_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    BooleanOperation op_;                  //!< boolean operation type
    std::shared_ptr<const Geometry> left_;  //!< left (first) operand
    std::shared_ptr<const Geometry> right_; //!< right (second) operand
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::BooleanGeometry, 0);
CEREAL_REGISTER_TYPE(siren::geometry::BooleanGeometry);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::BooleanGeometry);

CEREAL_FORCE_DYNAMIC_INIT(siren_BooleanGeometry);

#endif // SIREN_BooleanGeometry_H
