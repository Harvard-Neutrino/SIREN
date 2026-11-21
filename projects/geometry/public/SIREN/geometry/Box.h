#pragma once
#ifndef SIREN_Box_H
#define SIREN_Box_H

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

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Geometry.h"

namespace siren {
namespace geometry {

class Box : public Geometry {
public:
    Box();
    Box(double full_width_x, double full_width_y, double full_width_z);
    Box(Placement const &);
    Box(Placement const &, double x, double y, double z);
    Box(const Box&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("XFullWidth", full_width_x_));
            archive(::cereal::make_nvp("YFullWidth", full_width_y_));
            archive(::cereal::make_nvp("ZFullWidth", full_width_z_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Box only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Box(*this) ); };
    void swap(Geometry&) override;

    virtual ~Box() {}

    // Operators
    Box& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    double GetX() const { return full_width_x_; }
    double GetY() const { return full_width_y_; }
    double GetZ() const { return full_width_z_; }

    void SetX(double full_width_x) { full_width_x_ = full_width_x; RecomputeWorldAABB(); };
    void SetY(double full_width_y) { full_width_y_ = full_width_y; RecomputeWorldAABB(); };
    void SetZ(double full_width_z) { full_width_z_ = full_width_z; RecomputeWorldAABB(); };
protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double full_width_x_; //!< full width of box in x-direction
    double full_width_y_; //!< full width of box in y-direction
    double full_width_z_; //!< full width of box in z-direction
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Box, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Box)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Box);

#endif // SIREN_Box_H
