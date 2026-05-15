#pragma once
#ifndef SIREN_Torus_H
#define SIREN_Torus_H

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

// Torus of revolution about the z-axis.
//
// The torus is defined by:
//   rtor  - major radius (distance from z-axis to center of tube cross-section)
//   rmax  - outer radius of the tube cross-section
//   rmin  - inner radius of the tube cross-section (0 for solid torus)
//
// The surface equation is:
//   (sqrt(x^2 + y^2) - rtor)^2 + z^2 = r^2
//
// Full rotation only (no startphi / deltaphi support).
class Torus : public Geometry {
public:
    Torus();
    Torus(double rtor, double rmax, double rmin);
    Torus(Placement const &);
    Torus(Placement const &, double rtor, double rmax, double rmin);
    Torus(const Torus&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("MajorRadius", rtor_));
            archive(::cereal::make_nvp("OuterRadius", rmax_));
            archive(::cereal::make_nvp("InnerRadius", rmin_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Torus only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>(new Torus(*this)); };
    void swap(Geometry&) override;

    virtual ~Torus() {}

    // Operators
    Torus& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getters
    double GetMajorRadius() const { return rtor_; }
    double GetMinorRadius() const { return rmax_; }
    double GetInnerRadius() const { return rmin_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double rtor_; //!< major radius (z-axis to tube center)
    double rmax_; //!< outer tube radius
    double rmin_; //!< inner tube radius (0 for solid)
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Torus, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Torus)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Torus);

#endif // SIREN_Torus_H
