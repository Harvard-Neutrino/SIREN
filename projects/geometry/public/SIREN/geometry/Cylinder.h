#pragma once
#ifndef SIREN_Cylinder_H
#define SIREN_Cylinder_H

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

class Cylinder : public Geometry {
public:
    Cylinder();
    Cylinder(double radius, double inner_radius, double z);
    Cylinder(Placement const &);
    Cylinder(Placement const &, double radius, double inner_radius, double z);
    Cylinder(double radius, double inner_radius, double z, double start_phi, double delta_phi);
    Cylinder(Placement const &, double radius, double inner_radius, double z, double start_phi, double delta_phi);
    Cylinder(const Cylinder&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("OuterRadius", radius_));
            archive(::cereal::make_nvp("InnerRadius", inner_radius_));
            archive(::cereal::make_nvp("Height", z_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("DeltaPhi", delta_phi_));
            archive(cereal::virtual_base_class<Geometry>(this));
            has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
        } else {
            throw std::runtime_error("Cylinder only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Cylinder(*this) ); };
    void swap(Geometry&) override;

    virtual ~Cylinder() {}

    // Operators
    Cylinder& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }
    double GetZ() const { return z_; }
    double GetStartPhi() const { return start_phi_; }
    double GetDeltaPhi() const { return delta_phi_; }

    void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
    void SetRadius(double radius) { radius_ = radius; };
    void SetZ(double z) { z_ = z; };

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double radius_;       //!< the radius of the sphere/ cylinder
    double inner_radius_; //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
    double z_;            //!< height of box/cylinder
    double start_phi_;    //!< starting azimuthal angle (radians, default 0)
    double delta_phi_;    //!< azimuthal extent (radians, default 2*pi)
    bool has_phi_cut_;    //!< true if delta_phi < 2*pi (precomputed)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Cylinder, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Cylinder)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Cylinder);

#endif // SIREN_Cylinder_H
