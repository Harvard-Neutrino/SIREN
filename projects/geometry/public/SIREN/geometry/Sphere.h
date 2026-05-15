#pragma once
#ifndef SIREN_Sphere_H
#define SIREN_Sphere_H

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

class Sphere : public Geometry {
public:
    Sphere();
    Sphere(double radius, double inner_radius);
    Sphere(Placement const &);
    Sphere(Placement const &, double radius, double inner_radius);
    Sphere(double radius, double inner_radius,
           double start_phi, double delta_phi,
           double start_theta, double delta_theta);
    Sphere(Placement const &, double radius, double inner_radius,
           double start_phi, double delta_phi,
           double start_theta, double delta_theta);
    Sphere(const Sphere&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("OuterRadius", radius_));
            archive(::cereal::make_nvp("InnerRadius", inner_radius_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("DeltaPhi", delta_phi_));
            archive(::cereal::make_nvp("StartTheta", start_theta_));
            archive(::cereal::make_nvp("DeltaTheta", delta_theta_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Sphere only supports version <= 0!");
        }
    }

    /* Geometry* clone() const override{ return new Sphere(*this); }; */
    std::shared_ptr<Geometry> create() const override{ return std::shared_ptr<Geometry>( new Sphere(*this) ); }
    void swap(Geometry&) override;

    virtual ~Sphere() {}

    // Operators
    Sphere& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }
    double GetStartPhi() const { return start_phi_; }
    double GetDeltaPhi() const { return delta_phi_; }
    double GetStartTheta() const { return start_theta_; }
    double GetDeltaTheta() const { return delta_theta_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double radius_;       //!< the radius of the sphere/ cylinder
    double inner_radius_; //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
    double start_phi_;    //!< starting azimuthal angle (radians, default 0)
    double delta_phi_;    //!< azimuthal extent (radians, default 2*pi)
    double start_theta_;  //!< starting polar angle from +z (radians, default 0)
    double delta_theta_;  //!< polar extent (radians, default pi)
    bool has_phi_cut_;    //!< true if delta_phi < 2*pi (precomputed)
    bool has_theta_cut_;  //!< true if start_theta > 0 or delta_theta < pi (precomputed)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Sphere, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Sphere);
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Sphere);

CEREAL_FORCE_DYNAMIC_INIT(siren_Sphere);

#endif // SIREN_Sphere_H
