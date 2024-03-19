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
    Sphere(const Sphere&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("OuterRadius", radius_));
            archive(::cereal::make_nvp("InnerRadius", inner_radius_));
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
    std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;

    // Getter & Setter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }

    void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
    void SetRadius(double radius) { radius_ = radius; };

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double radius_;       //!< the radius of the sphere/ cylinder
    double inner_radius_; //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Sphere, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Sphere)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Sphere);

#endif // SIREN_Sphere_H
