#pragma once
#ifndef LI_Sphere_H
#define LI_Sphere_H

#include <map>
#include <memory>
#include <vector>
#include <math.h>
#include <float.h>
#include <iostream>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/geometry/Geometry.h"

namespace LI {
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
    std::shared_ptr<const Geometry> create() const override{ return std::shared_ptr<const Geometry>( new Sphere(*this) ); }
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
} // namespace LI

CEREAL_CLASS_VERSION(LI::geometry::Sphere, 0);
CEREAL_REGISTER_TYPE(LI::geometry::Sphere)
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::geometry::Geometry, LI::geometry::Sphere);

#endif // LI_Sphere_H
