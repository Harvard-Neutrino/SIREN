#pragma once
#ifndef LI_Cylinder_H
#define LI_Cylinder_H

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

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/geometry/Geometry.h"

namespace LI {
namespace geometry {

class Cylinder : public Geometry {
public:
    Cylinder();
    Cylinder(double radius, double inner_radius, double z);
    Cylinder(Placement const &);
    Cylinder(Placement const &, double radius, double inner_radius, double z);
    Cylinder(const Cylinder&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("OuterRadius", radius_));
            archive(::cereal::make_nvp("InnerRadius", inner_radius_));
            archive(::cereal::make_nvp("Height", z_));
            archive(cereal::virtual_base_class<Geometry>(this));
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
    std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;

    // Getter & Setter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }
    double GetZ() const { return z_; }

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
};


} // namespace geometry
} // namespace LI

CEREAL_CLASS_VERSION(LI::geometry::Cylinder, 0);
CEREAL_REGISTER_TYPE(LI::geometry::Cylinder)
CEREAL_REGISTER_POLYMORPHIC_RELATION(LI::geometry::Geometry, LI::geometry::Cylinder);

#endif // LI_Cylinder_H
