#pragma once
#ifndef SIREN_Cone_H
#define SIREN_Cone_H

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

class Cone : public Geometry {
public:
    Cone();
    Cone(double rmin1, double rmax1, double rmin2, double rmax2, double z);
    Cone(Placement const &);
    Cone(Placement const &, double rmin1, double rmax1, double rmin2, double rmax2, double z);
    Cone(const Cone&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InnerRadius1", rmin1_));
            archive(::cereal::make_nvp("OuterRadius1", rmax1_));
            archive(::cereal::make_nvp("InnerRadius2", rmin2_));
            archive(::cereal::make_nvp("OuterRadius2", rmax2_));
            archive(::cereal::make_nvp("Height", z_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Cone only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Cone(*this) ); };
    void swap(Geometry&) override;

    virtual ~Cone() {}

    // Operators
    Cone& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    double GetRmin1() const { return rmin1_; }
    double GetRmax1() const { return rmax1_; }
    double GetRmin2() const { return rmin2_; }
    double GetRmax2() const { return rmax2_; }
    double GetZ() const { return z_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double rmin1_; //!< inner radius at z = -z/2
    double rmax1_; //!< outer radius at z = -z/2
    double rmin2_; //!< inner radius at z = +z/2
    double rmax2_; //!< outer radius at z = +z/2
    double z_;     //!< full height along z-axis
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Cone, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Cone)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Cone);

#endif // SIREN_Cone_H
