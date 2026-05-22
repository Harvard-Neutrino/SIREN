#pragma once
#ifndef SIREN_EllipticalTube_H
#define SIREN_EllipticalTube_H

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

class EllipticalTube : public Geometry {
public:
    EllipticalTube();
    EllipticalTube(double dx, double dy, double dz);
    EllipticalTube(Placement const &);
    EllipticalTube(Placement const &, double dx, double dy, double dz);
    EllipticalTube(const EllipticalTube&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("SemiAxisX", dx_));
            archive(::cereal::make_nvp("SemiAxisY", dy_));
            archive(::cereal::make_nvp("HalfZ", dz_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("EllipticalTube only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>(new EllipticalTube(*this)); }
    void swap(Geometry&) override;

    virtual ~EllipticalTube() {}

    EllipticalTube& operator=(const Geometry&) override;

    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    double GetDx() const { return dx_; }
    double GetDy() const { return dy_; }
    double GetDz() const { return dz_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double dx_;
    double dy_;
    double dz_;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::EllipticalTube, 0);
CEREAL_REGISTER_TYPE(siren::geometry::EllipticalTube)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::EllipticalTube);

#endif // SIREN_EllipticalTube_H
