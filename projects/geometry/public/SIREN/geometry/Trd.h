#pragma once
#ifndef SIREN_Trd_H
#define SIREN_Trd_H

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

class Trd : public Geometry {
public:
    Trd();
    Trd(double dx1, double dx2, double dy1, double dy2, double dz);
    Trd(Placement const &);
    Trd(Placement const &, double dx1, double dx2, double dy1, double dy2, double dz);
    Trd(const Trd&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("HalfX1", dx1_));
            archive(::cereal::make_nvp("HalfX2", dx2_));
            archive(::cereal::make_nvp("HalfY1", dy1_));
            archive(::cereal::make_nvp("HalfY2", dy2_));
            archive(::cereal::make_nvp("HalfZ", dz_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Trd only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Trd(*this) ); };
    void swap(Geometry&) override;

    virtual ~Trd() {}

    // Operators
    Trd& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    double GetDx1() const { return dx1_; }
    double GetDx2() const { return dx2_; }
    double GetDy1() const { return dy1_; }
    double GetDy2() const { return dy2_; }
    double GetDz() const { return dz_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double dx1_; //!< half-width in x at z = -dz
    double dx2_; //!< half-width in x at z = +dz
    double dy1_; //!< half-width in y at z = -dz
    double dy2_; //!< half-width in y at z = +dz
    double dz_;  //!< half-height along z-axis
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Trd, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Trd)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Trd);

#endif // SIREN_Trd_H
