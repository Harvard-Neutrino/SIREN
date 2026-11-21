#pragma once
#ifndef SIREN_Ellipsoid_H
#define SIREN_Ellipsoid_H

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

class Ellipsoid : public Geometry {
public:
    Ellipsoid();
    Ellipsoid(double ax, double by, double cz);
    Ellipsoid(double ax, double by, double cz, double zcut1, double zcut2);
    Ellipsoid(Placement const &);
    Ellipsoid(Placement const &, double ax, double by, double cz);
    Ellipsoid(Placement const &, double ax, double by, double cz, double zcut1, double zcut2);
    Ellipsoid(const Ellipsoid&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("SemiAxisX", ax_));
            archive(::cereal::make_nvp("SemiAxisY", by_));
            archive(::cereal::make_nvp("SemiAxisZ", cz_));
            archive(::cereal::make_nvp("Zcut1", zcut1_));
            archive(::cereal::make_nvp("Zcut2", zcut2_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("Ellipsoid only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Ellipsoid(*this) ); };
    void swap(Geometry&) override;

    virtual ~Ellipsoid() {}

    // Operators
    Ellipsoid& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter
    double GetAx() const { return ax_; }
    double GetBy() const { return by_; }
    double GetCz() const { return cz_; }
    double GetZcut1() const { return zcut1_; }
    double GetZcut2() const { return zcut2_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double ax_;    //!< semi-axis along x
    double by_;    //!< semi-axis along y
    double cz_;    //!< semi-axis along z
    double zcut1_; //!< lower z-cut plane
    double zcut2_; //!< upper z-cut plane
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Ellipsoid, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Ellipsoid)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Ellipsoid);

#endif // SIREN_Ellipsoid_H
