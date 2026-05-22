#pragma once
#ifndef SIREN_Para_H
#define SIREN_Para_H

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

class Para : public Geometry {
public:
    Para();
    Para(double dx, double dy, double dz, double alpha, double theta, double phi);
    Para(Placement const &);
    Para(Placement const &, double dx, double dy, double dz, double alpha, double theta, double phi);
    Para(const Para&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("HalfX", dx_));
            archive(::cereal::make_nvp("HalfY", dy_));
            archive(::cereal::make_nvp("HalfZ", dz_));
            archive(::cereal::make_nvp("Alpha", alpha_));
            archive(::cereal::make_nvp("Theta", theta_));
            archive(::cereal::make_nvp("Phi", phi_));
            archive(cereal::virtual_base_class<Geometry>(this));
            ComputeFaceData();
        } else {
            throw std::runtime_error("Para only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Para(*this) ); };
    void swap(Geometry&) override;

    virtual ~Para() {}

    // Operators
    Para& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter
    double GetDx() const { return dx_; }
    double GetDy() const { return dy_; }
    double GetDz() const { return dz_; }
    double GetAlpha() const { return alpha_; }
    double GetTheta() const { return theta_; }
    double GetPhi() const { return phi_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;
    void ComputeFaceData();

    double dx_;    //!< half-length along x
    double dy_;    //!< half-length along y
    double dz_;    //!< half-length along z
    double alpha_; //!< angle of y-axis shear (radians)
    double theta_; //!< polar angle of z-axis tilt (radians)
    double phi_;   //!< azimuthal angle of z-axis tilt (radians)

    // Precomputed face normals and distances for the 3 slab pairs
    double face_normals_[3][3]; //!< face_normals_[i][j] = j-th component of i-th normal
    double face_dists_[3];      //!< half-slab distance for each pair
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Para, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Para)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Para);

#endif // SIREN_Para_H
