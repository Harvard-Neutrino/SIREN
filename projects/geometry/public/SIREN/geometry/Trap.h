#pragma once
#ifndef SIREN_Trap_H
#define SIREN_Trap_H

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

class Trap : public Geometry {
public:
    Trap();
    Trap(double dz, double theta, double phi,
         double dy1, double dx1, double dx2, double alpha1,
         double dy2, double dx3, double dx4, double alpha2);
    Trap(Placement const &);
    Trap(Placement const &,
         double dz, double theta, double phi,
         double dy1, double dx1, double dx2, double alpha1,
         double dy2, double dx3, double dx4, double alpha2);
    Trap(const Trap&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("HalfZ", dz_));
            archive(::cereal::make_nvp("Theta", theta_));
            archive(::cereal::make_nvp("Phi", phi_));
            archive(::cereal::make_nvp("HalfY1", dy1_));
            archive(::cereal::make_nvp("HalfX1", dx1_));
            archive(::cereal::make_nvp("HalfX2", dx2_));
            archive(::cereal::make_nvp("Alpha1", alpha1_));
            archive(::cereal::make_nvp("HalfY2", dy2_));
            archive(::cereal::make_nvp("HalfX3", dx3_));
            archive(::cereal::make_nvp("HalfX4", dx4_));
            archive(::cereal::make_nvp("Alpha2", alpha2_));
            archive(cereal::virtual_base_class<Geometry>(this));
            ComputeVertices();
            ComputePlanes();
        } else {
            throw std::runtime_error("Trap only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Trap(*this) ); };
    void swap(Geometry&) override;

    virtual ~Trap() {}

    // Operators
    Trap& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter
    double GetDz() const { return dz_; }
    double GetTheta() const { return theta_; }
    double GetPhi() const { return phi_; }
    double GetDy1() const { return dy1_; }
    double GetDx1() const { return dx1_; }
    double GetDx2() const { return dx2_; }
    double GetAlpha1() const { return alpha1_; }
    double GetDy2() const { return dy2_; }
    double GetDx3() const { return dx3_; }
    double GetDx4() const { return dx4_; }
    double GetAlpha2() const { return alpha2_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    void ComputeVertices();
    void ComputePlanes();
    void ComputePlane(const double v0[3], const double v1[3], const double v2[3], double plane[4]) const;

    double dz_;     //!< half-length along z
    double theta_;  //!< polar angle of line joining face centres
    double phi_;    //!< azimuthal angle of that line
    double dy1_;    //!< half-length along y of face at -dz
    double dx1_;    //!< half-length along x of side at y=-dy1 of face at -dz
    double dx2_;    //!< half-length along x of side at y=+dy1 of face at -dz
    double alpha1_; //!< angle wrt y-axis from centre of y-low to y-high at -dz
    double dy2_;    //!< half-length along y of face at +dz
    double dx3_;    //!< half-length along x of side at y=-dy2 of face at +dz
    double dx4_;    //!< half-length along x of side at y=+dy2 of face at +dz
    double alpha2_; //!< angle wrt y-axis from centre of y-low to y-high at +dz

    double vertices_[8][3]; //!< precomputed vertex positions
    double planes_[6][4];   //!< precomputed face planes (nx, ny, nz, d)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Trap, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Trap)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Trap);

#endif // SIREN_Trap_H
