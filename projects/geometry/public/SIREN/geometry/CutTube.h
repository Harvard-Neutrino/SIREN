#pragma once
#ifndef SIREN_CutTube_H
#define SIREN_CutTube_H

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

class CutTube : public Geometry {
public:
    CutTube();
    CutTube(double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm);
    CutTube(Placement const &);
    CutTube(Placement const &, double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm);
    CutTube(double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm,
            double start_phi, double delta_phi);
    CutTube(Placement const &, double rmin, double rmax, double dz, math::Vector3D low_norm, math::Vector3D high_norm,
            double start_phi, double delta_phi);
    CutTube(const CutTube&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("InnerRadius", rmin_));
            archive(::cereal::make_nvp("OuterRadius", rmax_));
            archive(::cereal::make_nvp("HalfZ", dz_));
            archive(::cereal::make_nvp("LowNorm", low_norm_));
            archive(::cereal::make_nvp("HighNorm", high_norm_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("DeltaPhi", delta_phi_));
            archive(cereal::virtual_base_class<Geometry>(this));
            has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
        } else {
            throw std::runtime_error("CutTube only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new CutTube(*this) ); };
    void swap(Geometry&) override;

    virtual ~CutTube() {}

    // Operators
    CutTube& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter
    double GetRmin() const { return rmin_; }
    double GetRmax() const { return rmax_; }
    double GetDz() const { return dz_; }
    math::Vector3D GetLowNorm() const { return low_norm_; }
    math::Vector3D GetHighNorm() const { return high_norm_; }
    double GetStartPhi() const { return start_phi_; }
    double GetDeltaPhi() const { return delta_phi_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    double rmin_;                //!< inner radius (0 for solid)
    double rmax_;                //!< outer radius
    double dz_;                  //!< half-height along z before cut
    math::Vector3D low_norm_;    //!< outward unit normal of low-z end cap (z < 0)
    math::Vector3D high_norm_;   //!< outward unit normal of high-z end cap (z > 0)
    double start_phi_;           //!< starting azimuthal angle (radians, default 0)
    double delta_phi_;           //!< azimuthal extent (radians, default 2*pi)
    bool has_phi_cut_;           //!< true if delta_phi < 2*pi (precomputed)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::CutTube, 0);
CEREAL_REGISTER_TYPE(siren::geometry::CutTube)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::CutTube);

#endif // SIREN_CutTube_H
