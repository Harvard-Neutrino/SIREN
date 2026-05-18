#pragma once
#ifndef SIREN_Polycone_H
#define SIREN_Polycone_H

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
#include <cereal/types/vector.hpp>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Geometry.h"

namespace siren {
namespace geometry {

class Polycone : public Geometry {
public:
    Polycone();
    Polycone(std::vector<double> const & z_planes,
             std::vector<double> const & rmin,
             std::vector<double> const & rmax);
    Polycone(Placement const &);
    Polycone(Placement const &,
             std::vector<double> const & z_planes,
             std::vector<double> const & rmin,
             std::vector<double> const & rmax);
    Polycone(std::vector<double> const & z_planes,
             std::vector<double> const & rmin,
             std::vector<double> const & rmax,
             double start_phi, double delta_phi);
    Polycone(Placement const &,
             std::vector<double> const & z_planes,
             std::vector<double> const & rmin,
             std::vector<double> const & rmax,
             double start_phi, double delta_phi);
    Polycone(const Polycone&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("ZPlanes", z_planes_));
            archive(::cereal::make_nvp("Rmin", rmin_));
            archive(::cereal::make_nvp("Rmax", rmax_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("DeltaPhi", delta_phi_));
            archive(cereal::virtual_base_class<Geometry>(this));
            has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
        } else {
            throw std::runtime_error("Polycone only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Polycone(*this) ); };
    void swap(Geometry&) override;

    virtual ~Polycone() {}

    // Operators
    Polycone& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    std::vector<double> const & GetZPlanes() const { return z_planes_; }
    std::vector<double> const & GetRmin() const { return rmin_; }
    std::vector<double> const & GetRmax() const { return rmax_; }
    double GetStartPhi() const { return start_phi_; }
    double GetDeltaPhi() const { return delta_phi_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    // Validate parameters and normalize z-plane order to ascending.
    void validate();

    std::vector<double> z_planes_; //!< z-positions of each plane (ascending)
    std::vector<double> rmin_;     //!< inner radius at each z-plane
    std::vector<double> rmax_;     //!< outer radius at each z-plane
    double start_phi_;    //!< starting azimuthal angle (radians, default 0)
    double delta_phi_;    //!< azimuthal extent (radians, default 2*pi)
    bool has_phi_cut_;    //!< true if delta_phi < 2*pi (precomputed)
};


} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Polycone, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Polycone)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Polycone);

#endif // SIREN_Polycone_H
