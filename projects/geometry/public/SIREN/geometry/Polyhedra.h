#pragma once
#ifndef SIREN_Polyhedra_H
#define SIREN_Polyhedra_H

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

class Polyhedra : public Geometry {
public:
    Polyhedra();
    Polyhedra(int num_sides,
              double start_phi,
              std::vector<double> const & z_planes,
              std::vector<double> const & rmin,
              std::vector<double> const & rmax);
    Polyhedra(Placement const &);
    Polyhedra(Placement const &,
              int num_sides,
              double start_phi,
              std::vector<double> const & z_planes,
              std::vector<double> const & rmin,
              std::vector<double> const & rmax);
    Polyhedra(const Polyhedra&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("NumSides", num_sides_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("ZPlanes", z_planes_));
            archive(::cereal::make_nvp("Rmin", rmin_));
            archive(::cereal::make_nvp("Rmax", rmax_));
            archive(cereal::virtual_base_class<Geometry>(this));
            precompute_trig();
        } else {
            throw std::runtime_error("Polyhedra only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new Polyhedra(*this) ); };
    void swap(Geometry&) override;

    virtual ~Polyhedra() {}

    // Operators
    Polyhedra& operator=(const Geometry&) override;

    // Methods
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    // Getter & Setter
    int GetNumSides() const { return num_sides_; }
    double GetStartPhi() const { return start_phi_; }
    std::vector<double> const & GetZPlanes() const { return z_planes_; }
    std::vector<double> const & GetRmin() const { return rmin_; }
    std::vector<double> const & GetRmax() const { return rmax_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    // Validate parameters
    void validate() const;

    // Precompute cached trig arrays from num_sides_ and start_phi_
    void precompute_trig();

    int num_sides_;            //!< number of polygon sides (>= 3)
    double start_phi_;         //!< starting angle of first polygon edge (radians)
    std::vector<double> z_planes_; //!< z-positions of each plane (ascending)
    std::vector<double> rmin_;     //!< inner radius at each z-plane (to edge center)
    std::vector<double> rmax_;     //!< outer radius at each z-plane (to edge center)

    // Precomputed cos/sin for polygon vertex angles.
    // cos_phi_[k] = cos(start_phi_ + k * dphi), k = 0..num_sides_
    // Entry [num_sides_] wraps around to [0].
    std::vector<double> cos_phi_;
    std::vector<double> sin_phi_;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Polyhedra, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Polyhedra)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::Polyhedra);

#endif // SIREN_Polyhedra_H
