#pragma once
#ifndef SIREN_GenericPolycone_H
#define SIREN_GenericPolycone_H

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

class GenericPolycone : public Geometry {
public:
    GenericPolycone();
    GenericPolycone(std::vector<double> const & r,
                    std::vector<double> const & z);
    GenericPolycone(Placement const &,
                    std::vector<double> const & r,
                    std::vector<double> const & z);
    GenericPolycone(std::vector<double> const & r,
                    std::vector<double> const & z,
                    double start_phi, double delta_phi);
    GenericPolycone(Placement const &,
                    std::vector<double> const & r,
                    std::vector<double> const & z,
                    double start_phi, double delta_phi);
    GenericPolycone(const GenericPolycone&);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("R", r_));
            archive(::cereal::make_nvp("Z", z_));
            archive(::cereal::make_nvp("StartPhi", start_phi_));
            archive(::cereal::make_nvp("DeltaPhi", delta_phi_));
            archive(cereal::virtual_base_class<Geometry>(this));
            has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
            ComputeWinding();
        } else {
            throw std::runtime_error("GenericPolycone only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>(new GenericPolycone(*this)); }
    void swap(Geometry&) override;

    virtual ~GenericPolycone() {}

    GenericPolycone& operator=(const Geometry&) override;

    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    std::vector<double> const & GetR() const { return r_; }
    std::vector<double> const & GetZ() const { return z_; }
    double GetStartPhi() const { return start_phi_; }
    double GetDeltaPhi() const { return delta_phi_; }

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;
    void validate() const;
    void ComputeWinding();

    std::vector<double> r_;
    std::vector<double> z_;
    double start_phi_;
    double delta_phi_;
    bool has_phi_cut_;
    int winding_;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::GenericPolycone, 0);
CEREAL_REGISTER_TYPE(siren::geometry::GenericPolycone)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::GenericPolycone);

#endif // SIREN_GenericPolycone_H
