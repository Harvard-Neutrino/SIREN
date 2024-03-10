#pragma once
#ifndef SIREN_GeometryMesh_H
#define SIREN_GeometryMesh_H

#include <memory>
#include <vector>
#include <cstdint>
#include <utility>
#include <iostream>
#include <stdexcept>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/serialization/array.h"

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/MeshBuilder.h"

namespace siren {
namespace geometry {

class TriangularMesh : public Geometry {
public:
    TriangularMesh();
    TriangularMesh(Mesh::TMesh const &);
    TriangularMesh(Placement const &);
    TriangularMesh(Placement const &, Mesh::TMesh const &);
    TriangularMesh(TriangularMesh const &);

    Mesh::VAttribute & GetVertex(Mesh::Vertex v);
    Mesh::EAttribute & GetEdge(Mesh::Edge e);
    Mesh::TAttribute & GetTriangle(Mesh::Triangle t);
    Mesh::VAttribute const & GetVertex(Mesh::Vertex v) const;
    Mesh::EAttribute const & GetEdge(Mesh::Edge v) const;
    Mesh::TAttribute const & GetTriangle(Mesh::Triangle t) const;

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            double data_;
            archive(::cereal::make_nvp("", data_));
            archive(cereal::virtual_base_class<Geometry>(this));
        } else {
            throw std::runtime_error("TriangularMesh only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>( new TriangularMesh(*this) ); };
    void swap(Geometry&) override;

    virtual ~TriangularMesh() {}

    // Operators
    TriangularMesh& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    Mesh::TMesh mesh;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::TriangularMesh, 0);
CEREAL_REGISTER_TYPE(siren::geometry::TriangularMesh)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::TriangularMesh);

#endif // SIREN_GeometryMesh_H

