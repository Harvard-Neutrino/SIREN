#pragma once
#ifndef SIREN_GeometryMesh_H
#define SIREN_GeometryMesh_H

#include <memory>
#include <vector>
#include <array>
#include <cstdint>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <string>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "SIREN/serialization/array.h"

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Geometry.h"

namespace siren {
namespace geometry {

class TriangularMesh : public Geometry {
public:
    using Triangle = std::array<std::array<double, 3>, 3>;

    TriangularMesh();
    TriangularMesh(std::vector<std::array<math::Vector3D, 3>> const & triangles);
    TriangularMesh(Placement const &, std::vector<std::array<math::Vector3D, 3>> const & triangles);
    TriangularMesh(Placement const &);
    TriangularMesh(TriangularMesh const &);

    template<typename Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("Triangles", triangles_));
            archive(cereal::virtual_base_class<Geometry>(this));
            BuildBVH();
        } else {
            throw std::runtime_error("TriangularMesh only supports version <= 0!");
        }
    }

    std::shared_ptr<Geometry> create() const override { return std::shared_ptr<Geometry>(new TriangularMesh(*this)); }
    void swap(Geometry&) override;

    virtual ~TriangularMesh() {}

    TriangularMesh& operator=(const Geometry&) override;

    std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const override;
    AABB GetBoundingBox() const override;

    size_t TriangleCount() const { return triangles_.size(); }

    // Validate that the mesh is a closed 2-manifold surface.
    // Returns empty string if valid, or a description of the defect.
    std::string ValidateClosed() const;

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;
    void BuildBVH();

    // Triangle storage: triangles_[i] = {{v0x,v0y,v0z}, {v1x,v1y,v1z}, {v2x,v2y,v2z}}
    std::vector<Triangle> triangles_;

    // BVH flat-array acceleration structure
    struct BVHNode {
        float bounds[6]; // min_x, min_y, min_z, max_x, max_y, max_z
        uint32_t offset; // leaf: first index into tri_indices_; internal: right child index
        uint16_t count;  // 0 = internal; >0 = leaf with count triangles
        uint16_t axis;   // split axis for internal nodes
    };
    std::vector<BVHNode> bvh_nodes_;
    std::vector<uint32_t> tri_indices_;
};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::TriangularMesh, 0);
CEREAL_REGISTER_TYPE(siren::geometry::TriangularMesh)
CEREAL_REGISTER_POLYMORPHIC_RELATION(siren::geometry::Geometry, siren::geometry::TriangularMesh);

#endif // SIREN_GeometryMesh_H
