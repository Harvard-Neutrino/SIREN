#pragma once
#ifndef LI_GeometryMesh_H
#define LI_GeometryMesh_H

#include <iostream>
#include <map>
#include <memory>
#include <math.h>
#include <float.h>

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
#include "serialization/array.h"

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Placement.h"
#include "earthmodel-service/Geometry.h"

namespace earthmodel {

class TriangularMesh : public Geometry {
public:

    TriangularMesh();
    TriangularMesh(TMesh const &);
    TriangularMesh(Placement const &);
    TriangularMesh(Placement const &, TMesh const &);
    TriangularMesh(TriangularMesh const &);

    VAttribute & GetVertex(Vertex v);
    EAttribute & GetEdge(Edge v);
    TAttribute & GetTriangle(Triangle t);
    VAttribute const & GetVertex(Vertex v) const;
    EAttribute const & GetEdge(Edge v) const;
    TAttribute const & GetTriangle(Triangle t) const;

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

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new TriangularMesh(*this) ); };
    void swap(Geometry&) override;

    virtual ~TriangularMesh() {}

    // Operators
    TriangularMesh& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const override;

protected:
    virtual bool equal(const Geometry&) const override;
    virtual bool less(const Geometry&) const override;
private:
    void print(std::ostream&) const override;

    TMesh mesh;
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::TriangularMesh, 0);
CEREAL_REGISTER_TYPE(earthmodel::TriangularMesh)
CEREAL_REGISTER_POLYMORPHIC_RELATION(earthmodel::Geometry, earthmodel::TriangularMesh);

#endif // LI_GeometryMesh_H

