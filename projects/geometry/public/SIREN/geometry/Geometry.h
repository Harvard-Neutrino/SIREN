
/******************************************************************************
 *                                                                            *
 * This file is part of the simulation tool PROPOSAL.                         *
 *                                                                            *
 * Copyright (C) 2017 TU Dortmund University, Department of Physics,          *
 *                    Chair Experimental Physics 5b                           *
 *                                                                            *
 * This software may be modified and distributed under the terms of a         *
 * modified GNU Lesser General Public Licence version 3 (LGPL),               *
 * copied verbatim in the file "LICENSE".                                     *
 *                                                                            *
 * Modifcations to the LGPL License:                                          *
 *                                                                            *
 *      1. The user shall acknowledge the use of PROPOSAL by citing the       *
 *         following reference:                                               *
 *                                                                            *
 *         J.H. Koehne et al.  Comput.Phys.Commun. 184 (2013) 2070-2090 DOI:  *
 *         10.1016/j.cpc.2013.04.001                                          *
 *                                                                            *
 *      2. The user should report any bugs/errors or improvments to the       *
 *         current maintainer of PROPOSAL or open an issue on the             *
 *         GitHub webpage                                                     *
 *                                                                            *
 *         "https://github.com/tudo-astroparticlephysics/PROPOSAL"            *
 *                                                                            *
 ******************************************************************************/

#pragma once
#ifndef SIREN_Geometry_H
#define SIREN_Geometry_H

#include <array>
#include <memory>
#include <string>
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
#include "SIREN/geometry/AABB.h"
#include "SIREN/geometry/Placement.h"
#include <NamedType/named_type.hpp>

namespace siren {
namespace geometry {

// Tagged types for local (body-frame) coordinates.
// Prevents accidentally passing global coordinates to local-frame methods.
using LocalPosition = fluent::NamedType<math::Vector3D, struct LocalPositionTag, fluent::Callable, fluent::Comparable>;
using LocalDirection = fluent::NamedType<math::Vector3D, struct LocalDirectionTag, fluent::Callable, fluent::Comparable>;

class Geometry {
friend cereal::access;
public:
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {};
    static constexpr const double GEOMETRY_PRECISION = 1.e-9;
    struct Intersection {
        double distance;
        int hierarchy;
        bool entering;
        int matID;
        math::Vector3D position;
        bool operator==(Intersection const & other) const {
            return other.distance == distance and other.hierarchy == hierarchy and other.position == position and other.entering == entering;
        }
        template<typename Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::make_nvp("Distance", distance));
                archive(cereal::make_nvp("Hierarchy", hierarchy));
                archive(cereal::make_nvp("Entering", entering));
                archive(cereal::make_nvp("Position", position));
            } else {
                throw std::runtime_error("Intersection only supports version <= 0!");
            }
        }
    };
    struct IntersectionList {
        math::Vector3D position;
        math::Vector3D direction;
        std::vector<Geometry::Intersection> intersections;
        bool operator==(IntersectionList const & other) const {
            return other.position == position and other.direction == direction and other.intersections == intersections;
        }
        template<typename Archive>
        void serialize(Archive & archive, std::uint32_t const version) {
            if(version == 0) {
                archive(cereal::make_nvp("Position", position));
                archive(cereal::make_nvp("Direction", direction));
                archive(cereal::make_nvp("Intersections", intersections));
            } else {
                throw std::runtime_error("IntersectionList only supports version <= 0!");
            }
        }
    };
public:
    Geometry();
    Geometry(const std::string);
    Geometry(const std::string, Placement const &);
    Geometry(Placement const &);
    //Geometry(const std::string, const math::Vector3D position);
    Geometry(const Geometry&);
    //Geometry(const nlohmann::json&);

    /* virtual Geometry* clone() const = 0; // virtual constructor idiom (used for deep copies) */
    virtual std::shared_ptr<Geometry> create() const = 0;
    virtual void swap(Geometry&);

    virtual ~Geometry(){};

    // Operators
    virtual Geometry& operator=(const Geometry&);
    bool operator==(const Geometry& geometry) const;
    bool operator<(const Geometry& geometry) const;
    bool operator!=(const Geometry& geometry) const;
    friend std::ostream& operator<<(std::ostream&, Geometry const&);

    // ----------------------------------------------------------------- //
    // Member functions
    // ----------------------------------------------------------------- //

    /*!
     * Tests whether a point (in global coordinates) is inside this geometry.
     * Uses ComputeIntersections and checks whether the first forward
     * intersection is an exit (meaning the point is inside).
     *
     * The direction parameter controls boundary disambiguation: a point
     * exactly on a surface is "inside" if the ray in that direction
     * exits the volume (i.e. the particle is moving outward from inside).
     * The directionless overload uses (0,0,1) for cases where boundary
     * behavior does not matter.
     */
    bool IsInside(const math::Vector3D& position, const math::Vector3D& direction) const;
    bool IsInside(const math::Vector3D& position) const;

    /*!
     * Calculates the intersections of a ray with the geometry surface.
     * Takes global coordinates; transforms to local internally.
     */
    std::vector<Intersection> Intersections(math::Vector3D const & position, math::Vector3D const & direction) const;

    /*!
     * Calculates intersections from pre-transformed local coordinates.
     * Skips the GlobalToLocal transform. Use with ToLocal() for
     * efficient pre-filtered intersection queries.
     */
    std::vector<Intersection> Intersections(LocalPosition const & position, LocalDirection const & direction) const;

    /*!
     * Transform global coordinates to this geometry's local frame.
     */
    std::pair<LocalPosition, LocalDirection> ToLocal(math::Vector3D const & position, math::Vector3D const & direction) const;

    /*!
     * Calculates the distance to the closest approch to the geometry center
     */
    double DistanceToClosestApproach(const math::Vector3D& position, const math::Vector3D& direction) const;

    // void swap(Geometry &geometry);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    //math::Vector3D GetPosition() const { return position_; }

    std::string GetName() const { return name_; }

    Placement GetPlacement() const { return placement_; }

    void SetPlacement(Placement const & placement) { placement_ = placement; RecomputeWorldAABB(); }

    //void SetPosition(const math::Vector3D& position) { position_ = position; };

    virtual std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const = 0;

    // Returns the axis-aligned bounding box in local coordinates.
    // Each subclass must implement this.
    virtual AABB GetBoundingBox() const = 0;

    // Returns a conservative axis-aligned bounding box in global coordinates.
    // Transforms the local AABB through the Placement (rotates the 8 corners
    // and takes the axis-aligned envelope).
    AABB GetWorldBoundingBox() const;

protected:
    // Implemented in child classes to be able to use equality operator
    virtual bool equal(const Geometry&) const = 0;
    virtual bool less(const Geometry&) const = 0;
    virtual void print(std::ostream&) const     = 0;

    //math::Vector3D position_; //!< x,y,z-coordinate of origin ( center of box, cylinder, sphere)

    void RecomputeWorldAABB();

    std::string name_; //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)
    Placement placement_;
    mutable AABB cached_world_aabb_;
    mutable bool world_aabb_valid_{false};

public:
    math::Vector3D LocalToGlobalPosition(math::Vector3D const & p) const;
    math::Vector3D LocalToGlobalDirection(math::Vector3D const & d) const;
    math::Vector3D GlobalToLocalPosition(math::Vector3D const & p) const;
    math::Vector3D GlobalToLocalDirection(math::Vector3D const & d) const;
};

enum Geometry_Type : int { SPHERE, BOX, CYLINDER, EXTRPOLY, TRIANGULARMESH, CONE, TRD, POLYCONE, POLYHEDRA, BOOLEAN, TORUS, ELTUBE, CUTTUBE, TRAP, ELLIPSOID, PARA};

const std::array<std::string, 16>  Geometry_Name = { "sphere", "box", "cylinder", "extrpoly", "triangularmesh", "cone", "trd", "polycone", "polyhedra", "boolean", "torus", "eltube", "cuttube", "trap", "ellipsoid", "para"};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Geometry, 0);
CEREAL_REGISTER_TYPE(siren::geometry::Geometry)
CEREAL_CLASS_VERSION(siren::geometry::Geometry::Intersection, 0);
CEREAL_CLASS_VERSION(siren::geometry::Geometry::IntersectionList, 0);

CEREAL_FORCE_DYNAMIC_INIT(siren_Geometry);

#endif // SIREN_Geometry_H
