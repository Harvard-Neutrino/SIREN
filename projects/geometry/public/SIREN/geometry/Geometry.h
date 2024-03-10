
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
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

class Geometry {
friend cereal::access;
public:
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {};
    static constexpr const double GEOMETRY_PRECISION = 1.e-9;
    struct ParticleLocation {
        enum Enum { InfrontGeometry= 0, InsideGeometry, BehindGeometry };
    };
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

    bool IsInside(const math::Vector3D& position, const math::Vector3D& direction) const;

    bool IsInfront(const math::Vector3D& position, const math::Vector3D& direction) const;

    bool IsBehind(const math::Vector3D& position, const math::Vector3D& direction) const;



    /*!
     * This function calculates the distance of the particle position
     * to the border of the geometry in direction of the particle trajectory.
     * If the particle trajectory does not have an intersection with the geometry
     * (-1 /-1) is returned
     * If the particle trajectory has two intersections (dist_1 /dist_2) is returned
     * If the particle has one intersection (dist_1 /-1) is returned
     * (one intersection means one intersection in direction of the particle trajectory
     * and one in the opposite direction. Cause we are not interested in this one. it is set to -1)
     * Note: If the particle is on the geometry border this is not treated as an intersection
     * A particle on the geometry border which moves inside has one intersection,
     * a particle on the geometry border which moves outside has no intersection.
     * Distances smaller then GEOMETRY_PRECISION (1e-9) are also set to -1
     */
    std::pair<double, double> DistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const;


    /*!
     * Calculates the intersections of a ray with the geometry surface
     */
    std::vector<Intersection> Intersections(math::Vector3D const & position, math::Vector3D const & direction) const;

    /*!
     * Calculates the distance to the closest approch to the geometry center
     */
    double DistanceToClosestApproach(const math::Vector3D& position, const math::Vector3D& direction) const;

    // void swap(Geometry &geometry);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation(const math::Vector3D& position, const math::Vector3D& direction) const;

    //math::Vector3D GetPosition() const { return position_; }

    std::string GetName() const { return name_; }

    Placement GetPlacement() const { return placement_; }

    void SetPlacement(Placement const & placement) { placement_ = placement; }

    //void SetPosition(const math::Vector3D& position) { position_ = position; };

    virtual std::vector<Intersection> ComputeIntersections(math::Vector3D const & position, math::Vector3D const & direction) const = 0;

protected:
    // Implemented in child classes to be able to use equality operator
    virtual bool equal(const Geometry&) const = 0;
    virtual bool less(const Geometry&) const = 0;
    virtual void print(std::ostream&) const     = 0;
    virtual std::pair<double, double> ComputeDistanceToBorder(const math::Vector3D& position, const math::Vector3D& direction) const = 0;

    //math::Vector3D position_; //!< x,y,z-coordinate of origin ( center of box, cylinder, sphere)

    std::string name_; //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)
    Placement placement_;

public:
    math::Vector3D LocalToGlobalPosition(math::Vector3D const & p) const;
    math::Vector3D LocalToGlobalDirection(math::Vector3D const & d) const;
    math::Vector3D GlobalToLocalPosition(math::Vector3D const & p) const;
    math::Vector3D GlobalToLocalDirection(math::Vector3D const & d) const;
};

enum Geometry_Type : int { SPHERE, BOX, CYLINDER, EXTRPOLY, TRIANGULARMESH};

const std::array<std::string, 5>  Geometry_Name = { "sphere", "box", "cylinder", "extrpoly", "triangularmesh"};

} // namespace geometry
} // namespace siren

CEREAL_CLASS_VERSION(siren::geometry::Geometry, 0);
CEREAL_CLASS_VERSION(siren::geometry::Geometry::Intersection, 0);
CEREAL_CLASS_VERSION(siren::geometry::Geometry::IntersectionList, 0);

#endif // SIREN_Geometry_H
