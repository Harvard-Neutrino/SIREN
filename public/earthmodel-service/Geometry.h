
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

#ifndef LI_Geometry_H
#define LI_Geometry_H

#include <iostream>
#include <map>
#include <memory>
#include <math.h>
#include <float.h>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Placement.h"

namespace earthmodel {

class Geometry
{
public:
    static constexpr const double GEOMETRY_PRECISION = 1.e-9;
    struct ParticleLocation {
        enum Enum { InfrontGeometry= 0, InsideGeometry, BehindGeometry };
    };
    struct Intersection {
        double distance;
        int hierarchy;
        bool entering;
        Vector3D position;
        bool operator==(Intersection const & other) const {
            return other.distance == distance and other.hierarchy == hierarchy and other.position == position and other.entering == entering;
        }
    };
    struct IntersectionList {
        Vector3D position;
        Vector3D direction;
        std::vector<Geometry::Intersection> intersections;
        bool operator==(IntersectionList const & other) const {
            return other.position == position and other.direction == direction and other.intersections == intersections;
        }
    };
public:
    Geometry(const std::string);
    Geometry(const std::string, Placement const &);
    Geometry(Placement const &);
    //Geometry(const std::string, const Vector3D position);
    Geometry(const Geometry&);
    //Geometry(const nlohmann::json&);

    /* virtual Geometry* clone() const = 0; // virtual constructor idiom (used for deep copies) */
    virtual std::shared_ptr<const Geometry> create() const = 0;
    virtual void swap(Geometry&);

    virtual ~Geometry(){};

    // Operators
    virtual Geometry& operator=(const Geometry&);
    bool operator==(const Geometry& geometry) const;
    bool operator!=(const Geometry& geometry) const;
    friend std::ostream& operator<<(std::ostream&, Geometry const&);

    // ----------------------------------------------------------------- //
    // Member functions
    // ----------------------------------------------------------------- //

    bool IsInside(const Vector3D& position, const Vector3D& direction) const;

    bool IsInfront(const Vector3D& position, const Vector3D& direction) const;

    bool IsBehind(const Vector3D& position, const Vector3D& direction) const;



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
    std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const;


    /*!
     * Calculates the intersections of a ray with the geometry surface
     */
    std::vector<Intersection> Intersections(Vector3D const & position, Vector3D const & direction) const;

    /*!
     * Calculates the distance to the closest approch to the geometry center
     */
    double DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction) const;

    // void swap(Geometry &geometry);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation(const Vector3D& position, const Vector3D& direction) const;

    //Vector3D GetPosition() const { return position_; }

    std::string GetName() const { return name_; }

    Placement GetPlacement() const { return placement_; }

    void SetPlacement(Placement const & placement) { placement_ = placement; }

    //void SetPosition(const Vector3D& position) { position_ = position; };

protected:
    // Implemented in child classes to be able to use equality operator
    virtual bool compare(const Geometry&) const = 0;
    virtual void print(std::ostream&) const     = 0;
    virtual std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const = 0;
    virtual std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const = 0;

    //Vector3D position_; //!< x,y,z-coordinate of origin ( center of box, cylinder, sphere)

    std::string name_; //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)
    Placement placement_;

    Vector3D LocalToGlobalPosition(Vector3D const & p) const;
    Vector3D LocalToGlobalDirection(Vector3D const & d) const;
    Vector3D GlobalToLocalPosition(Vector3D const & p) const;
    Vector3D GlobalToLocalDirection(Vector3D const & d) const;
};

class Box : public Geometry
{
public:
    Box();
    Box(double x, double y, double z);
    Box(Placement const &);
    Box(Placement const &, double x, double y, double z);
    Box(const Box&);
    //Box(const nlohmann::json& config);

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new Box(*this) ); };
    void swap(Geometry&) override;

    virtual ~Box() {}

    // Operators
    Box& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const override;

    // Getter & Setter
    double GetX() const { return x_; }
    double GetY() const { return y_; }
    double GetZ() const { return z_; }

    void SetX(double x) { x_ = x; };
    void SetY(double y) { y_ = y; };
    void SetZ(double z) { z_ = z; };

private:
    bool compare(const Geometry&) const override;
    void print(std::ostream&) const override;

    double x_; //!< width of box in x-direction
    double y_; //!< width of box in y-direction
    double z_; //!< width of box in z-direction
};

class Cylinder : public Geometry
{
public:
    Cylinder();
    Cylinder(double radius, double inner_radius, double z);
    Cylinder(Placement const &);
    Cylinder(Placement const &, double radius, double inner_radius, double z);
    Cylinder(const Cylinder&);
    //Cylinder(const nlohmann::json& config);

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new Cylinder(*this) ); };
    void swap(Geometry&) override;

    virtual ~Cylinder() {}

    // Operators
    Cylinder& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const override;

    // Getter & Setter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }
    double GetZ() const { return z_; }

    void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
    void SetRadius(double radius) { radius_ = radius; };
    void SetZ(double z) { z_ = z; };

private:
    bool compare(const Geometry&) const override;
    void print(std::ostream&) const override;

    double radius_;       //!< the radius of the sphere/ cylinder
    double inner_radius_; //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
    double z_;            //!< height of box/cylinder
};

class Sphere : public Geometry
{
public:
    Sphere();
    Sphere(double radius, double inner_radius);
    Sphere(Placement const &);
    Sphere(Placement const &, double radius, double inner_radius);
    Sphere(const Sphere&);
    //Sphere(const nlohmann::json& config);

    /* Geometry* clone() const override{ return new Sphere(*this); }; */
    std::shared_ptr<const Geometry> create() const override{ return std::shared_ptr<const Geometry>( new Sphere(*this) ); }
    void swap(Geometry&) override;

    virtual ~Sphere() {}

    // Operators
    Sphere& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const override;

    // Getter & Setter
    double GetInnerRadius() const { return inner_radius_; }
    double GetRadius() const { return radius_; }

    void SetInnerRadius(double inner_radius) { inner_radius_ = inner_radius; };
    void SetRadius(double radius) { radius_ = radius; };

private:
    bool compare(const Geometry&) const override;
    void print(std::ostream&) const override;

    double radius_;       //!< the radius of the sphere/ cylinder
    double inner_radius_; //!< for spherical shells or hollow cylinder (0 for sphere / cylinder)
};

class ExtrPoly : public Geometry
{
public:
    struct ZSection {
        ZSection(double zpos_, double offset_[2], double scale_)
            : zpos(zpos_), offset{offset_[0],offset_[1]}, scale(scale_) {}

        double zpos;
        double scale;
        double offset[2];
        void operator=(ZSection const & other) {
            zpos = other.zpos;
            scale = other.scale;
            offset[0] = other.offset[0];
            offset[1] = other.offset[1];
        }
        friend bool operator==(ZSection const & l, ZSection const & r) {
            return (l.zpos == r.zpos &&
                    l.scale == r.scale &&
                    l.offset[0] == r.offset[0] &&
                    l.offset[1] == r.offset[1]);
        }
    };

public:
    ExtrPoly();
    ExtrPoly(const std::vector<std::vector<double>>& polygon,
            const std::vector<ZSection>& zsections);
    ExtrPoly(Placement const &);
    ExtrPoly(Placement const &, const std::vector<std::vector<double>>& polygon,
            const std::vector<ZSection>& zsections);
    ExtrPoly(const ExtrPoly&);
    //ExtrPoly(const nlohmann::json& config);

    /* Geometry* clone() const override{ return new ExtrPoly(*this); }; */
    std::shared_ptr<const Geometry> create() const override{ return std::shared_ptr<const Geometry>( new ExtrPoly(*this) ); }
    void swap(Geometry&) override;

    virtual ~ExtrPoly() {}

    // Operators
    ExtrPoly& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> ComputeDistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> ComputeIntersections(Vector3D const & position, Vector3D const & direction) const override;

    // Getter & Setter
    std::vector<std::vector<double>> GetPolygon() const { return polygon_; }
    std::vector<ZSection> GetZSections() const { return zsections_; }

    void SetPolygon(std::vector<std::vector<double>> polygon ) { polygon_=polygon; }
    void SetZSections(std::vector<ZSection> zsections) { zsections_=zsections; }

    void ComputeLateralPlanes();

private:
    bool compare(const Geometry&) const override;
    void print(std::ostream&) const override;

    std::vector<std::vector<double>> polygon_; //!< vector of (x,y) pairs denoting vertices of polygon
    std::vector<ZSection> zsections_; //!< vector of z sections describing z extent of polygon
    struct plane { double a,b,c,d; }; // a*x + b*y + c*z + d = 0
    std::vector<plane> planes_;
};

} // namespace earthmodel

namespace earthmodel {
    enum Geometry_Type : int { SPHERE, BOX, CYLINDER, EXTRPOLY};
} // namespace earthmodel

namespace earthmodel {
    const std::array<std::string, 4>  Geometry_Name = { "sphere", "box", "cylinder", "extrpoly"};
} // namespace earthmodel

#endif // LI_Geometry_H

