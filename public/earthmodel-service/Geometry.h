
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

#include "earthmodel-service/Vector3D.h"

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
        Vector3D position;
        bool entering;
    };
public:
    Geometry(const std::string);
    Geometry(const std::string, const Vector3D position);
    Geometry(const std::string, const Vector3D position, int hierarchy);
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
    virtual std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const = 0;


    /*!
     * Calculates the intersections of a ray with the geometry surface
     */
    virtual std::vector<Intersection> Intersections(Vector3D const & position, Vector3D const & direction) const = 0;

    /*!
     * Calculates the distance to the closest approch to the geometry center
     */
    double DistanceToClosestApproach(const Vector3D& position, const Vector3D& direction) const;

    // void swap(Geometry &geometry);

    // ----------------------------------------------------------------- //
    // Getter & Setter
    // ----------------------------------------------------------------- //

    ParticleLocation::Enum GetLocation(const Vector3D& position, const Vector3D& direction) const;

    Vector3D GetPosition() const { return position_; }

    std::string GetName() const { return name_; }

    unsigned int GetHierarchy() const { return hierarchy_; }

    void SetPosition(const Vector3D& position) { position_ = position; };

    void SetHierarchy(unsigned int hierarchy) { hierarchy_ = hierarchy; };

protected:
    // Implemented in child classes to be able to use equality operator
    virtual bool compare(const Geometry&) const = 0;
    virtual void print(std::ostream&) const     = 0;

    Vector3D position_; //!< x,y,z-coordinate of origin ( center of box, cylinder, sphere)

    std::string name_; //!< "box" , "cylinder" , "sphere" (sphere and cylinder might be hollow)

    unsigned int hierarchy_; //!< adds a hierarchy of geometry objects to allow crossing geometries
};

class Box : public Geometry
{
public:
    Box();
    Box(const Vector3D position, double x, double y, double z);
    Box(const Vector3D position, double x, double y, double z, int hierarchy);
    Box(const Box&);
    //Box(const nlohmann::json& config);

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new Box(*this) ); };
    void swap(Geometry&) override;

    virtual ~Box() {}

    // Operators
    Box& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> Intersections(Vector3D const & position, Vector3D const & direction) const override;

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
    Cylinder(const Vector3D position, double radius, double inner_radius, double z);
    Cylinder(const Vector3D position, double radius, double inner_radius, double z, int hierarchy);
    Cylinder(const Cylinder&);
    //Cylinder(const nlohmann::json& config);

    std::shared_ptr<const Geometry> create() const override { return std::shared_ptr<const Geometry>( new Cylinder(*this) ); };
    void swap(Geometry&) override;

    virtual ~Cylinder() {}

    // Operators
    Cylinder& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> Intersections(Vector3D const & position, Vector3D const & direction) const override;

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
    Sphere(const Vector3D position, double radius, double inner_radius);
    Sphere(const Vector3D position, double radius, double inner_radius, int hierarchy);
    Sphere(const Sphere&);
    //Sphere(const nlohmann::json& config);

    /* Geometry* clone() const override{ return new Sphere(*this); }; */
    std::shared_ptr<const Geometry> create() const override{ return std::shared_ptr<const Geometry>( new Sphere(*this) ); }
    void swap(Geometry&) override;

    virtual ~Sphere() {}

    // Operators
    Sphere& operator=(const Geometry&) override;

    // Methods
    std::pair<double, double> DistanceToBorder(const Vector3D& position, const Vector3D& direction) const override;
    std::vector<Intersection> Intersections(Vector3D const & position, Vector3D const & direction) const override;

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

} // namespace earthmodel

namespace earthmodel {
    enum Geometry_Type : int { SPHERE, BOX, CYLINDER };
} // namespace earthmodel

namespace earthmodel {
    const std::array<std::string, 3>  Geometry_Name = { "sphere", "box", "cylinder" };
} // namespace earthmodel

#endif // LI_Geometry_H

