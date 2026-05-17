#include <iostream>
#include <stdlib.h>
#include <random>
#include <math.h>
#include <vector>
#include <utility>

#include <gtest/gtest.h>

#include "SIREN/math/Vector3D.h"

#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/Trap.h"
#include "SIREN/geometry/Trd.h"

using namespace siren::geometry;
using namespace siren::math;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

// Replaces the removed DistanceToBorder method.
// Returns pair<double,double> with the same semantics as the original.
std::pair<double, double> DistanceToBorder(
        Geometry const & geo,
        Vector3D const & position,
        Vector3D const & direction) {
    std::vector<Geometry::Intersection> isects = geo.Intersections(position, direction);
    std::pair<double, double> result(-1, -1);
    int count = 0;
    for(auto const & hit : isects) {
        if(hit.distance > Geometry::GEOMETRY_PRECISION) {
            if(count == 0) {
                result.first = hit.distance;
                if(!hit.entering) break;
                count++;
            } else {
                result.second = hit.distance;
                break;
            }
        }
    }
    return result;
}

TEST(Comparison, Comparison_equal)
{
    Sphere A;
    Sphere B;
    EXPECT_TRUE(A == B);

    Sphere* C = new Sphere();
    Sphere* D = new Sphere();
    EXPECT_TRUE(*C == *D);

    Geometry* E = new Sphere();
    Geometry* F = new Sphere();
    EXPECT_TRUE(*E == *F);

    delete C;
    delete D;
    delete E;
    delete F;
}

TEST(Comparison, Comparison_not_equal)
{
    Sphere A;
    Box B;
    EXPECT_TRUE(A != B);

    Sphere* C = new Sphere();
    Box* D    = new Box();
    EXPECT_TRUE(*C != *D);

    Geometry* E = new Sphere();
    Geometry* F = new Box();
    EXPECT_TRUE(*E != *F);

    delete C;
    delete D;
    delete E;
    delete F;
}

TEST(Assignment, Copyconstructor)
{
    Sphere A;
    Sphere B = A;

    EXPECT_TRUE(A == B);

    Geometry* C = new Sphere();
    Geometry* D = new Box();

    *D = *C;

    EXPECT_FALSE(*C == *D);

    Geometry* E = new Sphere(Vector3D(1.0, 0.0, 0.0), 20.0, 10.0);

    *C = *E;

    EXPECT_TRUE(*C == *E);
}

TEST(Assignment, Copyconstructor2)
{
    Sphere A;
    Sphere B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Operator)
{
    Sphere A;
    Sphere B(Vector3D(0.0, 0.0, 0.0), 2.0, 1.0);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    Sphere A;
    Sphere B;
    EXPECT_TRUE(A == B);
    Geometry* C = new Sphere(Vector3D(1.0, 2.0, 3.0), 4.0, 3.0);
    Geometry* D = new Sphere(Vector3D(1.0, 2.0, 3.0), 4.0, 3.0);

    EXPECT_TRUE(*C == *D);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

TEST(DistanceToClosestApproach, Method)
{
    Sphere A;
    Vector3D position = Vector3D(1, -1, 0);
    Vector3D direction = Vector3D(0, 1, 0);
    direction.CalculateSphericalCoordinates();

    double distance_closest_approach = A.DistanceToClosestApproach(position, direction);

    ASSERT_NEAR(distance_closest_approach, 1., 1e-9);
}

TEST(IsInside, Box)
{

    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;

    double rnd_theta;
    double rnd_phi;

    double width_x = 10;
    double width_y = 10;
    double height  = 10;

    double big_width_x = 4 * width_x;
    double big_width_y = 4 * width_y;
    double big_height  = 4 * height;

    Vector3D position_geometry(0, 0, 0);

    int is_inside  = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {

        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomDouble();
        rnd_y0 = RandomDouble();
        rnd_z0 = RandomDouble();

        position_geometry.SetCartesianCoordinates((2 * rnd_x0 - 1) * 0.5 * (big_width_x - width_x),
                                                  (2 * rnd_y0 - 1) * 0.5 * (big_width_y - width_y),
                                                  (2 * rnd_z0 - 1) * 0.5 * (big_height - height));


        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        Box A(position_geometry, width_x, width_y, height);
        volumia_ratio = width_x * width_y * height / (big_width_x * big_width_y * big_height);
        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomDouble();
            rnd_y = RandomDouble();
            rnd_z = RandomDouble();

            rnd_theta = RandomDouble();
            rnd_phi   = RandomDouble();

            particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
            particle_direction.CalculateCartesianFromSpherical();

            particle_position.SetCartesianCoordinates((2 * rnd_x - 1) * 0.5 * big_width_x,
                                                      (2 * rnd_y - 1) * 0.5 * big_width_y,
                                                      (2 * rnd_z - 1) * 0.5 * big_height);

            // if this constraints are true the particle is inside the box geometry
            if (particle_position.GetX() > position_geometry.GetX() - 0.5 * width_x &&
                particle_position.GetX() < position_geometry.GetX() + 0.5 * width_x &&
                particle_position.GetY() > position_geometry.GetY() - 0.5 * width_y &&
                particle_position.GetY() < position_geometry.GetY() + 0.5 * width_y &&
                particle_position.GetZ() > position_geometry.GetZ() - 0.5 * height &&
                particle_position.GetZ() < position_geometry.GetZ() + 0.5 * height)
            {
                is_inside++;
                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
    }
    // Check what happens if particles are on the border of the box

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Box A(Vector3D(0, 0, 0), width_x, width_y, height);

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle_position.SetCartesianCoordinates(0, 0, 0.5 * height);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        if (particle_direction.GetTheta() < 0.5 * M_PI) {
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetTheta() > 0.5 * M_PI) {
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        }
    }

    // Make this test for every surface of the box

    // bottom
    particle_position.SetCartesianCoordinates(0, 0, -0.5 * height);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        if (particle_direction.GetTheta() > 0.5 * M_PI) {
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetTheta() < 0.5 * M_PI) {
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        }
    }

    // Surface in positiv x direction
    particle_position.SetCartesianCoordinates(0.5 * width_x, 0, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // phi = 0 is in positive x direction
        if (particle_direction.GetPhi() < 0.5 * M_PI || particle_direction.GetPhi() > 1.5 * M_PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in negativ x direction
    particle_position.SetCartesianCoordinates(-0.5 * width_x, 0, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // phi = 0 is in positive x direction
        if (particle_direction.GetPhi() < 0.5 * M_PI || particle_direction.GetPhi() > 1.5 * M_PI)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in positiv y direction
    particle_position.SetCartesianCoordinates(0, 0.5 * width_y, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // phi = 0 is in positive x direction
        if (particle_direction.GetPhi() < M_PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface in negativ y direction
    particle_position.SetCartesianCoordinates(0, -0.5 * width_y, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // phi = 0 is in positive x direction
        if (particle_direction.GetPhi() < M_PI)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }

    // For completness check one corner
    particle_position.SetCartesianCoordinates(0.5 * width_x, 0.5 * width_y, 0.5 * height);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        if (particle_direction.GetTheta() < 0.5 * M_PI || particle_direction.GetPhi() < M_PI ||
            particle_direction.GetPhi() > 1.5 * M_PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
}

TEST(IsInside, Cylinder)
{

    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;

    double rnd_theta;
    double rnd_phi;
    double rnd_inner_radius;

    double radius       = 10;
    double inner_radius = 0;
    double height       = 10;

    double big_width_x = 4 * radius;
    double big_width_y = 4 * radius;
    double big_height  = 4 * height;

    Vector3D position_geometry(0, 0, 0);

    int is_inside  = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {

        // Chose the origin of the cylinder-geometry
        // This cylinder should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomDouble();
        rnd_y0 = RandomDouble();
        rnd_z0 = RandomDouble();

        position_geometry.SetCartesianCoordinates((2 * rnd_x0 - 1) * (0.5 * big_width_x - radius),
                                                  (2 * rnd_y0 - 1) * (0.5 * big_width_y - radius),
                                                  (2 * rnd_z0 - 1) * 0.5 * (big_height - height));

        // position_geometry.SetCartesianCoordinates(0,0,0);

        rnd_inner_radius = RandomDouble();

        inner_radius = radius * rnd_inner_radius;

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        Cylinder A(position_geometry, radius, inner_radius, height);

        volumia_ratio =
            height * M_PI * (pow(radius, 2) - pow(inner_radius, 2)) / (big_width_x * big_width_y * big_height);
        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomDouble();
            rnd_y = RandomDouble();
            rnd_z = RandomDouble();

            rnd_theta = RandomDouble();
            rnd_phi   = RandomDouble();

            particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
            particle_direction.CalculateCartesianFromSpherical();

            particle_position.SetCartesianCoordinates((2 * rnd_x - 1) * 0.5 * big_width_x,
                                                      (2 * rnd_y - 1) * 0.5 * big_width_y,
                                                      (2 * rnd_z - 1) * 0.5 * big_height);

            // if this constraints are true the particle is inside the cylinder geometry
            if (sqrt(pow((particle_position.GetX() - position_geometry.GetX()), 2) +
                     pow((particle_position.GetY() - position_geometry.GetY()), 2)) < radius &&
                sqrt(pow((particle_position.GetX() - position_geometry.GetX()), 2) +
                     pow((particle_position.GetY() - position_geometry.GetY()), 2)) > inner_radius &&
                particle_position.GetZ() > position_geometry.GetZ() - 0.5 * height &&
                particle_position.GetZ() < position_geometry.GetZ() + 0.5 * height)
            {

                is_inside++;

                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
    }

    // Test borders

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Sphere B(Vector3D(0, 0, 0), radius, 0);

    double cos;

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomDouble();

        particle_position.SetCartesianCoordinates(radius * rnd_x, radius * sqrt(1 - rnd_x * rnd_x), 0);

        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // cosine of angle between direction vector and position vector
        cos = -scalar_product(particle_position, particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_TRUE(B.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(B.IsInside(particle_position, particle_direction));
    }

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle_position.SetCartesianCoordinates(0, 0, 0.5 * height);

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Cylinder C(Vector3D(0, 0, 0), radius, 0, height);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();
        // Computer precision controll
        if (particle_position.GetX() * particle_position.GetX() + particle_position.GetY() * particle_position.GetY() -
                inner_radius * inner_radius ==
            0)
        {
            if (particle_direction.GetTheta() < M_PI / 2.) {
                EXPECT_FALSE(C.IsInside(particle_position, particle_direction));
            }
            if (particle_direction.GetTheta() > M_PI / 2.) {
                EXPECT_TRUE(C.IsInside(particle_position, particle_direction));
            }
        }
    }

    // Make this test for every surface of the box

    // bottom
    particle_position.SetCartesianCoordinates(0, 0, -0.5 * height);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        if (particle_direction.GetTheta() > M_PI / 2.) {
            EXPECT_FALSE(C.IsInside(particle_position, particle_direction));
        }
        if (particle_direction.GetTheta() < M_PI / 2.) {
            EXPECT_TRUE(C.IsInside(particle_position, particle_direction));
        }
    }

    // Test inner border
    inner_radius = 5;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Cylinder A(Vector3D(0, 0, 0), radius, inner_radius, height);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomDouble();

        particle_position.SetCartesianCoordinates(inner_radius * rnd_x, inner_radius * sqrt(1 - rnd_x * rnd_x), 0);

        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // cosine of angle between direction vector and position vector
        cos = -scalar_product(particle_position, particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
}

TEST(IsInside, Sphere)
{

    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double rnd_x;
    double rnd_y;
    double rnd_z;

    double rnd_x0;
    double rnd_y0;
    double rnd_z0;
    double rnd_inner_radius;

    double rnd_theta;
    double rnd_phi;

    double radius       = 10;
    double inner_radius = 0;

    double big_width_x = 4 * radius;
    double big_width_y = 4 * radius;
    double big_height  = 4 * radius;

    Vector3D position_geometry(0, 0, 0);

    int is_inside  = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e6;
    int number_volumina  = 1e1;

    for (int i = 0; i < number_volumina; i++)
    {
        // Chose the origin of the box-geometry
        // This box should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomDouble();
        rnd_y0 = RandomDouble();
        rnd_z0 = RandomDouble();

        rnd_inner_radius = RandomDouble();

        position_geometry.SetCartesianCoordinates((2 * rnd_x0 - 1) * (0.5 * big_width_x - radius),
                                                  (2 * rnd_y0 - 1) * (0.5 * big_width_y - radius),
                                                  (2 * rnd_z0 - 1) * (0.5 * big_height - radius));

        inner_radius = radius * rnd_inner_radius;

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        Sphere A(position_geometry, radius, inner_radius);

        volumia_ratio =
            (4. / 3. * M_PI * (pow(radius, 3) - pow(inner_radius, 3))) / (big_width_x * big_width_y * big_height);

        for (int j = 0; j < number_particles; j++)
        {

            // Chose particle location and angle
            rnd_x = RandomDouble();
            rnd_y = RandomDouble();
            rnd_z = RandomDouble();

            rnd_theta = RandomDouble();
            rnd_phi   = RandomDouble();

            particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
            particle_direction.CalculateCartesianFromSpherical();

            particle_position.SetCartesianCoordinates((2 * rnd_x - 1) * 0.5 * big_width_x,
                                                      (2 * rnd_y - 1) * 0.5 * big_width_y,
                                                      (2 * rnd_z - 1) * 0.5 * big_height);

            if ((particle_position - position_geometry).magnitude() < radius &&
                (particle_position - position_geometry).magnitude() > inner_radius)
            {
                is_inside++;
                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } else
            {
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
    }

    // Test borders

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Sphere A(Vector3D(0, 0, 0), radius, 0);

    double cos;

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomDouble();

        particle_position.SetCartesianCoordinates(radius * rnd_x, radius * sqrt(1 - rnd_x * rnd_x), 0);

        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // cosine of angle between direction vector and position vector
        cos = -scalar_product(particle_position, particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
    }

    // Test inner border
    inner_radius = 5;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Sphere B(Vector3D(0, 0, 0), radius, inner_radius);

    for (int i = 0; i < 1e4; i++)
    {
        rnd_x = RandomDouble();

        particle_position.SetCartesianCoordinates(inner_radius * rnd_x, inner_radius * sqrt(1 - rnd_x * rnd_x), 0);

        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // cosine of angle between direction vector and position vector
        cos = -scalar_product(particle_position, particle_direction) / radius;

        if (cos < 1 && cos > 0)
            EXPECT_FALSE(B.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(B.IsInside(particle_position, particle_direction));
    }
}

TEST(DistanceTo, Sphere)
{
    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double radius          = 10;
    double inner_radius    = 0;
    double particle_radius = 0;

    double rnd_phi;
    double rnd_theta;
    double rnd_inner_radius;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int i = 0; i < 11; i++)
    {

        particle_radius = 2. + i * 2.;

        for (int j = 0; j < number_particles; j++)
        {

            rnd_inner_radius = RandomDouble();
            inner_radius     = radius * rnd_inner_radius;

            // The values are divided by 100 to convert the units...
            // Init functions expects m but here everthing is in cm
            Sphere A(Vector3D(0, 0, 0), radius, inner_radius);

            rnd_phi   = RandomDouble();
            rnd_theta = RandomDouble();

            particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
            particle_direction.CalculateCartesianFromSpherical();

            // Chose particle location and angle
            particle_position.SetSphericalCoordinates(particle_radius, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
            particle_position.CalculateCartesianFromSpherical();
            particle_position = -particle_position;

            distance = DistanceToBorder(A, particle_position, particle_direction);

            if (particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius + inner_radius, 1e-8 * (particle_radius + inner_radius));
                ASSERT_NEAR(distance.second, particle_radius + radius, 1e-8 * (particle_radius + radius));
            }
            if (particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius - radius, 1e-8 * (particle_radius - radius));
                ASSERT_NEAR(distance.second, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius == radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }

            if (particle_radius >= radius)
            {
                particle_position = -1 * particle_position;

                // Now the particle is moving away from the sphere so we expect no intersection
                distance = DistanceToBorder(A, particle_position, particle_direction);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);
            }
            if (particle_radius > 20)
            {
                particle_position.SetSphericalCoordinates(inner_radius, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
                particle_position.CalculateCartesianFromSpherical();
                particle_position = -particle_position;

                distance = DistanceToBorder(A, particle_position, particle_direction);
                ASSERT_NEAR(distance.first, 2 * inner_radius, 1e-8 * (2 * inner_radius));
                ASSERT_NEAR(distance.second, inner_radius + radius, 1e-8 * (inner_radius + radius));
            }
        }
    }
}

//
TEST(DistanceTo, Cylinder)
{
    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double height          = 10;
    double radius          = 10;
    double inner_radius    = 0;
    double particle_radius = 0;

    double rnd_phi;
    double rnd_inner_radius;

    double z;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int i = 0; i < 10; i++)
    {

        particle_radius = 2. + i * 2.;

        for (int j = 0; j < number_particles; j++)
        {

            rnd_inner_radius = RandomDouble();
            inner_radius     = radius * rnd_inner_radius;

            // The values are divided by 100 to convert the units...
            // Init functions expects m but here everthing is in cm
            Cylinder A(Vector3D(0, 0, 0), radius, inner_radius, height);

            rnd_phi = RandomDouble();

            // Chose particle location and angle
            particle_position.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, 0.5 * M_PI);
            particle_position.CalculateCartesianFromSpherical();
            particle_position.SetCartesianCoordinates(particle_radius * particle_position.GetX(),
                                                      particle_radius * particle_position.GetY(),
                                                      0.5 * height * particle_position.GetZ());
            particle_position = -particle_position;

            particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, 0.5 * M_PI);
            particle_direction.CalculateCartesianFromSpherical();

            distance = DistanceToBorder(A, particle_position, particle_direction);

            if (particle_radius < radius && particle_radius > inner_radius)
            {
                EXPECT_EQ(distance.second, -1.);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius <= inner_radius)
            {
                ASSERT_NEAR(distance.first, particle_radius + inner_radius, 1e-8 * (particle_radius + inner_radius));
                ASSERT_NEAR(distance.second, particle_radius + radius, 1e-8 * (particle_radius + radius));
            }
            if (particle_radius > radius)
            {
                ASSERT_NEAR(distance.first, particle_radius - radius, 1e-8 * (particle_radius - radius));
                ASSERT_NEAR(distance.second, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }
            if (particle_radius == radius)
            {
                EXPECT_EQ(distance.second, -1);
                ASSERT_NEAR(distance.first, particle_radius - inner_radius, 1e-8 * (particle_radius - inner_radius));
            }

            if (particle_radius >= radius)
            {
                particle_position = -particle_position;

                // Now the particle is moving away from the sphere so we expect no intersection
                distance = DistanceToBorder(A, particle_position, particle_direction);
                EXPECT_EQ(distance.first, -1.);
                EXPECT_EQ(distance.second, -1.);
            }
            if (particle_radius > 20)
            {
                particle_radius = inner_radius;

                particle_position.SetSphericalCoordinates(particle_radius, rnd_phi * 2 * M_PI, 0.5 * M_PI);
                particle_position.CalculateCartesianFromSpherical();
                particle_position = -particle_position;

                distance = DistanceToBorder(A, particle_position, particle_direction);
                ASSERT_NEAR(distance.first, 2 * inner_radius, 1e-8 * (2 * inner_radius));
                ASSERT_NEAR(distance.second, inner_radius + radius, 1e-8 * (inner_radius + radius));
            }
        }
    }

    // One test for inner_radius =0
    inner_radius = 0;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Cylinder B(Vector3D(0, 0, 0), radius, inner_radius, height);

    // Chose particle location and angle

    particle_position.SetCartesianCoordinates(0, 0, height + 10);

    particle_direction.SetSphericalCoordinates(1, 0, M_PI);
    particle_direction.CalculateCartesianFromSpherical();

    z = particle_position.GetZ();

    distance = DistanceToBorder(B, particle_position, particle_direction);

    ASSERT_NEAR(distance.first, z - 0.5 * height, 1e-8 * (z - 0.5 * height));
    ASSERT_NEAR(distance.second, z + 0.5 * height, 1e-8 * (z + 0.5 * height));

    double rnd_alpha;
    double alpha;

    inner_radius = 6;

    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    Cylinder A(Vector3D(0, 0, 0), radius, inner_radius, height);

    for (int j = 0; j < number_particles; j++)
    {

        rnd_alpha = RandomDouble();

        rnd_phi = RandomDouble();

        alpha = 0.3 * M_PI * rnd_alpha;

        // Chose particle location and angle

        particle_position.SetCartesianCoordinates(0, 0, height + 0.5 * height);
        z = particle_position.GetZ();

        particle_direction.SetSphericalCoordinates(1, rnd_phi, M_PI - alpha);
        particle_direction.CalculateCartesianFromSpherical();

        double dist1 = inner_radius / std::sin(alpha);
        double dist2 = radius / std::sin(alpha);

        distance = DistanceToBorder(A, particle_position, particle_direction);

        //  case 1 throught inner cylinder => no intersection
        //  ___  x    ___
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |   | |   |   |
        // |___| |   |___|
        //       |
        if (alpha < std::atan(inner_radius / (z + 0.5 * height)))
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
        //  case 2 first inner cylinder then bottom surface
        //  ___  x      ___
        // |   |  \    |   |
        // |   |   \   |   |
        // |   |    \  |   |
        // |   |     \ |   |
        // |   |      *|   |
        // |___|       |\ _|
        //               *

        else if (alpha < std::atan(radius / (z + 0.5 * height)))
        {
            dist2 = (z + 0.5 * height) / std::cos(alpha);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        //  case 3 first inner cylinder then outer cylinder
        //  ___     x   ___
        // |   |     \ |   |
        // |   |      \|   |
        // |   |       *   |
        // |   |       |\  |
        // |   |       | \ |
        // |   |       |  \|
        // |   |       |   *
        // |___|       |___|\
        //

        else if (alpha < std::atan(inner_radius / height))
        {
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        //  case 4 first upper surface then outer cylinder
        //            x
        //             \
        //  ___         \__
        // |   |       | * |
        // |   |       |  \|
        // |   |       |   *
        // |   |       |   |\
        // |   |       |   |
        // |   |       |   |
        // |___|       |___|
        //
        else if (alpha < std::atan(radius / height))
        {
            dist1 = height / std::cos(alpha);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        //  case 5  no intersection
        //      x_____________
        //  ___      ___
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |   |    |   |
        // |___|    |___|
        //
        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }
}

TEST(DistanceTo, Box)
{
    Vector3D particle_position(0, 0, 0);
    Vector3D particle_direction(0, 0, 0);

    double width  = 10;
    double height = width;

    double rnd_phi;
    double rnd_theta;

    double phi;

    double dist;
    double dist1;
    double dist2;

    std::pair<double, double> distance;

    // MathModel M;
    int number_particles = 1e5;

    std::cout.precision(16);

    for (int j = 0; j < number_particles; j++)
    {
        rnd_phi = RandomDouble();

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        Box A(Vector3D(0, 0, 0), width, width, height);

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, 0.5 * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        distance = DistanceToBorder(A, particle_position, particle_direction);

        phi = particle_direction.GetPhi() * 180. / M_PI;
        if (phi < 45)
        {
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 90)
        {
            phi  = 90 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 135)
        {
            phi  = phi - 90;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 180)
        {
            phi  = 180 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 225)
        {
            phi  = phi - 180;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 270)
        {
            phi  = 270 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 315)
        {
            phi  = phi - 270;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        } else if (phi < 360)
        {
            phi  = 360 - phi;
            dist = 0.5 * width / std::cos(phi / 180 * M_PI);
            EXPECT_EQ(distance.second, -1.);
            ASSERT_NEAR(distance.first, dist, 1e-8 * (dist));
        }
    }

    Box A(Vector3D(0, 0, 0), width, width, height);

    for (int i = 0; i < number_particles; i++)
    {
        rnd_phi = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 0.5 * M_PI, 0.5 * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        particle_position.SetCartesianCoordinates(-1 * width, 0, 0);

        phi      = particle_direction.GetPhi();
        distance = DistanceToBorder(A, particle_position, particle_direction);

        //                       ________________           z|
        //                      |               |            |
        //                      |               |            |_____
        //                      |               |                  x
        //     x----------------*---------------*--------->
        //                      |               |
        //                      |               |
        //                      |               |
        //                      |_______________|

        if (phi < std::atan(0.5 * width / (0.5 * width - particle_position.GetX())))
        {
            dist1 = (-particle_position.GetX() - 0.5 * width) / std::cos(phi);
            dist2 = (0.5 * width - particle_position.GetX()) / std::cos(phi);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        //
        //                          ^
        //                         /                        z|
        //                       _*_____________             |
        //                      |/              |            |_____
        //                      *               |                  x
        //                     /|               |
        //                    / |               |
        //                   /  |               |
        //                  x   |               |
        //                      |               |
        //                      |_______________|

        else if (phi < std::atan(width * 0.5 / (-particle_position.GetX() - 0.5 * width)))
        {
            dist1 = (-particle_position.GetX() - 0.5 * width) / std::cos(phi);
            dist2 = 0.5 * width / std::sin(phi);
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        //                       ^
        //                      /
        //                     / _______________
        //                    / |               |  z|
        //                   /  |               |   |
        //                  /   |               |   |_____
        //                 /    |               |        x
        //                /     |               |
        //               x      |               |
        //                      |               |
        //                      |_______________|

        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }

    // and one test for z surfaces
    for (int i = 0; i < number_particles; i++)
    {
        rnd_theta = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, 0, rnd_theta * 0.5 * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        particle_position.SetCartesianCoordinates(0, 0, -1 * height);

        distance = DistanceToBorder(A, particle_position, particle_direction);

        //                       ________________       x|
        //                      |               |        |
        //                      |               |        |_____
        //                      |               |              z
        //     x----------------*---------------*--------->
        //                      |               |
        //                      |               |
        //                      |               |
        //                      |_______________|

        if (particle_direction.GetTheta() < std::atan(0.5 * height / (0.5 * height - particle_position.GetZ())))
        {
            dist1 = (-particle_position.GetZ() - 0.5 * height) / std::cos(particle_direction.GetTheta());
            dist2 = (0.5 * height - particle_position.GetZ()) / std::cos(particle_direction.GetTheta());
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));

        }
        //
        //                          ^
        //                         /                    x|
        //                       _*_____________         |
        //                      |/              |        |_____
        //                      *               |             z
        //                     /|               |
        //                    / |               |
        //                   /  |               |
        //                  x   |               |
        //                      |               |
        //                      |_______________|

        else if (particle_direction.GetTheta() < std::atan(height * 0.5 / (-particle_position.GetZ() - 0.5 * height)))
        {
            dist1 = (-particle_position.GetZ() - 0.5 * height) / std::cos(particle_direction.GetTheta());
            dist2 = 0.5 * height / std::sin(particle_direction.GetTheta());
            ASSERT_NEAR(distance.first, dist1, 1e-8 * (dist1));
            ASSERT_NEAR(distance.second, dist2, 1e-8 * (dist2));
        }
        //                       ^
        //                      /
        //                     / _______________        x|
        //                    / |               |        |
        //                   /  |               |        |_____
        //                  /   |               |             z
        //                 /    |               |
        //                /     |               |
        //               x      |               |
        //                      |               |
        //                      |_______________|

        else
        {
            EXPECT_EQ(distance.first, -1);
            EXPECT_EQ(distance.second, -1);
        }
    }
}

// =========================================================================
// Fix 1: operator< uses dynamic type, not pointer type
// Before the fix, typeid(this) compared pointer types (always equal),
// so cross-type ordering fell through to name comparison.
// After the fix, typeid(*this) compares the actual dynamic types.
// =========================================================================
TEST(Ordering, CrossTypeOrdering) {
    // Different geometry types must have a consistent total order
    Box box(10, 8, 6);
    Sphere sphere(5, 0);
    Cylinder cyl(4, 0, 10);

    // The key property: for different types, exactly one of a<b or b<a is true
    // (strict weak ordering requires this for non-equal elements of different types)
    bool box_lt_sphere = box < sphere;
    bool sphere_lt_box = sphere < box;
    EXPECT_NE(box_lt_sphere, sphere_lt_box)
        << "Cross-type comparison must be asymmetric (one true, one false)";

    bool box_lt_cyl = box < cyl;
    bool cyl_lt_box = cyl < box;
    EXPECT_NE(box_lt_cyl, cyl_lt_box)
        << "Cross-type comparison must be asymmetric";

    bool sphere_lt_cyl = sphere < cyl;
    bool cyl_lt_sphere = cyl < sphere;
    EXPECT_NE(sphere_lt_cyl, cyl_lt_sphere)
        << "Cross-type comparison must be asymmetric";

    // Irreflexivity: nothing is less than itself.
    EXPECT_FALSE(box < box);
    EXPECT_FALSE(sphere < sphere);
    EXPECT_FALSE(cyl < cyl);

    // Transitivity, verified from the independent pairwise results (not by
    // sorting with the comparator under test, which would be tautological).
    // A strict order over 3 elements must be acyclic: the tournament
    // formed by the three pairwise '<' relations has no 3-cycle.
    bool cyc1 = box_lt_sphere && sphere_lt_cyl && cyl_lt_box;
    bool cyc2 = sphere_lt_box && cyl_lt_sphere && box_lt_cyl;
    EXPECT_FALSE(cyc1) << "Ordering has a 3-cycle (box<sphere<cyl<box)";
    EXPECT_FALSE(cyc2) << "Ordering has a 3-cycle (sphere<box, cyl<sphere, box<cyl)";

    // Spot transitivity: derive the minimum from pairwise results and
    // confirm it precedes the maximum independently.
    Geometry const * lo = box_lt_sphere ? static_cast<Geometry const *>(&box)
                                         : static_cast<Geometry const *>(&sphere);
    Geometry const * hi = box_lt_sphere ? static_cast<Geometry const *>(&sphere)
                                         : static_cast<Geometry const *>(&box);
    bool lo_lt_cyl = *lo < cyl, cyl_lt_lo = cyl < *lo;
    bool hi_lt_cyl = *hi < cyl, cyl_lt_hi = cyl < *hi;
    if(lo_lt_cyl && cyl_lt_hi)
        EXPECT_TRUE(*lo < *hi) << "Transitivity: lo<cyl<hi implies lo<hi";
    if(cyl_lt_lo)
        EXPECT_TRUE(cyl_lt_hi) << "Transitivity: cyl<lo<hi implies cyl<hi";
    (void)hi_lt_cyl;
}

TEST(Ordering, SameTypeDifferentParams) {
    // Two spheres with different radii
    Sphere small_sphere(3, 0);
    Sphere big_sphere(10, 0);

    bool s_lt_b = small_sphere < big_sphere;
    bool b_lt_s = big_sphere < small_sphere;
    // They are not equal, so exactly one must be true
    EXPECT_NE(s_lt_b, b_lt_s);
    // Irreflexivity
    EXPECT_FALSE(small_sphere < small_sphere);
    EXPECT_FALSE(big_sphere < big_sphere);
}

TEST(Ordering, SameTypeEqualParams) {
    Sphere a(5, 2);
    Sphere b(5, 2);
    // Equal objects: neither is less than the other
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(b < a);
}

// =========================================================================
// Fix 10: Degenerate shape validation
// =========================================================================
TEST(Validation, ConeZeroOuterRadii) {
    // Both outer radii zero should throw
    EXPECT_THROW(Cone(0, 0, 0, 0, 10), std::invalid_argument);
    // One positive outer radius is fine
    EXPECT_NO_THROW(Cone(0, 5, 0, 0, 8));
    EXPECT_NO_THROW(Cone(0, 0, 0, 5, 8));
}

TEST(Validation, TrdZeroHeight) {
    EXPECT_THROW(Trd(5, 3, 4, 2, 0), std::invalid_argument);
    EXPECT_THROW(Trd(5, 3, 4, 2, -1), std::invalid_argument);
    // Positive height is fine
    EXPECT_NO_THROW(Trd(5, 3, 4, 2, 6));
}

// --- CutTube negative rmin validation (#11) ---

TEST(Validation, CutTubeNegativeRmin) {
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    // Negative rmin should throw (even if rmax is positive and swap would "fix" it)
    EXPECT_THROW(CutTube(-5.0, 10.0, 5.0, low_norm, high_norm), std::invalid_argument);
    // Negative rmax should also throw
    EXPECT_THROW(CutTube(0.0, -3.0, 5.0, low_norm, high_norm), std::invalid_argument);
    // Both negative
    EXPECT_THROW(CutTube(-2.0, -1.0, 5.0, low_norm, high_norm), std::invalid_argument);
    // Valid: rmin > rmax triggers swap, but both are non-negative
    EXPECT_NO_THROW(CutTube(8.0, 5.0, 3.0, low_norm, high_norm));
    // Valid: zero rmin is fine
    EXPECT_NO_THROW(CutTube(0.0, 5.0, 3.0, low_norm, high_norm));
}

TEST(Validation, CutTubeNegativeRminWithPlacement) {
    Placement p;
    Vector3D low_norm(0, 0, -1);
    Vector3D high_norm(0, 0, 1);
    EXPECT_THROW(CutTube(p, -5.0, 10.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_THROW(CutTube(p, 0.0, -3.0, 5.0, low_norm, high_norm), std::invalid_argument);
    EXPECT_NO_THROW(CutTube(p, 8.0, 5.0, 3.0, low_norm, high_norm));
}

// --- Trap infinity sentinel test (#12) ---

TEST(TrapIntersection, FarFieldRay) {
    // A ray originating far from the trap must still produce intersections.
    // At 1e6 distance the double ULP is ~1e6*2^-52 ~ 2e-10, so a 1e-6
    // positional tolerance is amply safe and still exercises the
    // far-field path (NearFieldBasic uses distance 100). The previous
    // 1e12 distance forced an unnecessarily loose 1e-3 tolerance.
    Trap trap(10, 0, 0, 8, 5, 5, 0, 6, 3, 3, 0);
    Vector3D origin(0, 0, -1e6);
    Vector3D dir(0, 0, 1);
    auto hits = trap.Intersections(origin, dir);
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_NEAR(hits[0].position.GetZ(), -10.0, 1e-6);
    EXPECT_NEAR(hits[1].position.GetZ(),  10.0, 1e-6);
}

TEST(TrapIntersection, NearFieldBasic) {
    // Basic sanity check that Trap intersections work for normal distances
    Trap trap(10, 0, 0, 8, 5, 5, 0, 6, 3, 3, 0);
    Vector3D origin(0, 0, -100);
    Vector3D dir(0, 0, 1);
    auto hits = trap.Intersections(origin, dir);
    ASSERT_EQ(hits.size(), 2u);
    EXPECT_NEAR(hits[0].position.GetZ(), -10.0, 1e-6);
    EXPECT_NEAR(hits[1].position.GetZ(),  10.0, 1e-6);
}

// --- Ellipsoid z-cut boundary dedup test (#18) ---

TEST(EllipsoidIntersection, ZCutBoundaryNoDuplicates) {
    // Ellipsoid with z-cuts at the equator: zcut1=-5, zcut2=5
    // with semi-axes (10, 10, 10) -- effectively a sphere with z-cuts.
    Ellipsoid ell(10.0, 10.0, 10.0, -5.0, 5.0);
    // Shoot a ray along z through the point where the ellipsoid surface
    // meets the z-cut boundary at z=5 (r = sqrt(100-25) ~ 8.66 from z-axis)
    double r_at_cut = std::sqrt(10.0*10.0 - 5.0*5.0);
    Vector3D origin(r_at_cut, 0, -20);
    Vector3D dir(0, 0, 1);
    auto hits = ell.Intersections(origin, dir);
    // Must have an even number of intersections and no duplicates
    EXPECT_EQ(hits.size() % 2, 0u);
    for(size_t i = 1; i < hits.size(); ++i) {
        EXPECT_GT(hits[i].distance - hits[i-1].distance, Geometry::GEOMETRY_PRECISION);
    }
}

TEST(EllipsoidIntersection, ZCutBoundaryEvenCount) {
    // Test with a ray that hits exactly at zcut2 boundary
    Ellipsoid ell(5.0, 5.0, 10.0, -8.0, 8.0);
    // At z=8, ellipsoid radius in xy = 5*sqrt(1 - 64/100) = 5*0.6 = 3
    double r_at_cut = 5.0 * std::sqrt(1.0 - (8.0*8.0)/(10.0*10.0));
    Vector3D origin(r_at_cut - 0.001, 0, -20);
    Vector3D dir(0, 0, 1);
    auto hits = ell.Intersections(origin, dir);
    // Deterministic: at x just inside the cut-plane rim radius the
    // ellipsoid surface (z=+/-8.001) lies outside the z-cuts (+/-8), so
    // the solid is capped by the two flat z-cut discs. The axial ray
    // enters the bottom disc and exits the top disc: exactly 2.
    EXPECT_EQ(hits.size(), 2u);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


