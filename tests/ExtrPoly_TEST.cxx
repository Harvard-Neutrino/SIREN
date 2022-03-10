#include <iostream>
#include <stdlib.h>
#include <random>
#include <math.h>

#include <gtest/gtest.h>

#include "earthmodel-service/Geometry.h"

using namespace earthmodel;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

TEST(Comparison, Comparison_equal)
{
    ExtrPoly A;
    ExtrPoly B;
    EXPECT_TRUE(A == B);

    ExtrPoly* C = new ExtrPoly();
    ExtrPoly* D = new ExtrPoly();
    EXPECT_TRUE(*C == *D);

    Geometry* E = new ExtrPoly();
    Geometry* F = new ExtrPoly();
    EXPECT_TRUE(*E == *F);

    delete C;
    delete D;
    delete E;
    delete F;
}

TEST(Comparison, Comparison_not_equal)
{
    ExtrPoly A;
    Box B;
    EXPECT_TRUE(A != B);

    ExtrPoly* C = new ExtrPoly();
    Box* D    = new Box();
    EXPECT_TRUE(*C != *D);

    Geometry* E = new ExtrPoly();
    Geometry* F = new Box();
    EXPECT_TRUE(*E != *F);

    delete C;
    delete D;
    delete E;
    delete F;
}

TEST(Assignment, Copyconstructor)
{
    ExtrPoly A;
    ExtrPoly B = A;

    EXPECT_TRUE(A == B);

    Geometry* C = new ExtrPoly();
    Geometry* D = new Box();

    *D = *C;

    EXPECT_FALSE(*C == *D);

    std::vector<std::vector<double>> poly;
    std::vector<double> polyVert;
    std::vector<ExtrPoly::ZSection> zsecs;
    for (double i = 0; i < 5.0; i+=1) {
			polyVert.push_back(i);
			polyVert.push_back(i*i);
			poly.push_back(polyVert);
			polyVert.clear();
    }
    Geometry* E = new ExtrPoly(poly,zsecs);

    *C = *E;

    EXPECT_TRUE(*C == *E);
}

TEST(Assignment, Copyconstructor2)
{
    ExtrPoly A;
    ExtrPoly B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Operator)
{
    ExtrPoly A;
    std::vector<std::vector<double>> poly;
    std::vector<double> polyVert;
    std::vector<ExtrPoly::ZSection> zsecs;
    for (double i = 0; i < 5.0; i+=1) {
			polyVert.push_back(i);
			polyVert.push_back(i*i);
			poly.push_back(polyVert);
			polyVert.clear();
    }
    double offset[2] = {0,0};
    zsecs.push_back(ExtrPoly::ZSection(0,offset,1));
    ExtrPoly B(poly,zsecs);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    ExtrPoly A;
    ExtrPoly B;
    EXPECT_TRUE(A == B);
    std::vector<std::vector<double>> poly;
    std::vector<double> polyVert;
    std::vector<ExtrPoly::ZSection> zsecs;
    for (double i = 0; i < 5.0; i+=1) {
			polyVert.push_back(i);
			polyVert.push_back(i*i);
			poly.push_back(polyVert);
			polyVert.clear();
    }
    double offset[2] = {0,0};
    zsecs.push_back(ExtrPoly::ZSection(0,offset,1));
    Geometry* C = new ExtrPoly(poly,zsecs);
    Geometry* D = new ExtrPoly(poly,zsecs);

    EXPECT_TRUE(*C == *D);

    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

TEST(DistanceToClosestApproach, Method)
{
    ExtrPoly A;
    Vector3D position = Vector3D(1, -1, 0);
    Vector3D direction = Vector3D(0, 1, 0);
    direction.CalculateSphericalCoordinates();

    double distance_closest_approach = A.DistanceToClosestApproach(position, direction);

    ASSERT_NEAR(distance_closest_approach, 1., 1e-9);
}

TEST(IsInside, ExtrPoly)
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

    double big_width_x = 100;
    double big_width_y = 100;
    double big_height  = 100;

    double extr_face = 5.0;
    double extr_extent = 5.0;

		// Choose extr poly dimensions. Here we make a triangular prism
    std::vector<std::vector<double>> poly;
    std::vector<ExtrPoly::ZSection> zsecs;
		poly.push_back(std::vector<double>{extr_face,-extr_face});
		poly.push_back(std::vector<double>{-extr_face,-extr_face});
		poly.push_back(std::vector<double>{0.0,extr_face});
    double offset[2] = {0,0};
    zsecs.push_back(ExtrPoly::ZSection(-extr_extent,offset,1));
    zsecs.push_back(ExtrPoly::ZSection(extr_extent,offset,1));


    Vector3D position_geometry(0, 0, 0);
    
    Vector3D Del(0, 0, 0);
    double s = 0;
    double t = 0;

    int is_inside  = 0;
    int is_outside = 0;

    double volumia_ratio = 0;

    // MathModel M;
    int number_particles = 1e1;
    int number_volumina  = 1e1;

    Placement p; 

    for (int i = 0; i < number_volumina; i++)
    {

        // Chose the origin of the extr_poly-geometry
        // This poly should be inside the big box in which the particle
        // will be located
        rnd_x0 = RandomDouble();
        rnd_y0 = RandomDouble();
        rnd_z0 = RandomDouble();

        position_geometry.SetCartesianCoordinates((2 * rnd_x0 - 1) * 0.5 * (big_width_x - extr_face),
                                                  (2 * rnd_y0 - 1) * 0.5 * (big_width_y - extr_face),
                                                  (2 * rnd_z0 - 1) * 0.5 * (big_height - extr_extent));


        p.SetPosition(position_geometry);

        // The values are divided by 100 to convert the units...
        // Init functions expects m but here everthing is in cm
        ExtrPoly A(p, poly,zsecs);
        volumia_ratio = 4 * extr_face * extr_face * extr_extent / (big_width_x * big_width_y * big_height);
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

            // if this constraints are true the particle is inside the triangle geometry
            // see https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
            // NOTE: maybe try to generalize this to arbitrary convex poly
            Del = particle_position - position_geometry;
            s = Del.GetX() / (2 * extr_face);
            t = (Del.GetY() - 0.5 * Del.GetX()) / (2 * extr_face); 

            if (0 <= s && s <= 1 &&
                0 <= t && t <= 1 &&
                s+t <=1 &&
                std::abs(Del.GetZ()) < extr_extent)
            {
                is_inside++;
                EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
            } 
            else
            {
                is_outside++;
                if(A.IsInside(particle_position, particle_direction))
                {
									std::cout << "\n\ns t DelZ: " << s << " " << t << " " << Del.GetZ() << "\n\n";
									std::cout << "poly position:\n" << position_geometry << "\n\n";
									std::cout << "particle position:\n" << particle_position << "\n\n";
									std::cout << "particle direction:\n" << particle_direction << "\n\n";
								}
                EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
            }
        }
        ASSERT_NEAR(1. * is_inside, volumia_ratio * number_particles, 3 * sqrt(volumia_ratio * number_particles));
        is_inside  = 0;
        is_outside = 0;
    }
    // Check what happens if particles are on the border of the poly
    // NICK: Implement next
    
    // The values are divided by 100 to convert the units...
    // Init functions expects m but here everthing is in cm
    ExtrPoly A(Vector3D(0, 0, 0), poly,zsecs);

    // Particle is on the top surface.
    // Theta 0° - 90° means particle is moving outside
    // This should be treated as outside
    // Theta 90° - 180° means particle is moving inside (should be treated as inside)
    // The value of phi does not matter
    particle_position.SetCartesianCoordinates(0, 0, extr_extent);
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
    particle_position.SetCartesianCoordinates(0, 0, -extr_extent);
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

    // Surface 1: along x axis
    particle_position.SetCartesianCoordinates(0, -extr_face, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        // phi = 0 is in positive x direction
        if (particle_direction.GetPhi() > M_PI)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface 2: right face
    particle_position.SetCartesianCoordinates(extr_face/2, 0, 0);
    double cosa;
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();
				cosa = (1./sqrt(5.)) * (2*particle_direction.GetX() + particle_direction.GetY());

        // phi = 0 is in positive x direction
        if (cosa>0)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
    }
    // Surface 3: left face
    particle_position.SetCartesianCoordinates(-extr_face/2, 0, 0);
    for (int i = 0; i < 1e4; i++)
    {
        rnd_theta = RandomDouble();
        rnd_phi   = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, rnd_theta * M_PI);
        particle_direction.CalculateCartesianFromSpherical();
				cosa = (1./sqrt(5.)) * (-2*particle_direction.GetX() + particle_direction.GetY());

        // phi = 0 is in positive x direction
        if (cosa>0)
            EXPECT_FALSE(A.IsInside(particle_position, particle_direction));
        else
            EXPECT_TRUE(A.IsInside(particle_position, particle_direction));
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
		
		// Choose extr poly dimensions. Here we make a box
    std::vector<std::vector<double>> poly;
    std::vector<ExtrPoly::ZSection> zsecs;
		poly.push_back(std::vector<double>{width/2,width/2});
		poly.push_back(std::vector<double>{width/2,-width/2});
		poly.push_back(std::vector<double>{-width/2,-width/2});
		poly.push_back(std::vector<double>{-width/2,width/2});
    double offset[2] = {0,0};
    zsecs.push_back(ExtrPoly::ZSection(-height/2,offset,1));
    zsecs.push_back(ExtrPoly::ZSection(height/2,offset,1));
        
		// The values are divided by 100 to convert the units...
		// Init functions expects m but here everthing is in cm
    ExtrPoly A(Vector3D(0, 0, 0), poly,zsecs);

    for (int j = 0; j < number_particles; j++)
    {
        rnd_phi = RandomDouble();


        particle_direction.SetSphericalCoordinates(1, rnd_phi * 2 * M_PI, 0.5 * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        distance = A.DistanceToBorder(particle_position, particle_direction);

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

    //ExtrPoly A(Vector3D(0, 0, 0), poly,zsecs);

    for (int i = 0; i < number_particles; i++)
    {
        rnd_phi = RandomDouble();

        particle_direction.SetSphericalCoordinates(1, rnd_phi * 0.5 * M_PI, 0.5 * M_PI);
        particle_direction.CalculateCartesianFromSpherical();

        particle_position.SetCartesianCoordinates(-1 * width, 0, 0);

        phi      = particle_direction.GetPhi();
        distance = A.DistanceToBorder(particle_position, particle_direction);

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

        distance = A.DistanceToBorder(particle_position, particle_direction);

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

