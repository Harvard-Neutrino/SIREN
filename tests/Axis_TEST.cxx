
#include <math.h>
#include <cmath>
#include <random>
#include <iostream>
#include <gtest/gtest.h>

#include "earthmodel-service/DensityDist.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Polynomial.h"

using namespace earthmodel;

std::mt19937 rng_;
std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);

double RandomDouble() {
    return uniform_distribution(rng_);
}

Vector3D RandomVector() {
    return Vector3D(RandomDouble()*4-2, RandomDouble()*4-2, RandomDouble()*4-2);
}

Vector3D RandomDirection() {
    double az = RandomDouble()*2*M_PI;
    double zen = RandomDouble()*M_PI;
    Vector3D res;
    res.SetSphericalCoordinates(1.0, az, zen);
    res.CalculateCartesianFromSpherical();
    return res;
}

TEST(Comparison, Comparison_equal)
{
    Axis1D* r_ax_A = new RadialAxis1D();
    Axis1D* r_ax_B = new RadialAxis1D();
    RadialAxis1D r_ax_C;
    RadialAxis1D r_ax_D;
    EXPECT_TRUE(*r_ax_A == *r_ax_B);
    EXPECT_TRUE(*r_ax_A == r_ax_C);
    EXPECT_TRUE(r_ax_C == r_ax_D);

    delete r_ax_A;
    delete r_ax_B;

    r_ax_A = new RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_B = new RadialAxis1D(Vector3D(1,0,0), Vector3D(1,0,0));
    r_ax_C = RadialAxis1D(Vector3D(1,0,0),Vector3D(0,0,0));
    r_ax_D = RadialAxis1D(Vector3D(1,0,0),Vector3D(1,0,0));
    EXPECT_FALSE(*r_ax_A == *r_ax_B);
    EXPECT_TRUE(*r_ax_A == r_ax_C);
    EXPECT_FALSE(*r_ax_A == r_ax_D);
    EXPECT_FALSE(r_ax_C == r_ax_D);

    delete r_ax_A;
    delete r_ax_B;

    r_ax_A = new RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_B = new RadialAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    r_ax_C = RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_D = RadialAxis1D(Vector3D(1,1,0), Vector3D(0,0,0));
    EXPECT_TRUE(*r_ax_A == *r_ax_B);
    EXPECT_TRUE(*r_ax_A == r_ax_C);
    EXPECT_TRUE(*r_ax_A == r_ax_D);
    EXPECT_TRUE(r_ax_C == r_ax_D);

    delete r_ax_A;
    delete r_ax_B;
    r_ax_A = new RadialAxis1D();
    r_ax_B = new RadialAxis1D();
    r_ax_C = RadialAxis1D();
    r_ax_D = RadialAxis1D();

    Axis1D* c_ax_A = new CartesianAxis1D();
    Axis1D* c_ax_B = new CartesianAxis1D();
    CartesianAxis1D c_ax_C;
    CartesianAxis1D c_ax_D;
    EXPECT_TRUE(*c_ax_A == *c_ax_B);
    EXPECT_TRUE(c_ax_C == c_ax_D);
    EXPECT_TRUE(*c_ax_A == c_ax_C);

    EXPECT_FALSE(*r_ax_A == *c_ax_A);
    EXPECT_FALSE(*r_ax_B == *c_ax_B);
    EXPECT_FALSE(r_ax_C == c_ax_C);
    EXPECT_FALSE(r_ax_D == c_ax_D);

    delete c_ax_A;
    delete c_ax_B;

    c_ax_A = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_B = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(1,0,0));
    c_ax_C = CartesianAxis1D(Vector3D(1,0,0),Vector3D(0,0,0));
    c_ax_D = CartesianAxis1D(Vector3D(1,0,0),Vector3D(1,0,0));
    EXPECT_FALSE(*c_ax_A == *c_ax_B);
    EXPECT_TRUE(*c_ax_A == c_ax_C);
    EXPECT_FALSE(*c_ax_A == c_ax_D);
    EXPECT_FALSE(c_ax_C == c_ax_D);

    EXPECT_FALSE(*r_ax_A == *c_ax_A);
    EXPECT_FALSE(*r_ax_B == *c_ax_B);
    EXPECT_FALSE(r_ax_C == c_ax_C);
    EXPECT_FALSE(r_ax_D == c_ax_D);

    delete c_ax_A;
    delete c_ax_B;

    c_ax_A = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_B = new CartesianAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    c_ax_C = CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_D = CartesianAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    EXPECT_FALSE(*c_ax_A == *c_ax_B);
    EXPECT_TRUE(*c_ax_A == c_ax_C);
    EXPECT_FALSE(*c_ax_A == c_ax_D);
    EXPECT_FALSE(c_ax_C == c_ax_D);

    EXPECT_FALSE(*r_ax_A == *c_ax_A);
    EXPECT_FALSE(*r_ax_B == *c_ax_B);
    EXPECT_FALSE(r_ax_C == c_ax_C);
    EXPECT_FALSE(r_ax_D == c_ax_D);

    delete r_ax_A;
    delete r_ax_B;
    delete c_ax_A;
    delete c_ax_B;
}

TEST(Comparison, Comparison_not_equal)
{
    Axis1D* r_ax_A = new RadialAxis1D();
    Axis1D* r_ax_B = new RadialAxis1D();
    RadialAxis1D r_ax_C;
    RadialAxis1D r_ax_D;
    EXPECT_FALSE(*r_ax_A != *r_ax_B);
    EXPECT_FALSE(*r_ax_A != r_ax_C);
    EXPECT_FALSE(r_ax_C != r_ax_D);

    delete r_ax_A;
    delete r_ax_B;

    r_ax_A = new RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_B = new RadialAxis1D(Vector3D(1,0,0), Vector3D(1,0,0));
    r_ax_C = RadialAxis1D(Vector3D(1,0,0),Vector3D(0,0,0));
    r_ax_D = RadialAxis1D(Vector3D(1,0,0),Vector3D(1,0,0));
    EXPECT_TRUE(*r_ax_A != *r_ax_B);
    EXPECT_FALSE(*r_ax_A != r_ax_C);
    EXPECT_TRUE(*r_ax_A != r_ax_D);
    EXPECT_TRUE(r_ax_C != r_ax_D);

    delete r_ax_A;
    delete r_ax_B;

    r_ax_A = new RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_B = new RadialAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    r_ax_C = RadialAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    r_ax_D = RadialAxis1D(Vector3D(1,1,0), Vector3D(0,0,0));
    EXPECT_FALSE(*r_ax_A != *r_ax_B);
    EXPECT_FALSE(*r_ax_A != r_ax_C);
    EXPECT_FALSE(*r_ax_A != r_ax_D);
    EXPECT_FALSE(r_ax_C != r_ax_D);

    delete r_ax_A;
    delete r_ax_B;
    r_ax_A = new RadialAxis1D();
    r_ax_B = new RadialAxis1D();
    r_ax_C = RadialAxis1D();
    r_ax_D = RadialAxis1D();

    Axis1D* c_ax_A = new CartesianAxis1D();
    Axis1D* c_ax_B = new CartesianAxis1D();
    CartesianAxis1D c_ax_C;
    CartesianAxis1D c_ax_D;
    EXPECT_FALSE(*c_ax_A != *c_ax_B);
    EXPECT_FALSE(c_ax_C != c_ax_D);
    EXPECT_FALSE(*c_ax_A != c_ax_C);

    EXPECT_TRUE(*r_ax_A != *c_ax_A);
    EXPECT_TRUE(*r_ax_B != *c_ax_B);
    EXPECT_TRUE(r_ax_C != c_ax_C);
    EXPECT_TRUE(r_ax_D != c_ax_D);

    delete c_ax_A;
    delete c_ax_B;

    c_ax_A = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_B = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(1,0,0));
    c_ax_C = CartesianAxis1D(Vector3D(1,0,0),Vector3D(0,0,0));
    c_ax_D = CartesianAxis1D(Vector3D(1,0,0),Vector3D(1,0,0));
    EXPECT_TRUE(*c_ax_A != *c_ax_B);
    EXPECT_FALSE(*c_ax_A != c_ax_C);
    EXPECT_TRUE(*c_ax_A != c_ax_D);
    EXPECT_TRUE(c_ax_C != c_ax_D);

    EXPECT_TRUE(*r_ax_A != *c_ax_A);
    EXPECT_TRUE(*r_ax_B != *c_ax_B);
    EXPECT_TRUE(r_ax_C != c_ax_C);
    EXPECT_TRUE(r_ax_D != c_ax_D);

    delete c_ax_A;
    delete c_ax_B;

    c_ax_A = new CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_B = new CartesianAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    c_ax_C = CartesianAxis1D(Vector3D(1,0,0), Vector3D(0,0,0));
    c_ax_D = CartesianAxis1D(Vector3D(0,1,0), Vector3D(0,0,0));
    EXPECT_TRUE(*c_ax_A != *c_ax_B);
    EXPECT_FALSE(*c_ax_A != c_ax_C);
    EXPECT_TRUE(*c_ax_A != c_ax_D);
    EXPECT_TRUE(c_ax_C != c_ax_D);

    EXPECT_TRUE(*r_ax_A != *c_ax_A);
    EXPECT_TRUE(*r_ax_B != *c_ax_B);
    EXPECT_TRUE(r_ax_C != c_ax_C);
    EXPECT_TRUE(r_ax_D != c_ax_D);

    delete r_ax_A;
    delete r_ax_B;
    delete c_ax_A;
    delete c_ax_B;
}

TEST(Assignment, Copyconstructor)
{
    CartesianAxis1D c_ax_A(Vector3D(0,0,1),Vector3D(1,1,1));
    CartesianAxis1D c_ax_B = c_ax_A;
    CartesianAxis1D c_ax_C(c_ax_A);
    EXPECT_TRUE(c_ax_A == c_ax_B);
    EXPECT_TRUE(c_ax_A == c_ax_C);

    RadialAxis1D r_ax_A(Vector3D(0,0,1),Vector3D(1,1,1));
    RadialAxis1D r_ax_B = r_ax_A;
    RadialAxis1D r_ax_C(r_ax_A);
    EXPECT_TRUE(r_ax_A == r_ax_B);
    EXPECT_TRUE(r_ax_A == r_ax_C);
}

TEST(Copy, clone)
{
    CartesianAxis1D c_ax_A;
    RadialAxis1D r_ax_A;

    Axis1D* c_ax_p = c_ax_A.clone();
    EXPECT_TRUE(c_ax_A == *c_ax_p);

    Axis1D* r_ax_p = r_ax_A.clone();
    EXPECT_TRUE(r_ax_A == *r_ax_p);

    EXPECT_TRUE((c_ax_A == r_ax_A) == (*c_ax_p == *r_ax_p));

    delete c_ax_p;
    delete r_ax_p;
}

TEST(Copy, create)
{
    CartesianAxis1D A;
    RadialAxis1D B;

    std::shared_ptr<const Axis1D> Ap = A.create();
    EXPECT_TRUE(A == *Ap);

    std::shared_ptr<const Axis1D> Bp = B.create();
    EXPECT_TRUE(B == *Bp);

    EXPECT_TRUE((A == B) == (*Ap == *Bp));

}

TEST(Evaluation, Cartesian)
{
    unsigned int N_RAND = 100;
    for(unsigned int i=0; i<N_RAND; ++i) {
        Vector3D center = RandomVector();
        Vector3D direction = RandomDirection();
        CartesianAxis1D cax(direction, center);
        Axis1D* ax = &cax;

        // Center is zero
        EXPECT_DOUBLE_EQ(ax->GetX(center), 0.0);
        for(unsigned int j=0; j<N_RAND; ++j) {
            // works along line
            double distance = RandomDouble()*20-10;
            EXPECT_NEAR(ax->GetX(center + distance*direction), distance, std::abs(distance)*1e-8);

            // X is constant wrt normal plane
            Vector3D point = center + distance*direction;
            double x_on_point = ax->GetX(point);
            for(unsigned int k=0; k<N_RAND; ++k) {
                double r = RandomDouble()*20-10;
                Vector3D new_dir = RandomDirection();
                new_dir = new_dir - (direction*new_dir)*direction;
                new_dir.normalize();
                new_dir = new_dir - (direction*new_dir)*direction;
                new_dir.normalize();
                new_dir = new_dir - (direction*new_dir)*direction;
                EXPECT_NEAR(direction*new_dir, 0.0, 1e-15);
                EXPECT_NEAR(ax->GetX(point+r*new_dir), x_on_point, std::abs(x_on_point)*1e-8);
            }
        }

        // dX is constant wrt different points
        for(unsigned int j=0; j<N_RAND; ++j) {
            Vector3D point = RandomVector();
            Vector3D new_dir = RandomDirection();
            double dX = new_dir*direction;
            EXPECT_DOUBLE_EQ(ax->GetdX(point, new_dir), dX);
            for(unsigned int k=0; k<N_RAND; ++k) {
                Vector3D new_point = RandomVector();
                EXPECT_DOUBLE_EQ(ax->GetdX(new_point, new_dir), dX);
            }
        }
    }
}

TEST(Evaluation, Radial)
{
    unsigned int N_RAND = 100;
    for(unsigned int i=0; i<N_RAND; ++i) {
        Vector3D center = RandomVector();
        RadialAxis1D rax(center);
        Axis1D* ax = &rax;

        // Center is zero
        EXPECT_DOUBLE_EQ(ax->GetX(center), 0.0);
        for(unsigned int j=0; j<N_RAND; ++j) {
            // works along line
            Vector3D direction = RandomDirection();
            for(unsigned int k=0; k<N_RAND; ++k) {
                double distance = RandomDouble()*20-10;
                EXPECT_NEAR(ax->GetX(center + distance*direction), std::abs(distance), std::abs(distance)*1e-8);
            }
        }

        // Constant wrt sphere surface
        for(unsigned int j=0; j<N_RAND; ++j) {
            double distance = RandomDouble()*20-10;
            for(unsigned int k=0; k<N_RAND; ++k) {
                Vector3D direction = RandomDirection();
                EXPECT_NEAR(ax->GetX(center + distance*direction), std::abs(distance), std::abs(distance)*1e-8);
            }
        }

        for(unsigned int j=0; j<N_RAND; ++j) {
            Vector3D point = RandomVector();
            Vector3D direction = RandomDirection();
            Vector3D r = (point-center);
            double R = r.magnitude();
            r.normalize();
            double dX = r*direction;

            // Check dX
            EXPECT_DOUBLE_EQ(ax->GetdX(point, direction), dX);

            // dX inverts
            EXPECT_NEAR(ax->GetdX(point - 2*R*r, direction), -dX, std::abs(dX)*1e-8);

            // dX constant along line
            for(unsigned int k=0; k<N_RAND; ++k) {
                double delta = (RandomDouble()*4-3)*R;
                Vector3D new_point = point+r*delta;
                if(delta > -R)
                    EXPECT_NEAR(ax->GetdX(new_point, direction), dX, std::abs(dX)*1e-8);
                else
                    EXPECT_NEAR(ax->GetdX(new_point, direction), -dX, std::abs(dX)*1e-8);
            }
        }
    }

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

