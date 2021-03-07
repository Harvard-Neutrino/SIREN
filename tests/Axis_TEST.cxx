
#include <gtest/gtest.h>

#include <earthmodel-service/DensityDist.h>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/Polynomial.h>

using namespace earthmodel;

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


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
