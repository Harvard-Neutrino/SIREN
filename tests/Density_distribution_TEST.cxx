
#include <gtest/gtest.h>

#include <earthmodel-service/DensityDist.h>
#include <earthmodel-service/Vector3D.h>
#include <earthmodel-service/Polynomial.h>

using namespace earthmodel;

TEST(Comparison, Comparison_equal)
{
    Axis1D* ax_A = new RadialAxis1D();
    Axis1D* ax_B = new RadialAxis1D();
    EXPECT_TRUE(*ax_A == *ax_B);

    Axis1D* ax_C = new CartesianAxis1D();
    CartesianAxis1D ax_D;
    EXPECT_TRUE(*ax_C == ax_D);

    CartesianAxis1D ax_E;
    EXPECT_TRUE(ax_E == ax_D);

    DensityDistribution* A = new DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(0.3);
    DensityDistribution* B = new DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(0.3);
    EXPECT_TRUE(*A == *B);

    double sigma = 1.;
    DensityDistribution* C = new DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_D, sigma);
    DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D> D(ax_D, sigma);
    EXPECT_TRUE(D == *C);

    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> G;
    DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D> H;
    EXPECT_TRUE(G == H);

    delete A;
    delete B;
    delete C;
    delete ax_A;
    delete ax_B;
    delete ax_C;
}

TEST(Comparison, Comparison_not_equal)
{
    Vector3D faxis_default(1,0,0);
    Vector3D faxis_new(0,1,0);
    Vector3D fp0_default(0,0,0);
    Vector3D fp0_new(1,1,1);
    CartesianAxis1D ax_A;
    CartesianAxis1D ax_B(faxis_default, fp0_new);
    CartesianAxis1D ax_C(faxis_new, fp0_default);
    EXPECT_TRUE(ax_A != ax_B);
    EXPECT_TRUE(ax_A != ax_C);

    Axis1D* ax_D = new CartesianAxis1D();
    Axis1D* ax_E = new RadialAxis1D();
    EXPECT_TRUE(*ax_D != *ax_E);

    DensityDistribution* A = new DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>();
    DensityDistribution* B = new DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, 1);
    DensityDistribution* D = new DensityDistribution1D<CartesianAxis1D,ExponentialDistribution1D>(ax_A, 2);
    EXPECT_TRUE(*A != *B);
    EXPECT_TRUE(*B != *D);

    std::vector<double> vecA = {1,2};
    std::vector<double> vecB = {2,3};
    Polynom poly_A(vecA);
    Polynom poly_B(vecB);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> E(ax_A, poly_A);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> F(ax_B, poly_A);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> G(ax_A, poly_B);
    EXPECT_TRUE(E != F);
    EXPECT_TRUE(E != G);

    delete ax_D;
    delete ax_E;
    delete A;
    delete B;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    CartesianAxis1D ax_A;
    CartesianAxis1D ax_B = ax_A;
    CartesianAxis1D ax_C(ax_A);
    EXPECT_TRUE(ax_A == ax_B);
    EXPECT_TRUE(ax_A == ax_C);

    std::vector<double> vecA = {1,2};
    Polynom poly_A(vecA);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> A(ax_A, poly_A);
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> B = A;
    DensityDistribution1D<CartesianAxis1D,PolynomialDistribution1D> C(A);
    EXPECT_TRUE(A == B);
    EXPECT_TRUE(A == C);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
