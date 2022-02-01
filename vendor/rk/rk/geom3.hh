//=========================================================================
// geom3.hh
//
// A collection of light-weight classes and utilities for calculations
// in 3d geometry. These classes work better than typical "schoolbook"
// implementations. In particular, points and vectors are distinct classes,
// rotations are implemented using quaternions, vector length remains
// constant under rotations, various numerical evaluations are performed
// without unnecessary loss of precision, short functions are inlined,
// everything lives on the stack (no dynamic memory allocation), the code
// is thread-safe, return value optimization is supported, etc.
//
// I. Volobouev
// December 2007
//=========================================================================

#ifndef GEOM3_GEOM3_HH_
#define GEOM3_GEOM3_HH_

#include <cmath>
#include <cassert>
#include <iostream>

// Forward declarations for the "rk" namespace
namespace rk {
    class P4;
    class Point4;
    class LT;
}

namespace geom3 {
    enum {
        X = 0,
        Y,
        Z
    };

    // Forward declarations
    class Point3;
    class UnitVector3;
    class Vector3;
    class Rotation3;
    class Matrix3x3;

    // None of the classes below have virtual destructors.
    // This is intentional -- keep these simple classes fast
    // and efficient, and preferably use them by composition.

    // The "UnitVector3" class is used to represent directions
    // in the 3d space. Note that polar angle values are limited
    // in range from 0 to pi, and the code will abort if the polar
    // angle argument is out of this range. Strict adherence to this
    // convention avoids some subtle bugs (e.g., when two directions
    // are compared).
    class UnitVector3
    {
    public:
        UnitVector3(double theta, double phi);
        UnitVector3(double x, double y, double z, bool mustNormalize=true);

        // Fast constructor for directions in the XY plane.
        // Z component will be 0.
        explicit UnitVector3(double phi);

        // There is no special constructor for making a unit vector
        // out of eta and phi. However, this can be done very easily
        // using the component-wise constructor:
        // UnitVector3(cos(phi), sin(phi), sinh(eta));

        // Modifiers of the direction angles
        UnitVector3& setTheta(double value);
        UnitVector3& setPhi(double value);
        UnitVector3& setEta(double value);

        // Basic accessors
        double x() const;
        double y() const;
        double z() const;

        // Polar angle
        double theta() const;
        double cosTheta() const;

        // Azimuthal angle
        double phi() const;

        // Pseudorapidity
        double eta() const;

        // Scalar product
        double dot(const UnitVector3& r) const;
        double dot(const Vector3& r) const;

        // Cross product
        Vector3 cross(const UnitVector3& r) const;
        Vector3 cross(const Vector3& r) const;

        // Angle between two directions
        double angle(const UnitVector3& r) const;
        double angle(const Vector3& r) const;

        // Unary operators
        UnitVector3 operator-() const;
        UnitVector3 operator+() const;

        // Binary operators
        friend bool operator==(const UnitVector3& l, const UnitVector3& r);
        friend bool operator!=(const UnitVector3& l, const UnitVector3& r);
        friend Vector3 operator*(const UnitVector3& l, const double& r);
        friend Vector3 operator*(const double& l, const UnitVector3& r);
        friend Vector3 operator/(const UnitVector3& l, const double& r);
        friend Vector3 operator+(const UnitVector3& l, const UnitVector3& r);
        friend Vector3 operator+(const Vector3& l, const UnitVector3& r);
        friend Vector3 operator+(const UnitVector3& l, const Vector3& r);
        friend Vector3 operator-(const UnitVector3& l, const UnitVector3& r);
        friend Vector3 operator-(const Vector3& l, const UnitVector3& r);
        friend Vector3 operator-(const UnitVector3& l, const Vector3& r);

        // Value lookup by index. Unlike normal array
        // lookup, this operator can not be used to
        // modify the unit vector (it returns the value,
        // not the reference).
        double operator[](unsigned i) const;

        // Standard directions
        static UnitVector3 xAxis();
        static UnitVector3 yAxis();
        static UnitVector3 zAxis();

        // Random direction uniformly distributed in the full 4 pi solid angle.
        // Needs two random numbers uniformly distributed on [0, 1).
        static UnitVector3 random(double rnd0, double rnd1);

    private:
        // Usually, it does not make sense to have a default direction.
        // Therefore, no default constructor.
        UnitVector3();
        UnitVector3(double x, double y, double z, double norm);
        void normalize();

        double x_;
        double y_;
        double z_;

        friend class Vector3;
        friend class Matrix3x3;
        friend class Rotation3;
    };

    // The "Vector3" class is used to represent vectors
    // in the 3d space
    class Vector3
    {
    public:
        // The default constructor makes a vector with
        // all components set to 0. Use with care.
        Vector3();
        Vector3(double x, double y, double z);

        // In the following constructor, negative length argument
        // is permitted. In this case the absolute value of the
        // length will be used, and the direction will be changed
        // to the opposite.
        Vector3(double length, const UnitVector3& direction);

        // Constructor which takes the difference of two points
        Vector3(const Point3& p, const Point3& origin);

        // Modifier of the vector length. Negative length argument
        // is permitted.
        Vector3& setLength(double value);

        // Modifier of the vector components
        Vector3& set(unsigned whichCoord, double value);

        // Basic accessors
        double x() const;
        double y() const;
        double z() const;

        double length() const;
        double lengthSquared() const;
        UnitVector3 direction() const;

        // Polar angle
        double theta() const;
        double cosTheta() const;

        // Azimuthal angle
        double phi() const;

        // Pseudorapidity
        double eta() const;

        // Scalar product
        double dot(const Vector3& r) const;
        double dot(const UnitVector3& r) const;

        // Cross product
        Vector3 cross(const Vector3& r) const;
        Vector3 cross(const UnitVector3& r) const;

        // Angle between two vectors
        double angle(const Vector3& r) const;
        double angle(const UnitVector3& r) const;

        // Unary operators
        Vector3 operator-() const;
        Vector3 operator+() const;

        // Binary operators
        friend bool operator==(const Vector3& l, const Vector3& r);
        friend bool operator!=(const Vector3& l, const Vector3& r);
        friend Vector3 operator*(const Vector3& l, const double& r);
        friend Vector3 operator*(const double& l, const Vector3& r);
        friend Vector3 operator/(const Vector3& l, const double& r);
        friend Vector3 operator+(const Vector3& l, const Vector3& r);
        friend Vector3 operator+(const Vector3& l, const UnitVector3& r);
        friend Vector3 operator+(const UnitVector3& l, const Vector3& r);
        friend Vector3 operator-(const Vector3& l, const Vector3& r);
        friend Vector3 operator-(const Vector3& l, const UnitVector3& r);
        friend Vector3 operator-(const UnitVector3& l, const Vector3& r);

        Vector3& operator+=(const Vector3& r);
        Vector3& operator-=(const Vector3& r);
        Vector3& operator*=(const double& d);
        Vector3& operator/=(const double& d);

        // Value lookup by index. Unlike normal array
        // lookup, this operator can not be used to
        // modify the vector (it returns the value,
        // not the reference).
        double operator[](unsigned i) const;

    private:
        double x_;
        double y_;
        double z_;
        mutable double l_;

        friend class UnitVector3;
        friend class Point3;
        friend class Matrix3x3;
        friend class Rotation3;
        friend class rk::P4;

        friend Point3 operator-(const Point3& l, const Vector3& r);
        friend Point3 operator+(const Point3& l, const Vector3& r);
        friend Point3 operator+(const Vector3& l, const Point3& r);
        friend Matrix3x3 operator*(const Matrix3x3& l, const double& r);
        friend Matrix3x3 operator/(const Matrix3x3& l, const double& r);
        friend Matrix3x3 operator*(const Matrix3x3& l, const Matrix3x3& r);
        friend Matrix3x3 operator+(const Matrix3x3& l, const Matrix3x3& r);
        friend Matrix3x3 operator-(const Matrix3x3& l, const Matrix3x3& r);
    };

    // The "Point3" class can be used to represent locations in 3d
    class Point3
    {
    public:
        // The default constructor make a point
        // with all coordinates set to 0
        Point3();
        Point3(double x, double y, double z);

        // Modifier of the components
        Point3& set(unsigned whichCoord, double value);

        double x() const;
        double y() const;
        double z() const;

        // Binary operators
        friend bool operator==(const Point3& l, const Point3& r);
        friend bool operator!=(const Point3& l, const Point3& r);
        friend Vector3 operator-(const Point3& l, const Point3& r);
        friend Point3 operator-(const Point3& l, const Vector3& r);
        friend Point3 operator+(const Point3& l, const Vector3& r);
        friend Point3 operator+(const Vector3& l, const Point3& r);

        // The following operators shift the point location by a given vector
        Point3& operator+=(const Vector3& r);
        Point3& operator-=(const Vector3& r);

        // Value lookup by index. Unlike normal array
        // lookup, this operator can not be used to
        // modify the point (it returns the value,
        // not the reference).
        double operator[](unsigned i) const;

    private:
        double x_;
        double y_;
        double z_;        

        friend class Vector3;
        friend class rk::Point4;
    };

    // The "Matrix3x3" class can be used to represent linear
    // transformations in 3d
    class Matrix3x3
    {
    public:
        // The default constructor makes a unit matrix
        Matrix3x3();

        // Constructor from a C-style array (row-by-row).
        // The array must have at least 9 elements.
        explicit Matrix3x3(const double*);

        // Explicitly specify all elements
        Matrix3x3(double m00, double m01, double m02,
                  double m10, double m11, double m12,
                  double m20, double m21, double m22);

        // The following constructors build the matrix row-by-row
        Matrix3x3(const Vector3& row0,
                  const Vector3& row1,
                  const Vector3& row2);
        Matrix3x3(const UnitVector3& row0,
                  const UnitVector3& row1,
                  const UnitVector3& row2);

        // Modifier for each element by index
        Matrix3x3& set(unsigned i, unsigned j, double value);

        // Row lookup by number. Element inspection looks
        // just like that for a 2d array: m[i][j]
        const Vector3& operator[](unsigned i) const;

        // Fill a C-style array using contents of this matrix (row-by-row).
        // The array must have at least 9 elements.
        void fillCArray(double* array) const;

        // Fill a Fortran-style array using contents of this matrix
        // (column-by-column). The array must have at least 9 elements.
        void fillFArray(double* array) const;

        // Transposed matrix
        Matrix3x3 T() const;

        // Inverse matrix
        Matrix3x3 inverse() const;

        // Determinant
        double det() const;

        // Trace
        double tr() const;

        // Unary operators
        Matrix3x3 operator-() const;
        Matrix3x3 operator+() const;

        // Binary operators
        friend bool operator==(const Matrix3x3& l, const Matrix3x3& r);
        friend bool operator!=(const Matrix3x3& l, const Matrix3x3& r);
        friend Matrix3x3 operator*(const Matrix3x3& l, const double& r);
        friend Matrix3x3 operator*(const double& l, const Matrix3x3& r);
        friend Matrix3x3 operator*(const Matrix3x3& l, const Matrix3x3& r);
        friend Vector3 operator*(const Matrix3x3& l, const Vector3& r);
        friend Vector3 operator*(const Matrix3x3& l, const UnitVector3& r);
        friend Matrix3x3 operator/(const Matrix3x3& l, const double& r);
        friend Matrix3x3 operator/(const Matrix3x3& l, const Matrix3x3& r);
        friend Matrix3x3 operator+(const Matrix3x3& l, const Matrix3x3& r);
        friend Matrix3x3 operator-(const Matrix3x3& l, const Matrix3x3& r);

        Matrix3x3& operator+=(const Matrix3x3& r);
        Matrix3x3& operator-=(const Matrix3x3& r);
        Matrix3x3& operator*=(const double);
        Matrix3x3& operator/=(const double);
        Matrix3x3& operator/=(const Matrix3x3& r);

        // The following works like this: *self = r*(*self).
        // That is, r is multiplied from the left, not
        // from the right!
        Matrix3x3& operator*=(const Matrix3x3& r);

        // Null matrix
        static Matrix3x3 null();

        // Diagonal matrix
        static Matrix3x3 diag(double m00, double m11, double m22);

    private:
        Vector3 rx_;
        Vector3 ry_;
        Vector3 rz_;

        friend class Rotation3;
    };

    // The "Rotation3" class can be used to rotate vectors in 3d.
    // For combining rotations, it works much better than representation
    // by matrices. For rotating vectors, it is slower than using matrices.
    // If you need to rotate many vectors using the same rotation object
    // and do not care very much about the effects of roundoff errors
    // on the vector length, it can be useful to obtain the matrix
    // representation first using the "matrix()" function and then perform
    // all the rotations.
    class Rotation3
    {
    public:
        // The default constructor builds the identity transformation
        // (that is, no rotation)
        Rotation3();

        Rotation3(const UnitVector3& axis, double angle);
        explicit Rotation3(const Matrix3x3& matrix);

        const UnitVector3& axis() const;
        double angle() const;

        // The following operators rotate various objects.
        // The length of the vector is preserved.
        Vector3 operator*(const Vector3&) const;
        UnitVector3 operator*(const UnitVector3&) const;

        // Two rotations are multiplied like this: r1*(r2*v) == (r1*r2)*v
        Rotation3 operator*(const Rotation3&) const;

        // Fast rotation function which does not correct the effects
        // of roundoff errors on the vector length. It is still
        // about two times slower than rotation by matrices. The "*"
        // operator is about 3.5 times slower than this function
        // (timing measurements were performed on a 64-bit Linux system).
        Vector3 rotate(const Vector3&) const;

        // The following works like this: *self = r2*(*self).
        // That is, r2 is multiplied from the left, not
        // from the right! This order is better suited
        // for implementing rotation sequences.
        Rotation3& operator*=(const Rotation3& r2);

        // Inverse rotation
        Rotation3 inverse() const;

        friend bool operator==(const Rotation3& l, const Rotation3& r);
        friend bool operator!=(const Rotation3& l, const Rotation3& r);

        // A measure of how far apart two rotations are.
        // This distance is not terribly meaningful in the 3d space.
        // However, it is symmetric, non-negative, 0 only when
        // the two rotations are identical, and it does conform
        // to the triangle inequality -- that is, this distance
        // can be used as a metric in the space of rotations.
        double distance(const Rotation3& r) const;

        // Returns a matrix which corresponds to this rotation.
        // This is sometimes useful when rotations are combined
        // with other transformations.
        Matrix3x3 matrix() const;

        // Linear rotation interpolation. Arguments are times
        // t0 and t1, rotations at times t0 and t1, and the time
        // to which the rotation should be interpolated.
        static Rotation3 interpolate(double t0, double t1,
                                     const Rotation3& r0,
                                     const Rotation3& r1,
                                     double t);

        // Cubic rotation interpolation. Arguments are times t0 and t1,
        // rotations at times t0, t0 + (t1-t0)/3, t0 + 2 (t1-t0)/3, and t1
        // (that is, rotations at the beginning, 1/3, 2/3, and the end of
        // a time interval), and the time to which the rotation should be
        // interpolated. This interpolation works best in the middle third
        // of the time interval.
        static Rotation3 interpolate(double t0, double t1,
                                     const Rotation3& r0,
                                     const Rotation3& r1_3,
                                     const Rotation3& r2_3,
                                     const Rotation3& r1,
                                     double t);

        // Random rotation. The arguments are three random numbers on [0, 1).
        static Rotation3 random(double rnd0, double rnd1, double rnd2);

    private:
        // The real role of the rotation class is to keep
        // the underlying quaternion properly normalized
        struct Quaternion
        {
            Quaternion(const UnitVector3& q, double s);
            Quaternion(const Vector3& q, double s);
            Quaternion(double x, double y, double z, double s);

            double norm() const;
            double dot(const Quaternion& r) const;
            Quaternion operator-() const;
            Quaternion operator*(const double& d) const;

            // The following function is used to normalize quaternions
            const Quaternion& normalize();

            Vector3 v_;
            double s_;
        };
        friend Quaternion operator*(const Quaternion& l,
                                    const Quaternion& r);
        friend Quaternion operator+(const Quaternion& l,
                                    const Quaternion& r);
        friend Quaternion operator-(const Quaternion& l,
                                    const Quaternion& r);

        // The following constructor assumes that
        // the quaternion is normalized
        Rotation3(const Quaternion &);

        UnitVector3 axis_;
        double angle_;
        Quaternion q_;
        Quaternion qbar_;

        friend Vector3 angularVelocity(const Rotation3&,
                                       const Rotation3&, double);
        friend Rotation3 slerp(const Rotation3&, const Rotation3&, double);

        friend class rk::LT;
    };

    // Some useful utility functions follow.
    // Conversion from degrees to radians
    inline double deg2rad(double deg) {return deg/180.0*M_PI;}

    // Conversion from radians to degrees
    inline double rad2deg(double rad) {return rad/M_PI*180.0;}

    // Average angular velocity between two orientations
    Vector3 angularVelocity(const Rotation3& r1, const Rotation3& r2,
                            double delta_t);

    // Phi difference between two directions.
    // Arguments could be either vectors or unit vectors or a mix of them.
    // The arguments can also be two doubles in which case they are treated
    // as azimuthal angles in radians.
    template <typename V1, typename V2>
    double deltaPhi(const V1& v1, const V2& v2);

    // Delta R distance in the eta-phi space (eta is pseudorapidity).
    // Arguments could be either vectors or unit vectors or a mix of them.
    template <typename V1, typename V2>
    double deltaR(const V1& v1, const V2& v2);

    // Outer product of two vectors or unit vectors or a mix of them
    template <typename V1, typename V2>
    Matrix3x3 outerProduct(const V1& v1, const V2& v2);

    // Outer product of a vector type with itself. Can be useful,
    // for example, in calculating moments of inertia.
    template <typename V>
    Matrix3x3 outerSquared(const V& v);
}

std::ostream& operator<<(std::ostream& os, const geom3::UnitVector3& u);
std::ostream& operator<<(std::ostream& os, const geom3::Vector3& v);
std::ostream& operator<<(std::ostream& os, const geom3::Rotation3& r);
std::ostream& operator<<(std::ostream& os, const geom3::Point3& p);
std::ostream& operator<<(std::ostream& os, const geom3::Matrix3x3& m);

#include "rk/geom3_UnitVector3.icc"
#include "rk/geom3_Vector3.icc"
#include "rk/geom3_Point3.icc"
#include "rk/geom3_Matrix3x3.icc"
#include "rk/geom3_Rotation3.icc"

#endif // GEOM3_GEOM3_HH_
