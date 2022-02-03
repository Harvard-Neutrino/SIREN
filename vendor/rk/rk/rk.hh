//=========================================================================
// rk.hh
//
// A collection of light-weight classes and utilities for calculations
// in relativistic kinematics. Boosts are implemented in such a way that
// particle masses are conserved exactly.
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef RK_RK_HH_
#define RK_RK_HH_

#include <utility>

#include "rk/geom3.hh"

namespace rk {
    // Forward declarations
    class Boost;

    // The "P4" class can be used to represent time-like and light-like
    // vectors in 4-dimensional Minkowski space-time. The main purpose
    // of this class is to represent 4-momenta of relativistic particles.
    //
    // With some restrictions, it can also be used to represent space-like
    // vectors (e.g., results of intermediate 4-vector subtractions).
    // However, calling the "m()" method on a space-like vector will result
    // in a run-time error. This strongly limits usability of this code
    // for representing space-like vectors (for example, such vectors
    // can not be boosted). The advantage is that the mass can be memoized
    // for time-like and light-like vectors. This leads to a significant
    // speedup in a variety of calulations.
    //
    class P4
    {
    public:
        // The default constructor makes a 4-vector
        // with all components set to 0. Use with care.
        P4();

        // Constructor from 3-momentum and mass
        P4(const geom3::Vector3& p, double m, bool isEnergyNegative=false);

        // Do not use the following constructor to make light-like
        // vectors, use the one with the explicit mass argument
        // instead. For example, construction of massless 4-vectors
        // with known p_T, eta, and phi can be accomplished as follows:
        //
        // rk::P4(p_T*geom3::Vector3(cos(phi), sin(phi), sinh(eta)), 0.0);
        //
        P4(double e, const geom3::Vector3& p);

        // Basic accessors
        const geom3::Vector3& momentum() const;
        double e() const;
        double m() const;

        // Other useful accessors with obvious meaning
        double px() const;
        double py() const;
        double pz() const;
        double p() const;

        // Magnitude of the transverse momentum
        double pt() const;

        // Transverse momentum as a 3-vector with z-component set to 0
        geom3::Vector3 transverse() const;

        // Transverse energy
        double et() const;

        // Rapidity
        double rapidity() const;

        // Pseudorapidity
        double eta() const;

        // Various angular quantities. They are the same as
        // the corresponding functions of the 3-vectors.
        double theta() const;
        double cosTheta() const;
        double phi() const;

        // Speed and gamma-factor
        double beta() const;
        double gamma() const;

        // The following function is faster and has better precision
        // than the product beta()*gamma()
        double betaGamma() const;

        // 3-velocity in units of c
        geom3::Vector3 velocity() const;

        // 4-velocity
        P4 fourVelocity() const;

        // Scalar product : E1*E2 - p1.dot(p2)
        double dot(const P4&) const;

        // Scalar product of a vector with itself. This can be negative
        // and, unlike m()*m(), can be calculated for spacelike vectors.
        double squared() const;

        // Boost this particle. If you want to get a new P4
        // object instead of boosting this one, use the "*"
        // operator of the boost object.
        P4& boost(const Boost& b);

        // Rotate the momentum of this particle. If you want to get
        // a new P4 object instead of changing this one, use the "*"
        // operator of the rotation object on the momentum of this
        // particle and construct another P4 using the result.
        P4& rotate(const geom3::Rotation3& r);

        // Lorentz transform this particle. If you want to get a new
        // P4 object instead of transforming this one, use the "*"
        // operator of the transform object.
        P4& transform(const LT& t);

        // Boost into the rest frame of this particle
        Boost restBoost() const;

        // Boost from the rest frame of this particle into the lab.
        // Equivalent to but slightly faster than restBoost().inverse().
        Boost labBoost() const;

        // Energy and momentum reversal. Use unary "-"
        // to reverse everything.
        P4 reverseE() const;
        P4 reverseP() const;

        // Unary operators
        P4 operator-() const;
        P4 operator+() const;

        // Binary operators
        friend bool operator==(const P4& l, const P4& r);
        friend bool operator!=(const P4& l, const P4& r);
        friend P4 operator+(const P4& l, const P4& r);
        friend P4 operator-(const P4& l, const P4& r);
        friend P4 operator*(const P4& l, const double& r);
        friend P4 operator*(const double& l, const P4& r);
        friend P4 operator/(const P4& l, const double& r);
        friend P4 operator*(const Boost& l, const P4& r);

        P4& operator+=(const P4& r);
        P4& operator-=(const P4& r);
        P4& operator*=(const double& d);
        P4& operator/=(const double& d);

    private:
        P4(double e, const geom3::Vector3& p, double m);
        P4(double e, const geom3::Vector3& p, bool nonNegative);
        P4(double e, const geom3::Vector3& p, double m, bool nonNegative);

        // The following function ensures m_ > 0.0 and causes
        // a run-time error if it is not
        void ensureMass_(void) const;

        geom3::Vector3 p_;
        double e_;
        mutable double m_;
        mutable bool msqIsNonNegative_;

        friend class Boost;
        friend class Point4;

        friend Point4 operator-(const Point4& l, const P4& r);
        friend Point4 operator+(const Point4& l, const P4& r);
        friend Point4 operator+(const P4& l, const Point4& r);
    };

    // The "Boost" class can be used to perform Lorentz transformations.
    // Here, the code assumes that the particle 4-momentum is known
    // in some system S and we need to know it in some system S' which
    // is moving w.r.t. S in the given direction with the speed
    // s = c*tanh(rapidity).
    //
    // The most convenient ways to create boosts are the "restBoost"
    // and "labBoost" methods of the "P4" class.
    //
    class Boost
    {
    public:
        // The default constructor builds a boost with 0 relative
        // velocity (unit transformation)
        Boost();

        explicit Boost(const geom3::Vector3& velocity_in_units_of_c);
        Boost(const geom3::UnitVector3& direction, double rapidity);

        // Basic accessors
        const geom3::UnitVector3& direction() const;
        double rapidity() const;
        double beta() const;
        double gamma() const;
        double betaGamma() const;
        geom3::Vector3 velocity() const;

        // Inverse boost
        Boost inverse() const;

        // Boosts are performed by the "*" operator. 3-vectors are
        // transformed as spatial components of lightlike 4-vectors.
        friend P4 operator*(const Boost& l, const P4& r);
        friend geom3::Vector3 operator*(const Boost& l,
                                        const geom3::Vector3& r);

        // The following operator can be used to figure out
        // how a direction changes under boost. Note that
        // the change in the vector magnitude is ignored.
        friend geom3::UnitVector3 operator*(const Boost& l,
                                            const geom3::UnitVector3& r);
        // Binary comparison operators
        friend bool operator==(const Boost& l, const Boost& r);
        friend bool operator!=(const Boost& l, const Boost& r);

    private:
        Boost(const geom3::UnitVector3& direction, double ch, double sh);

        geom3::UnitVector3 dir_;
        mutable double rapidity_;
        double c_;
        double s_;

        friend class P4;
        friend class LT;
    };

    // The "Point4" class can be used to represent space-time locations
    class Point4
    {
    public:
        Point4();
        Point4(double t, const geom3::Point3& location);

        // Basic accessors
        const geom3::Point3& location() const;
        double t() const;
        double x() const;
        double y() const;
        double z() const;

        // Binary operators
        friend bool operator==(const Point4& l, const Point4& r);
        friend bool operator!=(const Point4& l, const Point4& r);
        friend P4 operator-(const Point4& l, const Point4& r);
        friend Point4 operator-(const Point4& l, const P4& r);
        friend Point4 operator+(const Point4& l, const P4& r);
        friend Point4 operator+(const P4& l, const Point4& r);

        // The following operators shift the point location by a given vector
        Point4& operator+=(const P4& r);
        Point4& operator-=(const P4& r);

    private:
        geom3::Point3 location_;
        double t_;

        friend class P4;
    };

    // Return the velocity of a particle moving with velocity "v" (in units
    // of c) in some coordinate system S w.r.t. another system S'. The motion
    // of S' itself w.r.t. S is described by the given boost. Naturally, if
    // you already have a P4 object describing the particle, you can also
    // apply the boost to that object and then call the "velocity" method
    // of the result.
    //
    // Use the inverse boost in order to obtain the relativistic "velocity
    // addition formula".
    //
    geom3::Vector3 transformVelocity(const Boost& b, const geom3::Vector3& v);

    // Invariant mass of 2 and 3 particles. More precise
    // (but slightly slower) than adding the 4-momenta
    // and then calling the m() function.
    double invMass(const P4& p1, const P4& p2);
    double invMass(const P4& p1, const P4& p2, const P4& p3);

    // Kinematic lambda function (in the literature on relativistic
    // kinematics it is usually denoted by the square root of greek
    // lambda). This function makes certain assumptions about its usage:
    // the square root of one of the arguments must be larger than
    // the sum of the square roots of the two other.
    //
    // The center-of-mass momentum of the daughters in a two-body
    // decay is lambda(M_mother**2, M_dau1**2, M_dau2**2) / (2*M_mother).
    //
    double lambda(double x, double y, double z);

    // Random phase space particle decay into two daughters
    // with given masses (or "s-channel"). rnd0 and rnd1 should be
    // random numbers between 0 and 1. The first member of the pair
    // will have mass m1 and the second will have mass m2.
    // The mass of the parent must not be less than m1 + m2.
    std::pair<P4,P4> phaseSpaceDecay(const P4& parent, double m1, double m2,
                                     double rnd0, double rnd1);

    // The same function with a slightly different signature.
    // The results will be placed in *dau1 and *dau2.
    // Use the signature that is more appropriate for your code.
    // If your code can store the info in pairs of P4 objects and
    // can avoid unnecessary copying, the previous function will be
    // slightly faster to use.
    void phaseSpaceDecay(const P4& parent, double m1, double m2,
                         double rnd0, double rnd1, P4* dau1, P4* dau2);

    // Phase space generation via peripheral split (or "t-channel"):
    // pa + pb -> dau1 + dau2. The 4-momentum pb is allowed to be
    // spacelike (this can be useful for recursive phase space splitting),
    // but the sum pa + pb must be timelike. On exit, the 4-momentum *dau1
    // will have the mass m1 and the 4-momentum *dau2 will have the mass m2.
    //
    // For the purpose of this function, the Mandelstam t variable (the
    // square of the momentum transfer) is defined as (pa - dau1) squared.
    // t is calculated from the first random number argument, "rnd0", by
    // linear mapping within t kinematic limits. The limits and the value
    // of t are returned in *t0, *t1, and *t (assuming that the argument
    // pointers are not NULL). The kinematic limits are calculated in such
    // a way that *t0 corresponds to a small momentum transfer and *t1
    // corresponds to a large momentum transfer in a scattering process
    // in which masses of dau1, dau2 are the same as masses of pa and pb,
    // respectively. *t will be equal to *t0 for rnd0 = 0 and to *t1 for
    // rnd0 = 1. The second random number argument, "rnd1", will be used
    // to generate phi angle in the pa + pb CMS.
    void peripheralSplit(const P4& pa, const P4& pb, double m1, double m2,
                         double rnd0, double rnd1, P4* dau1, P4* dau2,
                         double* t0 = 0, double* t1 = 0, double* t = 0);
}

std::ostream& operator<<(std::ostream& os, const rk::P4& v);
std::ostream& operator<<(std::ostream& os, const rk::Point4& p);
std::ostream& operator<<(std::ostream& os, const rk::Boost& b);

#include "rk/rk_P4.icc"
#include "rk/rk_Boost.icc"
#include "rk/rk_Point4.icc"

#endif // RK_RK_HH_
