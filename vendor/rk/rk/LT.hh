//=========================================================================
// LT.hh
//
// Representation of arbitrary Lorentz transformations using biquaternions.
// Particle masses are conserved exactly in these transformations and
// a combination of Lorentz transformations is a Lorentz transformation
// (that is, round-off errors are corrected).
//
// I. Volobouev
// February 2013
//=========================================================================

#ifndef RK_LT_HH_
#define RK_LT_HH_

#include <complex>

#include "rk/rk.hh"

namespace rk {
    class LT
    {
    public:
        // The default constructor builds the identity transformation
        LT();

        // Constructors from rotations and boosts
        LT(const geom3::Rotation3& r);
        LT(const Boost& b);

        // The following operator transforms 4-vectors
        P4 operator*(const P4&) const;

        // Two transforms are multiplied like this: t1*(t2*P4) == (t1*t2)*P4
        LT operator*(const LT&) const;

        // The following works like this: *self = t2*(*self).
        // That is, t2 is multiplied from the left, not
        // from the right! This order is better suited
        // for implementing transformation sequences.
        LT& operator*=(const LT& t2);

        // Inverse transformation
        LT inverse() const;

        friend bool operator==(const LT& l, const LT& r);
        friend bool operator!=(const LT& l, const LT& r);

        // A measure of how far apart two transformations are.
        // This distance is not terribly meaningful in the 4-d space.
        // However, it is symmetric, non-negative, 0 only when
        // the two transformations are identical, and it does conform
        // to the triangle inequality -- that is, this distance can be
        // used as a metric in the space of Lorentz transformations.
        double distance(const LT& r) const;

        // Decompose the general transformation into a boost followed
        // by a rotation. Note the order of arguments -- it is as in the
        // corresponding product of transformations. Any of the argument
        // pointers can be NULL in which case the corresponding result
        // is not filled.
        void decompose(geom3::Rotation3* r, Boost* b) const;

        // Decompose the general transformation into a rotation followed
        // by a boost. Any of the argument pointers can be NULL in
        // which case the corresponding result is not filled.
        void decompose(Boost* b, geom3::Rotation3* r) const;

    private:
        typedef std::complex<double> Complex;

        struct Biquaternion
        {
            Biquaternion(const Complex& q0, const Complex& q1,
                         const Complex& q2, const Complex& q3);

            // Special constructor from a boost. It is more efficient
            // to put it here instead of making it from the LT class.
            Biquaternion(const Boost& b);

            // Unary minus
            Biquaternion operator-() const;

            // Quaternion conjugate
            Biquaternion biconjugate() const;

            // Complex conjugate
            Biquaternion conj() const;

            // Hermitian conjugate (biconjugate followed by complex conjugate)
            Biquaternion hconjugate() const;

            // Biquaternion "length"
            double length() const;

            // Normalize this biquaternion (so that it represents a valid
            // Lorentz transformation)
            const Biquaternion& normalize();

            // Extract 3-momentum of a 4-vector from a biquaternion
            // representation. Then 4-momentum can be constructed
            // either using mass (if known) or energy.
            geom3::Vector3 momentum() const;

            // Extract energy of a 4-vector from a biquaternion
            // representation
            double e() const;

            Complex q0_, q1_, q2_, q3_;
        };

        friend Biquaternion operator*(const Biquaternion& l,
                                      const Biquaternion& r);
        friend Biquaternion operator+(const Biquaternion& l,
                                      const Biquaternion& r);
        friend Biquaternion operator-(const Biquaternion& l,
                                      const Biquaternion& r);
        friend bool operator==(const Biquaternion& l,
                               const Biquaternion& r);
        friend bool operator!=(const Biquaternion& l,
                               const Biquaternion& r);

        // The following constructor assumes that
        // the biquaternion is normalized
        LT(const Biquaternion&);

        Biquaternion q_;
        mutable Biquaternion qdagger_;
        mutable bool updated_;

    public:
        // Code to enable serialization. Typically,
        // should not be used directly by applications.
        // Use "rkIO.hh" header instead with "Geners".
        void save(double data[8]) const;
        LT& load(const double data[8]);
    };
}

#include "rk/rk_LT.icc"

#endif // RK_LT_HH_
