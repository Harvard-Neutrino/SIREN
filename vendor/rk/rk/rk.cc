#include "rk/LT.hh"

//
// Large functions (or functions whose implementation can be expected
// to change) from the rk namespace are implemented below
//
namespace rk {
    P4& P4::boost(const Boost& b)
    {
        if (m_ < 0.0)
        {
            const double msq(e_*e_ - p_.lengthSquared());
            assert(msq >= 0.0);
            msqIsNonNegative_ = true;
            m_ = sqrt(msq);
        }
        const double par(b.dir_.dot(p_));
        p_ += b.dir_*(b.c_*par - e_*b.s_);
        e_ = sqrt(m_*m_ + p_.lengthSquared()) * (e_ < 0.0 ? -1.0 : 1.0);
        return *this;
    }

    P4& P4::rotate(const geom3::Rotation3& r)
    {
        p_ = r*p_;
        return *this;
    }

    P4& P4::transform(const LT& t)
    {
        *this = (t * *this);
        return *this;
    }

    LT::Biquaternion::Biquaternion(const Boost& b)
    {
        const double halfAngleCh = sqrt(1.0 + b.c_/2.0);
        const double mHalfAngleSh = -b.s_/halfAngleCh/2.0;
        q0_ = Complex(halfAngleCh, 0.0);
        q1_ = Complex(0.0, mHalfAngleSh*b.dir_.x());
        q2_ = Complex(0.0, mHalfAngleSh*b.dir_.y());
        q3_ = Complex(0.0, mHalfAngleSh*b.dir_.z());
    }

    void LT::save(double data[8]) const
    {
        data[0] = q_.q0_.real();
        data[1] = q_.q0_.imag();
        data[2] = q_.q1_.real();
        data[3] = q_.q1_.imag();
        data[4] = q_.q2_.real();
        data[5] = q_.q2_.imag();
        data[6] = q_.q3_.real();
        data[7] = q_.q3_.imag();
    }

    LT& LT::load(const double data[8])
    {
        q_ = Biquaternion(Complex(data[0], data[1]), Complex(data[2], data[3]),
                          Complex(data[4], data[5]), Complex(data[6], data[7]));
        updated_ = false;
        return *this;
    }

    LT::Biquaternion operator*(const LT::Biquaternion& l,
                               const LT::Biquaternion& r)
    {
        return LT::Biquaternion(
            l.q0_*r.q0_ - l.q1_*r.q1_ - l.q2_*r.q2_ - l.q3_*r.q3_,
            r.q0_*l.q1_ + r.q1_*l.q0_ + r.q3_*l.q2_ - r.q2_*l.q3_,
            r.q0_*l.q2_ + r.q1_*l.q3_ + r.q2_*l.q0_ - r.q3_*l.q1_,
            r.q0_*l.q3_ + r.q2_*l.q1_ + r.q3_*l.q0_ - r.q1_*l.q2_);
    }

    P4 LT::operator*(const P4& p) const
    {
        if (!updated_)
        {
            qdagger_ = q_.hconjugate();
            updated_ = true;
        }
        const Biquaternion& pbq = Biquaternion(Complex(p.e(), 0.0),
                                               Complex(0.0, p.px()),
                                               Complex(0.0, p.py()),
                                               Complex(0.0, p.pz()));
        return P4((q_*pbq*qdagger_).momentum(), p.m(), p.e() < 0.0);
    }

    const LT::Biquaternion& LT::Biquaternion::normalize()
    {
        const double rq0 = q0_.real();
        const double iq0 = q0_.imag();
        const double rq1 = q1_.real();
        const double iq1 = q1_.imag();
        const double rq2 = q2_.real();
        const double iq2 = q2_.imag();
        const double rq3 = q3_.real();
        const double iq3 = q3_.imag();

        // We have to ensure two things in this normalization:
        //
        // 1)  rq0*iq0 + rq1*iq1 + rq2*iq2 + rq3*iq3 = 0
        //
        // 2)  rq0*rq0 + rq1*rq1 + rq2*rq2 + rq3*rq3 - 
        //     (iq0*iq0 + iq1*iq1 + iq2*iq2 + iq3*iq3) = 0
        //
        // The procedure which follows is, of course, somewhat
        // arbitrary. We are going to adjust iqN first by rotating
        // it in 4-d to satisfy the first condition and then we
        // will scale all rqN values to satisfy the second.
        //
        const double rnormsq = rq0*rq0 + rq1*rq1 + rq2*rq2 + rq3*rq3;
        assert(rnormsq > 0.0);
        const double rat = (rq0*iq0 + rq1*iq1 + rq2*iq2 + rq3*iq3)/rnormsq;
        const double if0 = iq0 - rat*rq0;
        const double if1 = iq1 - rat*rq1;
        const double if2 = iq2 - rat*rq2;
        const double if3 = iq3 - rat*rq3;

        const double inormsq = iq0*iq0 + iq1*iq1 + iq2*iq2 + iq3*iq3;
        const double rfactor = sqrt((inormsq + 1.0)/rnormsq);

        const double fnormsq = if0*if0 + if1*if1 + if2*if2 + if3*if3;
        const double ffactor = fnormsq > 0.0 ? sqrt(inormsq/fnormsq) : 1.0;

        q0_ = Complex(rq0*rfactor, if0*ffactor);
        q1_ = Complex(rq1*rfactor, if1*ffactor);
        q2_ = Complex(rq2*rfactor, if2*ffactor);
        q3_ = Complex(rq3*rfactor, if3*ffactor);

        return *this;
    }

    // Boost followed by a rotation
    void LT::decompose(geom3::Rotation3* r, Boost* b) const
    {
        if (r == 0 && b == 0)
            return;

        geom3::Rotation3::Quaternion q(q_.q1_.real(), q_.q2_.real(),
                                       q_.q3_.real(), q_.q0_.real());
        q.normalize();
        if (r)
            *r = geom3::Rotation3(q);

        if (b)
        {
            Biquaternion qbar(Complex(q.s_, 0.0), Complex(-q.v_.x(), 0.0),
                              Complex(-q.v_.y(), 0.0), Complex(-q.v_.z(), 0.0));
            const Biquaternion& qb = qbar*q_;
            const double iq0 = qb.q0_.imag();
            const double iq1 = -qb.q1_.imag();
            const double iq2 = -qb.q2_.imag();
            const double iq3 = -qb.q3_.imag();
            const double inormsq = iq0*iq0 + iq1*iq1 + iq2*iq2 + iq3*iq3;
            const double halfAngleSh = sqrt(inormsq);
            if (halfAngleSh > 0.0)
                *b = Boost(geom3::UnitVector3(iq1, iq2, iq3),
                           2.0*asinh(halfAngleSh));
            else
                *b = Boost();
        }
    }

    // Rotation followed by a boost
    void LT::decompose(Boost* pb, geom3::Rotation3* r) const
    {
        if (r == 0 && pb == 0)
            return;

        geom3::Rotation3::Quaternion q(q_.q1_.real(), q_.q2_.real(),
                                       q_.q3_.real(), q_.q0_.real());
        q.normalize();
        Biquaternion qbar(Complex(q.s_, 0.0), Complex(-q.v_.x(), 0.0),
                          Complex(-q.v_.y(), 0.0), Complex(-q.v_.z(), 0.0));
        const Biquaternion& qb = q_*qbar;
        const double iq0 = qb.q0_.imag();
        const double iq1 = -qb.q1_.imag();
        const double iq2 = -qb.q2_.imag();
        const double iq3 = -qb.q3_.imag();
        const double inormsq = iq0*iq0 + iq1*iq1 + iq2*iq2 + iq3*iq3;
        const double halfAngleSh = sqrt(inormsq);
        const Boost& b = halfAngleSh ? Boost(geom3::UnitVector3(iq1, iq2, iq3),
                                             2.0*asinh(halfAngleSh)) : Boost();
        if (pb)
            *pb = b;

        if (r)
        {
            const Biquaternion& qr = Biquaternion(b.inverse())*q_;
            geom3::Rotation3::Quaternion q2(qr.q1_.real(), qr.q2_.real(),
                                            qr.q3_.real(), qr.q0_.real());
            q2.normalize();
            *r = geom3::Rotation3(q2);
        }
    }

    double lambda(const double x, const double y, const double z)
    {
        if (x == 0.0)
            return fabs(y - z);
        else if (y == 0.0)
            return fabs(x - z);
        else if (z == 0.0)
            return fabs(x - y);
        else
        {
            const double ymz(y - z);
            const double dtmp(ymz*ymz + x*(x - 2.0*(y + z)));
            assert(dtmp >= 0.0);
            return sqrt(dtmp);
        }
    }

    std::pair<P4,P4> phaseSpaceDecay(const P4& parent,
                                     const double m1,
                                     const double m2,
                                     const double rnd0,
                                     const double rnd1)
    {
        assert(m1 >= 0.0 && m2 >= 0.0);
        const double parentM(parent.m());
        assert(parentM >= m1 + m2);
        const Boost& lb = parent.labBoost();
        if (parentM == m1 + m2)
            return std::pair<P4,P4>(lb*P4(geom3::Vector3(), m1),
                                    lb*P4(geom3::Vector3(), m2));
        else
        {
            const geom3::UnitVector3& dir = 
                geom3::UnitVector3::random(rnd0, rnd1);
            const double p(lambda(parentM*parentM, m1*m1, m2*m2)/2.0/parentM);
            return std::pair<P4,P4>(lb*P4(dir*p, m1),
                                    lb*P4(dir*(-p), m2));            
        }
    }

    void phaseSpaceDecay(const P4& parent, const double m1, const double m2,
                         const double rnd0, const double rnd1,
                         P4* dau1, P4* dau2)
    {
        assert(dau1);
        assert(dau2);
        assert(m1 >= 0.0 && m2 >= 0.0);
        const double parentM(parent.m());
        assert(parentM >= m1 + m2);
        const Boost& lb = parent.labBoost();
        if (parentM == m1 + m2)
        {
            *dau1 = lb*P4(geom3::Vector3(), m1);
            *dau2 = lb*P4(geom3::Vector3(), m2);
        }
        else
        {
            const geom3::UnitVector3& dir = 
                geom3::UnitVector3::random(rnd0, rnd1);
            const double p(lambda(parentM*parentM, m1*m1, m2*m2)/2.0/parentM);
            *dau1 = lb*P4(dir*p, m1);
            *dau2 = lb*P4(dir*(-p), m2);
        }
    }

    void peripheralSplit(const P4& pa, const P4& pb,
                         const double m1, const double m2,
                         const double rnd0, const double rnd1,
                         P4* dau1, P4* dau2,
                         double* pt0, double *pt1, double* pt)
    {
        assert(dau1);
        assert(dau2);
        assert(m1 >= 0.0 && m2 >= 0.0);
        const P4& parent = pa + pb;
        const double parentM(parent.m());
        assert(parentM >= m1 + m2);
        const double s(parentM*parentM);
        const double ma(pa.m());
        const double mbsq(pb.squared());
        const double tmidpoint = ma*ma + m1*m1 - 
            (s + ma*ma - mbsq)*(s + m1*m1 - m2*m2)/2.0/s;
        double tdelta = 0.0;
        const Boost& lb = parent.labBoost();

        // Make random cosTheta so that random value 0 maps
        // into the cosine value of 1
        const double cosTheta = 1.0 - rnd0*2.0;

        // Check for various degenerate situations
        if (parentM == m1 + m2)
        {
            *dau1 = lb*P4(geom3::Vector3(), m1);
            *dau2 = lb*P4(geom3::Vector3(), m2);
        }
        else if (mbsq >= 0.0 && parentM <= ma + pb.m())
        {
            // Momentum of particle a is 0 in the CMS,
            // so the Z axis can be chosen as the direction
            // of this particle (just as well as any other axis).
            // Then we are back to the s-channel case.
            const geom3::UnitVector3& dir = 
                geom3::UnitVector3::random(rnd0, rnd1);
            const double p(lambda(s, m1*m1, m2*m2)/2.0/parentM);
            *dau1 = lb*P4(dir*p, m1);
            *dau2 = lb*P4(dir*(-p), m2);
        }
        else
        {
            // Perhaps, the simplest way to do the generation
            // is to figure out the direction of the pa 3-momentum
            // in the CMS. In the CMS, t is related linearly to the
            // cosine of the angle between pa and dau1.
            const P4& pa4cms = lb.inverse()*pa;
            const geom3::UnitVector3& zdir = pa4cms.momentum().direction();

            // We now need to come up with an orthonormal system in which
            // "zdir" works as the Z axis
            const geom3::UnitVector3& xdir = 
                fabs(zdir.x()) > fabs(zdir.y()) ?
                zdir.cross(geom3::UnitVector3::yAxis()).direction() :
                zdir.cross(geom3::UnitVector3::xAxis()).direction();
            const geom3::Vector3& ydir = zdir.cross(xdir);

            // Magnitudes of pa and dau1 momenta in the CMS
            const double pacms = pa4cms.p();
            double p1cms = 0.0;
            if (mbsq >= 0.0 && m1 == ma && m2 == pb.m())
            {
                p1cms = pacms;
                tdelta = -tmidpoint;
            }
            else
            {
                p1cms = lambda(s, m1*m1, m2*m2)/2.0/parentM;
                tdelta = 2.0*pacms*p1cms;
            }

            // Generate px, py, pz for dau1 in the rotated system
            const double sinThetaSq(1.0 - cosTheta*cosTheta);
            assert(sinThetaSq >= 0.0 && sinThetaSq <= 1.0);
            const double sinTheta(sqrt(sinThetaSq));
            const double phi(2.0*M_PI*rnd1);
            const double px = p1cms*sinTheta*cos(phi);
            const double py = p1cms*sinTheta*sin(phi);
            const double pz = p1cms*cosTheta;

            // Figure out the dau1 momentum in the original system
            geom3::Vector3 p1mom(xdir*px + ydir*py + zdir*pz);

            // Boost daus back to the lab
            *dau1 = lb*P4(p1mom, m1);
            p1mom *= -1.0;
            *dau2 = lb*P4(p1mom, m2);
        }

        // Return Mandelstam t and integration limits
        if (pt)
            *pt = tmidpoint + tdelta*cosTheta;
        if (pt0)
            *pt0 = tmidpoint + tdelta;
        if (pt1)
            *pt1 = tmidpoint - tdelta;
    }

    geom3::Vector3 transformVelocity(const Boost& b, const geom3::Vector3& v)
    {
        assert(v.length() <= 1.0);
        const geom3::Vector3& vb = b.velocity();
        const geom3::Vector3& vpar = b.direction()*(b.direction().dot(v));
        return ((v - vpar)/b.gamma() + vpar - vb)/(1.0 - vb.dot(v));
    }
}
