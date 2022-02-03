#include "rk/geom3.hh"

// Numerically stable implementation of sin(x)/x
inline static double sinx_over_x(const double x)
{
    if (fabs(x) < 1.0e-8)
        return 1.0;
    else
        return sin(x)/x;
}

// Large functions from the geom3::Rotation3 class are implemented below
namespace geom3 {
    Rotation3::Rotation3(const Matrix3x3& m)
        : axis_(1.0, 0.0, 0.0, false),
          angle_(0.0),
          q_(0.0, 0.0, 0.0, 1.0),
          qbar_(q_)
    {
        // Problem is, what do we do in case the input matrix
        // is not a rotation matrix? How do we make sense
        // of the arguments? What is done below is, of course,
        // somewhat arbitrary.
        //
        // 1) Assume that the transformed z direction
        //    is correct. Renormalize the z scale.
        const UnitVector3 zdir(m.rx_.z_, m.ry_.z_, m.rz_.z_);

        // 2) Renormalize the y direction and scale
        const UnitVector3 ydir(zdir.cross(Vector3(
            m.rx_.x_, m.ry_.x_, m.rz_.x_)).direction());

        // 3) Renormalize the x direction
        const UnitVector3 xdir(ydir.cross(zdir).direction());

        // 4) Convert the normalized rotation matrix (xdir, ydir, zdir)
        //    (xdir, etc. are columns) into the quaternion representation
        double x, y, z, w;
        const double t = 1.0 + xdir.x_ + ydir.y_ + zdir.z_;
        if (t > 0.01)
        {
            // Standard formula will be numerically stable
            x = ydir.z_ - zdir.y_;
            y = zdir.x_ - xdir.z_;
            z = xdir.y_ - ydir.x_;
            w = t;
        }
        else if (xdir.x_ >= ydir.y_ && xdir.x_ >= zdir.z_)
        {
            x = 1.0 + xdir.x_ - ydir.y_ - zdir.z_;
            y = ydir.x_ + xdir.y_;
            z = zdir.x_ + xdir.z_;
            w = ydir.z_ - zdir.y_;
        }
        else if (ydir.y_ >= xdir.x_ && ydir.y_ >= zdir.z_)
        {
            x = ydir.x_ + xdir.y_;
            y = 1.0 + ydir.y_ - xdir.x_ - zdir.z_;
            z = zdir.y_ + ydir.z_;
            w = zdir.x_ - xdir.z_;
        }
        else
        {
            x = zdir.x_ + xdir.z_;
            y = zdir.y_ + ydir.z_;
            z = 1.0 + zdir.z_ - xdir.x_ - ydir.y_;
            w = xdir.y_ - ydir.x_;
        }
        *this = Rotation3(Quaternion(x, y, z, w).normalize());
    }

    // The code below is, essentially, a precise and numerically
    // stable implementation of the Ken Shoemake's "slerp" quaternion
    // interpolation formula (Proceedings of SIGGRAPH 85, pp 245-254).
    // This formula works much faster than determination of the average
    // angular velocity and subsequent rotation from the start time
    // to time t with that velocity.
    Rotation3 Rotation3::interpolate(const double t1, const double t2,
                                     const Rotation3& r1, const Rotation3& r2,
                                     double t)
    {
        if (t1 == t2)
        {
            assert(t == t2);
            assert(r1 == r2);
            return r1;
        }
        else
        {
            t = (t - t1)/(t2 - t1);
            double cosa = r1.q_.dot(r2.q_);
            Rotation3::Quaternion q2(cosa >= 0.0 ? r2.q_ : -r2.q_);
            cosa = fabs(cosa);
            double a;
            if (cosa < 0.99)
                a = acos(cosa);
            else
                a = 2.0*asin((r1.q_ - q2).norm()/2.0);
            const double onemt = 1.0 - t;
            const double c0 = sinx_over_x(a);
            const double c1 = onemt*sinx_over_x(a*onemt)/c0;
            const double c2 = t*sinx_over_x(a*t)/c0;
            return Rotation3(r1.q_*c1 + q2*c2);
        }
    }

    // The code below is different from the popular "squad" formula
    // of Shoemake. Shoemake's formula is a Bezier curve. This one
    // passes through the control points.
    Rotation3 Rotation3::interpolate(const double t0, const double t1,
                                     const Rotation3& r0, const Rotation3& r1_3,
                                     const Rotation3& r2_3, const Rotation3& r1,
                                     double t)
    {
        if (t0 == t1)
        {
            assert(t == t1);
            assert(r0 == r1);
            assert(r0 == r1_3);
            assert(r0 == r2_3);
            return r1;
        }
        else
        {
            t = (t - t0)/(t1 - t0);
            return interpolate(0.0, 1.0, interpolate(0.0,1.0,r0,r1,t),
                               interpolate(0.0,1.0,r1_3,r2_3,3.0*t-1.0),
                               4.5*t*(1.0-t));
        }
    }

    Rotation3& Rotation3::operator*=(const Rotation3& r)
    {
        q_ = (r.q_*q_).normalize();
        qbar_.v_ = -q_.v_;
        qbar_.s_ = q_.s_;
        axis_ = q_.v_.direction();
        angle_ = 2.0*atan2(q_.v_.length(), q_.s_);
        return *this;
    }

    // See Ken Shoemake's article "Uniform Random Rotations" in
    // "Graphics Gems III", section III.6.
    Rotation3 Rotation3::random(const double rnd0, const double rnd1,
                                const double rnd2)
    {
        assert(rnd0 >= 0.0 && rnd0 <= 1.0);
        const double theta1(2*M_PI*rnd1);
        const double theta2(2*M_PI*rnd2);
        const double r1(sqrt(1.0 - rnd0));
        const double r2(sqrt(rnd0));
        return Rotation3(Rotation3::Quaternion(
                             r1*sin(theta1), r1*cos(theta1),
                             r2*sin(theta2), r2*cos(theta2)));
    }
}

// Large functions from the geom3::UnitVector3 class
namespace geom3 {
    double UnitVector3::angle(const UnitVector3& r) const
    {
        const double cosangle = dot(r);
        if (fabs(cosangle) < 0.99)
            return acos(cosangle);
        else
        {
	    /* acos would loose too much precision */
            if (cosangle > 0.0)
                return 2.0*asin((*this - r).length()/2.0);
            else
                return M_PI - 2.0*asin((*this + r).length()/2.0);
        }
    }

    double UnitVector3::theta() const
    {
        if (fabs(z_) < 0.99)
            return acos(z_);
        else
        {
	    /* acos would loose too much precision */
            const double t0 = asin(sqrt(x_*x_ + y_*y_));
            if (z_ > 0.0)
                return t0;
            else
                return M_PI - t0;
        }
    }

    double UnitVector3::operator[](const unsigned i) const
    {
        switch (i)
        {
        case X:
            return x_;
        case Y:
            return y_;
        case Z:
            return z_;
        default:
            assert(!"geom3::UnitVector3::[] index out of range");
            return 0.0;
        }
    }

    UnitVector3 UnitVector3::random(const double rnd0, const double rnd1)
    {
        const double cosTheta(rnd0*2.0 - 1.0);
        const double sinThetaSq(1.0 - cosTheta*cosTheta);
        assert(sinThetaSq >= 0.0 && sinThetaSq <= 1.0);
        const double s(sqrt(sinThetaSq));
        const double phi(2.0*M_PI*rnd1);
        return UnitVector3(s*cos(phi), s*sin(phi), cosTheta, false);
    }

    UnitVector3& UnitVector3::setTheta(const double theta)
    {
        assert(theta >= 0.0 && theta <= M_PI);
        const double t(sqrt(x_*x_ + y_*y_));
        if (t > 0.0)
        {
            double ts(sin(theta)/t);
            x_ *= ts;
            y_ *= ts;
        }
        else
        {
            x_ = sin(theta);
            y_ = 0.0;
        }
        z_ = cos(theta);
        return *this;
    }

    UnitVector3& UnitVector3::setEta(const double eta)
    {
        const double t(sqrt(x_*x_ + y_*y_));
        const double z_over_t(sinh(eta));
        const double s(1.0/sqrt(1.0 + z_over_t*z_over_t));
        if (t > 0.0)
        {
            const double ts(s/t);
            x_ *= ts;
            y_ *= ts;
        }
        else
        {
            x_ = s;
            y_ = 0.0;
        }
        z_ = s*z_over_t;
        return *this;
    }
}

// Large functions from the geom3::Vector3 class
namespace geom3 {
    double Vector3::theta() const
    {
        if (l_ < 0.0)
            l_ = sqrt(x_*x_ + y_*y_ + z_*z_);
        if (l_ == 0.0)
            return M_PI/2.0;
        else
        {
            const double cosangle = z_/l_;
            if (fabs(cosangle) < 0.99)
                return acos(cosangle);
            else
            {
                /* acos would loose too much precision */
                const double t0 = asin(sqrt(x_*x_ + y_*y_)/l_);
                if (z_ > 0.0)
                    return t0;
                else
                    return M_PI - t0;
            }
        }
    }

    Vector3& Vector3::setLength(double value)
    {
        if (l_ < 0.0)
            l_ = sqrt(x_*x_ + y_*y_ + z_*z_);
        if (l_ > 0.0)
        {
            const double rat(value/l_);
            x_ *= rat;
            y_ *= rat;
            z_ *= rat;
        }
        else
        {
            x_ = value;
            y_ = 0.0;
            z_ = 0.0;
        }
        l_ = fabs(value);
        return *this;
    }

    double Vector3::operator[](const unsigned i) const
    {
        switch (i)
        {
        case X:
            return x_;
        case Y:
            return y_;
        case Z:
            return z_;
        default:
            assert(!"geom3::Vector3::[] index out of range");
            return 0.0;
        }
    }

    Vector3& Vector3::set(const unsigned i, const double value)
    {
        switch (i)
        {
        case X:
            x_ = value;
            break;
        case Y:
            y_ = value;
            break;
        case Z:
            z_ = value;
            break;
        default:
            assert(!"geom3::Vector3::set index out of range");
        }
        l_ = -1.0;
        return *this;
    }
}

// Large functions from the geom3::Matrix3x3 class
namespace geom3 {
    Matrix3x3 Matrix3x3::inverse() const
    {
        const double d(det());
        assert(d != 0.0);
        return Matrix3x3(ry_.cross(rz_)/=d,
                         rz_.cross(rx_)/=d,
                         rx_.cross(ry_)/=d).T();
    }

    const Vector3& Matrix3x3::operator[](const unsigned i) const
    {
        switch (i)
        {
        case X:
            return rx_;
        case Y:
            return ry_;
        case Z:
            return rz_;
        default:
            assert(!"geom3::Matrix3x3::[] index out of range");
            return rz_;
        }
    }

    Matrix3x3& Matrix3x3::set(const unsigned i, const unsigned j,
                              const double value)
    {
        switch (i)
        {
        case X:
            rx_.set(j, value);
            break;
        case Y:
            ry_.set(j, value);
            break;
        case Z:
            rz_.set(j, value);
            break;
        default:
            assert(!"geom3::Matrix3x3::set index out of range");
        }
        return *this;
    }
}

// Large functions from the geom3::Point3 class
namespace geom3 {
    double Point3::operator[](const unsigned i) const
    {
        switch (i)
        {
        case X:
            return x_;
        case Y:
            return y_;
        case Z:
            return z_;
        default:
            assert(!"geom3::Point3::[] index out of range");
            return 0.0;
        }
    }

    Point3& Point3::set(const unsigned i, const double value)
    {
        switch (i)
        {
        case X:
            x_ = value;
            break;
        case Y:
            y_ = value;
            break;
        case Z:
            z_ = value;
            break;
        default:
            assert(!"geom3::Point3::set index out of range");
        }
        return *this;
    }
}

// Large functions from the geom3 namespace which do not belong
// to any class
namespace geom3 {
    Vector3 angularVelocity(const Rotation3& r1, const Rotation3& r2,
                            const double dt)
    {
        assert(dt != 0.0);
        Rotation3::Quaternion q(r1.q_.dot(r2.q_) >= 0.0 ? 
                                r2.q_*r1.qbar_ : 
                                (-r2.q_)*r1.qbar_);
        q.normalize();
        return Vector3(2.0*atan2(q.v_.length(), q.s_)/dt, q.v_.direction());
    }
}
