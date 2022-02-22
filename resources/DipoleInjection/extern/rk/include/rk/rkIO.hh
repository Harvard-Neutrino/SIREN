//=========================================================================
// rkIO.hh
//
// The following wrapper code makes classes in the "rk" package
// work with generic I/O mechanisms of "Geners" serialization
// software (see http://geners.hepforge.org/).
//
// I. Volobouev
// April 2012
//=========================================================================

#ifndef RK_RKIO_HH_
#define RK_RKIO_HH_

#include "rk/LT.hh"
#include "geners/GenericIO.hh"

gs_declare_type_external(geom3::Point3)
gs_declare_type_external(geom3::UnitVector3)
gs_declare_type_external(geom3::Vector3)
gs_declare_type_external(geom3::Rotation3)
gs_declare_type_external(geom3::Matrix3x3)
gs_declare_type_external(rk::P4)
gs_declare_type_external(rk::Boost)
gs_declare_type_external(rk::Point4)
gs_declare_type_external(rk::LT)

gs_specialize_class_id(geom3::Point3, 1)
gs_specialize_class_id(geom3::UnitVector3, 1)
gs_specialize_class_id(geom3::Vector3, 1)
gs_specialize_class_id(geom3::Rotation3, 1)
gs_specialize_class_id(geom3::Matrix3x3, 1)
gs_specialize_class_id(rk::P4, 1)
gs_specialize_class_id(rk::Boost, 1)
gs_specialize_class_id(rk::Point4, 1)
gs_specialize_class_id(rk::LT, 1)

namespace gs {
    template <class Stream, class State>
    struct GenericWriter<Stream, State, geom3::Point3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const geom3::Point3& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<geom3::Point3>().write(os) : true;
            if (status)
            {
                write_pod(os, s.x());
                write_pod(os, s.y());
                write_pod(os, s.z());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, geom3::UnitVector3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const geom3::UnitVector3& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<geom3::UnitVector3>().write(os) : true;
            if (status)
            {
                write_pod(os, s.x());
                write_pod(os, s.y());
                write_pod(os, s.z());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, geom3::Vector3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const geom3::Vector3& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<geom3::Vector3>().write(os) : true;
            if (status)
            {
                write_pod(os, s.x());
                write_pod(os, s.y());
                write_pod(os, s.z());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, geom3::Rotation3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const geom3::Rotation3& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<geom3::Rotation3>().write(os) : true;
            if (status)
            {
                const geom3::UnitVector3& a = s.axis();
                write_pod(os, a.x());
                write_pod(os, a.y());
                write_pod(os, a.z());
                write_pod(os, s.angle());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, geom3::Matrix3x3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const geom3::Matrix3x3& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<geom3::Matrix3x3>().write(os) : true;
            if (status)
            {
                for (unsigned i=0; i<3; ++i)
                    for (unsigned j=0; j<3; ++j)
                    {
                        double d = s[i][j];
                        write_pod(os, d);
                    }
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, rk::P4,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const rk::P4& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<rk::P4>().write(os) : true;
            if (status)
            {
                const geom3::Vector3& p = s.momentum();
                write_pod(os, p.x());
                write_pod(os, p.y());
                write_pod(os, p.z());
                write_pod(os, s.m());
                unsigned char negativee = (s.e() < 0.0);
                write_pod(os, negativee);
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, rk::Boost,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const rk::Boost& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<rk::Boost>().write(os) : true;
            if (status)
            {
                const geom3::UnitVector3& dir = s.direction();
                write_pod(os, dir.x());
                write_pod(os, dir.y());
                write_pod(os, dir.z());
                write_pod(os, s.rapidity());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, rk::LT,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const rk::LT& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<rk::LT>().write(os) : true;
            if (status)
            {
                double data[8];
                s.save(data);
                write_pod_array(os, &data[0], 8);
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericWriter<Stream, State, rk::Point4,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool process(const rk::Point4& s, Stream& os,
                                   State*, const bool processClassId)
        {
            const bool status = processClassId ? 
                ClassId::makeId<rk::Point4>().write(os) : true;
            if (status)
            {
                write_pod(os, s.t());
                write_pod(os, s.x());
                write_pod(os, s.y());
                write_pod(os, s.z());
            }
            return status && !os.fail();
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, geom3::Point3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(geom3::Point3*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<geom3::Point3>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new geom3::Point3(x, y, z);
            else
                *ptr = geom3::Point3(x, y, z);
            return true;
        }

        inline static bool process(geom3::Point3& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            geom3::Point3* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, geom3::Vector3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(geom3::Vector3*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<geom3::Vector3>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new geom3::Vector3(x, y, z);
            else
                *ptr = geom3::Vector3(x, y, z);
            return true;
        }

        inline static bool process(geom3::Vector3& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            geom3::Vector3* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, geom3::UnitVector3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(geom3::UnitVector3*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<geom3::UnitVector3>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new geom3::UnitVector3(x, y, z, false);
            else
                *ptr = geom3::UnitVector3(x, y, z, false);
            return true;
        }

        inline static bool process(geom3::UnitVector3& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            geom3::UnitVector3* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, geom3::Rotation3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(geom3::Rotation3*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<geom3::Rotation3>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0, ang=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            read_pod(is, &ang);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new geom3::Rotation3(geom3::UnitVector3(
                                               x, y, z, false), ang);
            else
                *ptr = geom3::Rotation3(geom3::UnitVector3(
                                            x, y, z, false), ang);
            return true;
        }

        inline static bool process(geom3::Rotation3& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            geom3::Rotation3* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, geom3::Matrix3x3,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(geom3::Matrix3x3*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<geom3::Matrix3x3>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double m00=0, m01=0, m02=0,
                   m10=0, m11=0, m12=0,
                   m20=0, m21=0, m22=0;
            read_pod(is, &m00);
            read_pod(is, &m01);
            read_pod(is, &m02);
            read_pod(is, &m10);
            read_pod(is, &m11);
            read_pod(is, &m12);
            read_pod(is, &m20);
            read_pod(is, &m21);
            read_pod(is, &m22);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new geom3::Matrix3x3(m00, m01, m02, m10, m11,
                                           m12, m20, m21, m22);
            else
                *ptr = geom3::Matrix3x3(m00, m01, m02, m10, m11,
                                        m12, m20, m21, m22);
            return true;
        }

        inline static bool process(geom3::Matrix3x3& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            geom3::Matrix3x3* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, rk::P4,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(rk::P4*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<rk::P4>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0, m=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            read_pod(is, &m);
            unsigned char negativee = 0;
            read_pod(is, &negativee);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new rk::P4(geom3::Vector3(x, y, z), m, negativee);
            else
                *ptr = rk::P4(geom3::Vector3(x, y, z), m, negativee);
            return true;
        }

        inline static bool process(rk::P4& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            rk::P4* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, rk::Point4,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(rk::Point4*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<rk::Point4>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0, t=0;
            read_pod(is, &t);
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new rk::Point4(t, geom3::Point3(x, y, z));
            else
                *ptr = rk::Point4(t, geom3::Point3(x, y, z));
            return true;
        }

        inline static bool process(rk::Point4& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            rk::Point4* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, rk::Boost,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(rk::Boost*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<rk::Boost>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double x=0, y=0, z=0, rap=0;
            read_pod(is, &x);
            read_pod(is, &y);
            read_pod(is, &z);
            read_pod(is, &rap);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new rk::Boost(geom3::UnitVector3(x, y, z, false), rap);
            else
                *ptr = rk::Boost(geom3::UnitVector3(x, y, z, false), rap);
            return true;
        }

        inline static bool process(rk::Boost& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            rk::Boost* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };

    template <class Stream, class State>
    struct GenericReader<Stream, State, rk::LT,
                         Int2Type<IOTraits<int>::ISEXTERNAL> >
    {
        inline static bool readIntoPtr(rk::LT*& ptr, Stream& is,
                                       State*, const bool processClassId)
        {
            if (processClassId)
            {
                static const ClassId current(ClassId::makeId<rk::LT>());
                ClassId id(is, 1);
                current.ensureSameName(id);
            }

            double data[8];
            read_pod_array(is, &data[0], 8);
            if (is.fail())
                return false;

            if (ptr == 0)
                ptr = new rk::LT();
            ptr->load(data);
            return true;
        }

        inline static bool process(rk::LT& s, Stream& is,
                                   State* st, const bool processClassId)
        {
            rk::LT* ps = &s;
            return readIntoPtr(ps, is, st, processClassId);
        }
    };
}

#endif // RK_RKIO_HH_
