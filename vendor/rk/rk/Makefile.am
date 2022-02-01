AM_CXXFLAGS = -O3

lib_LTLIBRARIES = librk.la

librk_la_SOURCES = rk.cc geom3.cc

INCLUDES = -I@top_srcdir@/

includedir = ${prefix}/include/rk

include_HEADERS = rk.hh geom3.hh geom3_Matrix3x3.icc geom3_Point3.icc \
    geom3_Rotation3.icc geom3_UnitVector3.icc geom3_Vector3.icc rk_Boost.icc \
    rk_P4.icc rk_Point4.icc LT.hh rk_LT.icc rkIO.hh
