#include <vector>
#include <string>
#include <map>
#include <utility>
#include <tuple>

#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Coordinates.h>
#include <LeptonInjector/Controller.h>
#include <LeptonInjector/Random.h>
#include <LeptonInjector/Constants.h>
#include <earthmodel-service/EarthModelCalculator.h>
#include <earthmodel-service/EarthModel.h>

// #include <converter/LeptonInjectionConfigurationConverter.h>
#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>


#include "container_conversions.h"

using namespace boost::python;

template<class T>
struct VecToList
{
	static PyObject* convert(const std::vector<T>& vec){
		boost::python::list* l = new boost::python::list();
		for(size_t i =0; i < vec.size(); i++)
			(*l).append(vec[i]);

		return l->ptr();
	}
};


template<typename T, typename U, typename P>
struct ThreeTupleToPyTuple
{
	static PyObject* convert(const std::tuple<T,U,P>& tup){
		boost::python::list* l = new boost::python::list();
		(*l).append(std::get<0>(tup));
		(*l).append(std::get<1>(tup));
		(*l).append(std::get<2>(tup));

		return l->ptr();
	}
};

template<class T>
struct DeqToList
{
	static PyObject* convert(const std::deque<T>& vec){
		boost::python::list* l = new boost::python::list();
		for(size_t i =0; i < vec.size(); i++)
			(*l).append(vec[i]);

		return l->ptr();
	}
};

void ListToVec(std::vector<unsigned int> &ret, boost::python::list l){
	for(int i=0;i<boost::python::len(l);i++)
		ret.push_back(boost::python::extract<unsigned int>(l[i]));
}

namespace earthmodel{
	struct LIEarthModelCalculator {
	};
}

BOOST_PYTHON_MODULE(EarthModelService){
	using namespace earthmodel;
    typedef std::vector<EarthSector> Sectors;

    class_<Sectors>("Sectors")
        .def(vector_indexing_suite<Sectors>());

    class_<Vector3D>("Vector3D", init<>())
        .def(init<const double, const double, const double>())
        .def(init<const Vector3D&>())
        .def("magnitude",&Vector3D::magnitude)
        .def("normalize",&Vector3D::normalize)
        .def("deflect",&Vector3D::deflect)
        .def("GetZ", &Vector3D::GetZ)
        .def("GetY", &Vector3D::GetY)
        .def("GetX", &Vector3D::GetX)
        .def("GetRadius",&Vector3D::GetRadius)
        .def("GetPhi",&Vector3D::GetPhi)
        .def("GetTheta",&Vector3D::GetTheta)
        .def(self + self)
        .def(self - self)
        .def(self += self)
        .def(self == self)
        .def(self * double())
        .def(double() * self)
        .def(self * self)
        .def("CalculateCartesianFromSpherical",&Vector3D::CalculateCartesianFromSpherical)
        .def("CalculateSphericalCoordinates",&Vector3D::CalculateSphericalCoordinates)
        .def("GetCartesianCoordinates",&Vector3D::GetCartesianCoordinates)
        .def("GetSphericalCoordinates",&Vector3D::GetSphericalCoordinates)
        .def("SetCartesianCoordinates",&Vector3D::SetCartesianCoordinates)
        .def("SetSphericalCoordinates",&Vector3D::SetSphericalCoordinates)
        ;

	class_<EarthSector>("EarthSector", init<>())
        .def_readwrite("name", &EarthSector::name)
        .def_readwrite("material_id", &EarthSector::material_id)
        .def_readwrite("level", &EarthSector::level)
        .def_readwrite("geo", &EarthSector::geo)
        .def_readwrite("density", &EarthSector::density)
		;

    double (EarthModel::*GetMassDensity_cached)(Geometry::IntersectionList const & intersections, Vector3D const & p0) const = &EarthModel::GetMassDensity;
    double (EarthModel::*GetMassDensity)(Vector3D const & p0) const = &EarthModel::GetMassDensity;

    double (EarthModel::*GetColumnDepthInCGS_cached)(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1) const = &EarthModel::GetColumnDepthInCGS;
    double (EarthModel::*GetColumnDepthInCGS)(Vector3D const & p0, Vector3D const & p1) const = &EarthModel::GetColumnDepthInCGS;

    double (EarthModel::*DistanceForColumnDepthToPoint_cached)(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double column_depth) const = &EarthModel::DistanceForColumnDepthToPoint;
    double (EarthModel::*DistanceForColumnDepthToPoint)(Vector3D const & end_point, Vector3D const & direction, double column_depth) const = &EarthModel::DistanceForColumnDepthToPoint;


	class_<EarthModel>("EarthModel", init<>())
        .def(init<const std::string&, const std::string&>())
        .def(init<const std::string &, const std::string &, const std::string &>())
        .def("LoadEarthModel",&EarthModel::LoadEarthModel)
        .def("LoadMaterialModel",&EarthModel::LoadMaterialModel)
        //.def("GetColumnDepthInCGS",&EarthModel::GetColumnDepthInCGS)
        .def("GetColumnDepthInCGS", GetMassDensity)
        .def("GetColumnDepthInCGS", GetMassDensity_cached)
        .def("DistanceForColumnDepthToPoint", DistanceForColumnDepthToPoint)
        .def("DistanceForColumnDepthToPoint", DistanceForColumnDepthToPoint_cached)
        .def("GetEarthCoordPosFromDetCoordPos",&EarthModel::GetEarthCoordPosFromDetCoordPos)
        .def("GetEarthCoordDirFromDetCoordDir",&EarthModel::GetEarthCoordDirFromDetCoordDir)
        .def("GetDetCoordPosFromEarthCoordPos",&EarthModel::GetDetCoordPosFromEarthCoordPos)
        .def("GetDetCoordDirFromEarthCoordDir",&EarthModel::GetDetCoordDirFromEarthCoordDir)
        .def("GetPath",&EarthModel::GetPath)
        .def("SetPath",&EarthModel::SetPath)
        .def("GetMaterials",&EarthModel::GetMaterials, return_value_policy<copy_const_reference>())
        .def("SetMaterials",&EarthModel::SetMaterials)
        .def("GetSectors",&EarthModel::GetSectors, return_value_policy<copy_const_reference>())
        .def("SetSectors",&EarthModel::SetSectors)
        .def("GetDetectorOrigin",&EarthModel::GetDetectorOrigin)
        .def("SetDetectorOrigin",&EarthModel::SetDetectorOrigin)
        .def("AddSector",&EarthModel::AddSector)
        .def("GetSector",&EarthModel::GetSector)
        .def("ClearSectors",&EarthModel::ClearSectors)
        ;

	{
		scope earthmodel = class_<LIEarthModelCalculator>("EarthModelCalculator");

		def("GetImpactParameter", &EarthModelCalculator::GetImpactParameter);
		def("GetIntersectionsWithSphere", &EarthModelCalculator::GetIntersectionsWithSphere);
		def("GetDistsToIntersectionsWithSphere", &EarthModelCalculator::GetDistsToIntersectionsWithSphere);
		def("GetLeptonRange", &EarthModelCalculator::GetLeptonRange);
		def("ColumnDepthCGStoMWE",&EarthModelCalculator::ColumnDepthCGStoMWE);
		def("MWEtoColumnDepthCGS",&EarthModelCalculator::MWEtoColumnDepthCGS);
	}

	using namespace scitbx::boost_python::container_conversions;
	from_python_sequence< std::vector<double>, variable_capacity_policy >();
	to_python_converter< std::vector<double, class std::allocator<double> >, VecToList<double> > ();

    //from_python_sequence< std::tuple<double,double,double>, variable_capacity_policy >();
	to_python_converter< std::tuple<double, double,double>, ThreeTupleToPyTuple<double,double,double> > ();

	from_python_sequence< std::vector<std::tuple<double,double,double>>, variable_capacity_policy >();
	to_python_converter< std::vector<std::tuple<double,double,double>, class std::allocator<std::tuple<double,double,double>> >, VecToList<std::tuple<double,double,double>> > ();

	from_python_sequence< std::vector<int>, variable_capacity_policy >();
	to_python_converter< std::vector<int, class std::allocator<int> >, VecToList<int> > ();

	from_python_sequence< std::vector<unsigned int>, variable_capacity_policy >();
	to_python_converter< std::vector<unsigned int, class std::allocator<unsigned int> >, VecToList<unsigned int> > ();

	from_python_sequence< std::vector<std::string>, variable_capacity_policy >();
	to_python_converter< std::vector<std::string, class std::allocator<std::string> >, VecToList<std::string> > ();
}
