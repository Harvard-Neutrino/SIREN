
#include <vector>

#include "../../public/LeptonInjector/detector/EarthModel.h"
#include "../../public/LeptonInjector/detector/Path.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>



using namespace pybind11;

PYBIND11_MODULE(detector,m) {
  using namespace LI::detector;

  class_<EarthModel, std::shared_ptr<EarthModel>>(m, "EarthModel")
    .def(init<>())
    .def(init<std::string const &, std::string const &>())
    .def(init<std::string const &, std::string const &, std::string const &>())
    .def("LoadEarthModel",&EarthModel::LoadEarthModel)
    .def("LoadMaterialModel",&EarthModel::LoadMaterialModel)
    .def("GetSectors",&EarthModel::GetSectors)
    // .def("GetColumnDepthInCGS",overload_cast<LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::math::Vector3D const &>(&EarthModel::GetColumnDepthInCGS))
    // .def("GetColumnDepthInCGS",overload_cast<LI::math::Vector3D const &, LI::math::Vector3D const &>(&EarthModel::GetColumnDepthInCGS))
    // .def("GetInteractionDepthInCGS",overload_cast<LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::math::Vector3D const &, std::vector< LI::dataclasses::Particle::ParticleType > const &, std::vector< double > const &, double const &>(&EarthModel::GetInteractionDepthInCGS))
    // .def("GetInteractionDepthInCGS",overload_cast<LI::math::Vector3D const &, LI::math::Vector3D const &, std::vector< LI::dataclasses::Particle::ParticleType > const &, std::vector< double > const &, double const &>(&EarthModel::GetInteractionDepthInCGS))
    // .def("DistanceForColumnDepthFromPoint",&EarthModel::DistanceForColumnDepthFromPoint)
    // .def("DistanceForColumnDepthToPoint",&EarthModel::DistanceForColumnDepthToPoint)
    // .def("DistanceForInteractionDepthFromPoint",&EarthModel::DistanceForInteractionDepthFromPoint)
    // .def("DistanceForInteractionDepthToPoint",&EarthModel::DistanceForInteractionDepthToPoint) 
    // .def("DistanceForColumnDepthFromPoint",overload_cast<LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::math::Vector3D const &, double>(&EarthModel::DistanceForColumnDepthFromPoint))
    // .def("DistanceForColumnDepthFromPoint",overload_cast<LI::math::Vector3D const &, LI::math::Vector3D const &, double>(&EarthModel::DistanceForColumnDepthFromPoint))
      
    // .def("DistanceForColumnDepthFromPoint", static_cast<double (EarthModel::*)(LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::math::Vector3D const &, double)>(&EarthModel::DistanceForColumnDepthFromPoint))
    // .def("DistanceForColumnDepthFromPoint", static_cast<double (EarthModel::*)(LI::math::Vector3D const &, LI::math::Vector3D const &, double)>(&EarthModel::DistanceForColumnDepthFromPoint))
    
    .def("GetMaterials",&EarthModel::GetMaterials)
    .def("GetDetectorOrigin",&EarthModel::GetDetectorOrigin);

  class_<EarthSector, std::shared_ptr<EarthSector>>(m, "EarthSector")
    .def_readwrite("name",&EarthSector::name)
    .def_readwrite("material_id",&EarthSector::material_id)
    .def_readwrite("geo",&EarthSector::geo);
    
  class_<Path, std::shared_ptr<Path>>(m, "Path")
    .def(init<std::shared_ptr<const EarthModel>>())
    .def(init<std::shared_ptr<const EarthModel>, LI::math::Vector3D const &, LI::math::Vector3D const &>())
    .def(init<std::shared_ptr<const EarthModel>, LI::math::Vector3D const &, LI::math::Vector3D const &, double>())
      
    .def("HasEarthModel",&Path::HasEarthModel)
    .def("HasPoints",&Path::HasPoints)
    .def("HasIntersections",&Path::HasIntersections)
    .def("HasColumnDepth",&Path::HasColumnDepth)
    .def("GetEarthModel",&Path::GetEarthModel)
    .def("HasEarthModel",&Path::HasEarthModel)
    .def("GetFirstPoint",&Path::GetFirstPoint)
    .def("GetLastPoint",&Path::GetLastPoint)
    .def("GetDirection",&Path::GetDirection)
    .def("GetDistance",&Path::GetDistance)
      
    .def("ClipToOuterBounds",&Path::ClipToOuterBounds)
      
    .def("ExtendFromEndByDistance",&Path::ExtendFromEndByDistance)
    .def("ExtendFromStartByDistance",&Path::ExtendFromStartByDistance)
    .def("ShrinkFromEndByDistance",&Path::ShrinkFromEndByDistance)
    .def("ShrinkFromStartByDistance",&Path::ShrinkFromStartByDistance)
    .def("ExtendFromEndByColumnDepth",&Path::ExtendFromEndByColumnDepth)
    .def("ExtendFromStartByColumnDepth",&Path::ExtendFromStartByColumnDepth)
    .def("ShrinkFromEndByColumnDepth",&Path::ShrinkFromEndByColumnDepth)
    .def("ShrinkFromStartByColumnDepth",&Path::ShrinkFromStartByColumnDepth)
    .def("ExtendFromEndByInteractionDepth",&Path::ExtendFromEndByInteractionDepth)
    .def("ExtendFromStartByInteractionDepth",&Path::ExtendFromStartByInteractionDepth)
    .def("ShrinkFromEndByInteractionDepth",&Path::ShrinkFromEndByInteractionDepth)
    .def("ShrinkFromStartByInteractionDepth",&Path::ShrinkFromStartByInteractionDepth)
    .def("ExtendFromEndToDistance",&Path::ExtendFromEndToDistance)
    .def("ExtendFromStartToDistance",&Path::ExtendFromStartToDistance)
    .def("ShrinkFromEndToDistance",&Path::ShrinkFromEndToDistance)
    .def("ShrinkFromStartToDistance",&Path::ShrinkFromStartToDistance)
    .def("ExtendFromEndToColumnDepth",&Path::ExtendFromEndToColumnDepth)
    .def("ExtendFromStartToColumnDepth",&Path::ExtendFromStartToColumnDepth)
    .def("ShrinkFromEndToColumnDepth",&Path::ShrinkFromEndToColumnDepth)
    .def("ShrinkFromStartToColumnDepth",&Path::ShrinkFromStartToColumnDepth)
    .def("ExtendFromEndToInteractionDepth",&Path::ExtendFromEndToInteractionDepth)
    .def("ExtendFromStartToInteractionDepth",&Path::ExtendFromStartToInteractionDepth)
    .def("ShrinkFromEndToInteractionDepth",&Path::ShrinkFromEndToInteractionDepth)
    .def("ShrinkFromStartToInteractionDepth",&Path::ShrinkFromStartToInteractionDepth)
    .def("GetColumnDepthInBounds",&Path::GetColumnDepthInBounds)
    .def("GetInteractionDepthInBounds",&Path::GetInteractionDepthInBounds)
    
    .def("GetDistanceFromStartAlongPath",overload_cast<double, std::vector< LI::dataclasses::Particle::ParticleType > const &, std::vector< double > const &, double const &>(&Path::GetDistanceFromStartAlongPath))
    .def("GetDistanceFromStartAlongPath",overload_cast<double>(&Path::GetDistanceFromStartAlongPath))
    
    .def("GetDistanceFromStartInReverse",overload_cast<double, std::vector< LI::dataclasses::Particle::ParticleType > const &, std::vector< double > const &, double const &>(&Path::GetDistanceFromStartInReverse))
    .def("GetDistanceFromStartInReverse",overload_cast<double>(&Path::GetDistanceFromStartInReverse));

}
