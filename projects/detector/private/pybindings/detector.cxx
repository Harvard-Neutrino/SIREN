
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
    .def("GetMassDensity", (
                double (EarthModel::*)(LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &) const
                )(&EarthModel::GetMassDensity))
    .def("GetMassDensity", (
                double (EarthModel::*)(LI::math::Vector3D const &) const
                )(&EarthModel::GetMassDensity))
    .def("GetParticleDensity", (
                double (EarthModel::*)(LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::dataclasses::Particle::ParticleType) const
                )(&EarthModel::GetParticleDensity))
    .def("GetParticleDensity", (
                double (EarthModel::*)(LI::math::Vector3D const &, LI::dataclasses::Particle::ParticleType) const
                )(&EarthModel::GetParticleDensity))
    .def("GetInteractionDensity", (
                double (EarthModel::*)(LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::GetInteractionDensity))
    .def("GetInteractionDensity", (
                double (EarthModel::*)(LI::math::Vector3D const &,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::GetInteractionDensity))
    .def("GetColumnDepthInCGS", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetColumnDepthInCGS))
    .def("GetColumnDepthInCGS", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetColumnDepthInCGS))
    .def("DistanceForColumnDepthFromPoint", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double) const
                )(&EarthModel::DistanceForColumnDepthFromPoint))
    .def("DistanceForColumnDepthFromPoint", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double) const
                )(&EarthModel::DistanceForColumnDepthFromPoint))
    .def("DistanceForColumnDepthToPoint", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double) const
                )(&EarthModel::DistanceForColumnDepthToPoint))
    .def("DistanceForColumnDepthToPoint", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double) const
                )(&EarthModel::DistanceForColumnDepthToPoint))
    .def("GetMassDensity", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    std::set<LI::dataclasses::Particle::ParticleType>) const
                )(&EarthModel::GetMassDensity))
    .def("GetMassDensity", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    std::set<LI::dataclasses::Particle::ParticleType>) const
                )(&EarthModel::GetMassDensity))
    .def("GetParticleDensity", (
                std::vector<double> (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    std::set<LI::dataclasses::Particle::ParticleType>) const
                )(&EarthModel::GetParticleDensity))
    .def("GetParticleDensity", (
                std::vector<double> (EarthModel::*)(
                    LI::math::Vector3D const &,
                    std::set<LI::dataclasses::Particle::ParticleType>) const
                )(&EarthModel::GetParticleDensity))
    .def("GetInteractionDepthInCGS", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::GetInteractionDepthInCGS))
    .def("GetInteractionDepthInCGS", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::GetInteractionDepthInCGS))
    .def("DistanceForInteractionDepthFromPoint", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::DistanceForInteractionDepthFromPoint))
    .def("DistanceForInteractionDepthFromPoint", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::DistanceForInteractionDepthFromPoint))
    .def("DistanceForInteractionDepthToPoint", (
                double (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::DistanceForInteractionDepthToPoint))
    .def("DistanceForInteractionDepthToPoint", (
                double (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    double,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &,
                    std::vector<double> const &,
                    double const &) const
                )(&EarthModel::DistanceForInteractionDepthToPoint))
    .def("GetParticleColumnDepth", (
                std::vector<double> (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &,
                    std::vector<LI::dataclasses::Particle::ParticleType> const &) const
                )(&EarthModel::GetParticleColumnDepth))
    .def("GetContainingSector", (
                EarthSector (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    LI::math::Vector3D const & p0) const
                )(&EarthModel::GetContainingSector))
    .def("GetContainingSector", (
                EarthSector (EarthModel::*)(
                    LI::math::Vector3D const & p0) const
                )(&EarthModel::GetContainingSector))
    .def("GetEarthCoordPosFromDetCoordPos", (
                LI::math::Vector3D (EarthModel::*)(
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetEarthCoordPosFromDetCoordPos))
    .def("GetEarthCoordDirFromDetCoordDir", (
                LI::math::Vector3D (EarthModel::*)(
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetEarthCoordDirFromDetCoordDir))
    .def("GetDetCoordPosFromEarthCoordPos", (
                LI::math::Vector3D (EarthModel::*)(
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetDetCoordPosFromEarthCoordPos))
    .def("GetDetCoordDirFromEarthCoordDir", (
                LI::math::Vector3D (EarthModel::*)(
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetDetCoordDirFromEarthCoordDir))
    .def_property("Path", &EarthModel::GetPath, &EarthModel::SetPath)
    .def_property("Materials", &EarthModel::GetMaterials, &EarthModel::SetMaterials)
    .def_property("Sectors", &EarthModel::GetSectors, &EarthModel::SetSectors)
    .def_property("DetectorOrigin", &EarthModel::GetDetectorOrigin, &EarthModel::SetDetectorOrigin)
    .def("AddSector", &EarthModel::AddSector)
    .def("GetSector", &EarthModel::GetSector)
    .def("ClearSectors", &EarthModel::ClearSectors)
    .def("GetIntersections", &EarthModel::GetIntersections)
    .def_static("SortIntersections", (
                void (*)(LI::geometry::Geometry::IntersectionList &)
                )(&EarthModel::SortIntersections))
    .def_static("SortIntersections", (
                void (*)(std::vector<LI::geometry::Geometry::Intersection> &)
                )(&EarthModel::SortIntersections))
    .def_static("SectorLoop", &EarthModel::SectorLoop)
    .def_static("GetOuterBounds", (
                LI::geometry::Geometry::IntersectionList (*)(
                    LI::geometry::Geometry::IntersectionList const &)
                )(&EarthModel::GetOuterBounds))
    .def("GetOuterBounds", (
                LI::geometry::Geometry::IntersectionList (EarthModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &) const
                )(&EarthModel::GetOuterBounds))
    .def("GetAvailableTargets", (
                std::set<LI::dataclasses::Particle::ParticleType> (EarthModel::*)(
                    LI::geometry::Geometry::IntersectionList const &,
                    std::array<double,3> const &) const
                )(&EarthModel::GetAvailableTargets))
    .def("GetAvailableTargets", (
                std::set<LI::dataclasses::Particle::ParticleType> (EarthModel::*)(
                    std::array<double,3> const &) const
                )(&EarthModel::GetAvailableTargets))
    .def("GetTargetMass", &EarthModel::GetTargetMass)
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
    .def_readwrite("level", &EarthSector::level)
    .def_readwrite("geo",&EarthSector::geo)
    .def_readwrite("density",&EarthSector::density);

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
