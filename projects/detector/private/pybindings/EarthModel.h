#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

void register_EarthModel(pybind11::module_ & m) {
    using namespace pybind11;
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
        .def_static("GetOuterBoundsFromIntersections", (
                    LI::geometry::Geometry::IntersectionList (*)(
                        LI::geometry::Geometry::IntersectionList const &)
                    )(&EarthModel::GetOuterBounds))
        .def("GetOuterBounds", [](LI::detector::EarthModel const & earth, LI::geometry::Geometry::IntersectionList const & intersections){
                return earth.GetOuterBounds(intersections);
                })
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
        .def("GetMaterials",&EarthModel::GetMaterials)
        .def("GetDetectorOrigin",&EarthModel::GetDetectorOrigin);

}
