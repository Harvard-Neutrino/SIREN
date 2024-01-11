#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"

void register_DetectorModel(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<DetectorModel, std::shared_ptr<DetectorModel>>(m, "DetectorModel")
        .def(init<>())
        .def(init<std::string const &, std::string const &>())
        .def(init<std::string const &, std::string const &, std::string const &>())
        .def("LoadDetectorModel",&DetectorModel::LoadDetectorModel)
        .def("LoadMaterialModel",&DetectorModel::LoadMaterialModel)
        .def("GetMassDensity", (
                    double (DetectorModel::*)(LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(LI::math::Vector3D const &) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetParticleDensity", (
                    double (DetectorModel::*)(LI::geometry::Geometry::IntersectionList const &, LI::math::Vector3D const &, LI::dataclasses::Particle::ParticleType) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetParticleDensity", (
                    double (DetectorModel::*)(LI::math::Vector3D const &, LI::dataclasses::Particle::ParticleType) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetInteractionDensity", (
                    double (DetectorModel::*)(LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDensity))
        .def("GetInteractionDensity", (
                    double (DetectorModel::*)(LI::math::Vector3D const &,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDensity))
        .def("GetColumnDepthInCGS", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetColumnDepthInCGS))
        .def("GetColumnDepthInCGS", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetColumnDepthInCGS))
        .def("DistanceForColumnDepthFromPoint", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthFromPoint))
        .def("DistanceForColumnDepthFromPoint", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthFromPoint))
        .def("DistanceForColumnDepthToPoint", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthToPoint))
        .def("DistanceForColumnDepthToPoint", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthToPoint))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        std::set<LI::dataclasses::Particle::ParticleType>) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        std::set<LI::dataclasses::Particle::ParticleType>) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetParticleDensity", (
                    std::vector<double> (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        std::set<LI::dataclasses::Particle::ParticleType>) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetParticleDensity", (
                    std::vector<double> (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        std::set<LI::dataclasses::Particle::ParticleType>) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetInteractionDepthInCGS", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDepthInCGS))
        .def("GetInteractionDepthInCGS", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDepthInCGS))
        .def("DistanceForInteractionDepthFromPoint", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthFromPoint))
        .def("DistanceForInteractionDepthFromPoint", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthFromPoint))
        .def("DistanceForInteractionDepthToPoint", (
                    double (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthToPoint))
        .def("DistanceForInteractionDepthToPoint", (
                    double (DetectorModel::*)(
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        double,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthToPoint))
        .def("GetParticleColumnDepth", (
                    std::vector<double> (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const &,
                        LI::math::Vector3D const &,
                        std::vector<LI::dataclasses::Particle::ParticleType> const &) const
                    )(&DetectorModel::GetParticleColumnDepth))
        .def("GetContainingSector", (
                    DetectorSector (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        LI::math::Vector3D const & p0) const
                    )(&DetectorModel::GetContainingSector))
        .def("GetContainingSector", (
                    DetectorSector (DetectorModel::*)(
                        LI::math::Vector3D const & p0) const
                    )(&DetectorModel::GetContainingSector))
        .def("GetEarthCoordPosFromDetCoordPos", (
                    LI::math::Vector3D (DetectorModel::*)(
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetEarthCoordPosFromDetCoordPos))
        .def("GetEarthCoordDirFromDetCoordDir", (
                    LI::math::Vector3D (DetectorModel::*)(
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetEarthCoordDirFromDetCoordDir))
        .def("GetDetCoordPosFromEarthCoordPos", (
                    LI::math::Vector3D (DetectorModel::*)(
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetDetCoordPosFromEarthCoordPos))
        .def("GetDetCoordDirFromEarthCoordDir", (
                    LI::math::Vector3D (DetectorModel::*)(
                        LI::math::Vector3D const &) const
                    )(&DetectorModel::GetDetCoordDirFromEarthCoordDir))
        .def_property("Path", &DetectorModel::GetPath, &DetectorModel::SetPath)
        .def_property("Materials", &DetectorModel::GetMaterials, &DetectorModel::SetMaterials)
        .def_property("Sectors", &DetectorModel::GetSectors, &DetectorModel::SetSectors)
        .def_property("DetectorOrigin", &DetectorModel::GetDetectorOrigin, &DetectorModel::SetDetectorOrigin)
        .def("AddSector", &DetectorModel::AddSector)
        .def("GetSector", &DetectorModel::GetSector)
        .def("ClearSectors", &DetectorModel::ClearSectors)
        .def("GetIntersections", &DetectorModel::GetIntersections)
        .def_static("SortIntersections", (
                    void (*)(LI::geometry::Geometry::IntersectionList &)
                    )(&DetectorModel::SortIntersections))
        .def_static("SortIntersections", (
                    void (*)(std::vector<LI::geometry::Geometry::Intersection> &)
                    )(&DetectorModel::SortIntersections))
        .def_static("SectorLoop", &DetectorModel::SectorLoop)
        .def_static("GetOuterBoundsFromIntersections", (
                    LI::geometry::Geometry::IntersectionList (*)(
                        LI::geometry::Geometry::IntersectionList const &)
                    )(&DetectorModel::GetOuterBounds))
        .def("GetOuterBounds", [](LI::detector::DetectorModel const & detector, LI::geometry::Geometry::IntersectionList const & intersections){
                return detector.GetOuterBounds(intersections);
                })
        .def("GetOuterBounds", (
                LI::geometry::Geometry::IntersectionList (DetectorModel::*)(
                    LI::math::Vector3D const &,
                    LI::math::Vector3D const &) const
                )(&DetectorModel::GetOuterBounds))
        .def("GetAvailableTargets", (
                    std::set<LI::dataclasses::Particle::ParticleType> (DetectorModel::*)(
                        LI::geometry::Geometry::IntersectionList const &,
                        std::array<double,3> const &) const
                    )(&DetectorModel::GetAvailableTargets))
        .def("GetAvailableTargets", (
                    std::set<LI::dataclasses::Particle::ParticleType> (DetectorModel::*)(
                        std::array<double,3> const &) const
                    )(&DetectorModel::GetAvailableTargets))
        .def("GetTargetMass", &DetectorModel::GetTargetMass)
        .def("GetMaterials",&DetectorModel::GetMaterials)
        .def("GetDetectorOrigin",&DetectorModel::GetDetectorOrigin);

}
