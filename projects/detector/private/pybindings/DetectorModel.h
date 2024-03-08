#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/detector/DetectorModel.h"
#include "../../public/SIREN/detector/Coordinates.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"

void register_DetectorModel(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace SI::detector;

    class_<DetectorModel, std::shared_ptr<DetectorModel>>(m, "DetectorModel")
        .def(init<>())
        .def(init<std::string const &, std::string const &>())
        .def(init<std::string const &, std::string const &, std::string const &>())
        .def("LoadDetectorModel",&DetectorModel::LoadDetectorModel)
        .def("LoadMaterialModel",&DetectorModel::LoadMaterialModel)
        .def("GetMassDensity", (
                    double (DetectorModel::*)(SI::geometry::Geometry::IntersectionList const &, DetectorPosition const &) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(DetectorPosition const &) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetParticleDensity", (
                    double (DetectorModel::*)(SI::geometry::Geometry::IntersectionList const &, DetectorPosition const &, SI::dataclasses::ParticleType) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetParticleDensity", (
                    double (DetectorModel::*)(DetectorPosition const &, SI::dataclasses::ParticleType) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetInteractionDensity", (
                    double (DetectorModel::*)(SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDensity))
        .def("GetInteractionDensity", (
                    double (DetectorModel::*)(DetectorPosition const &,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDensity))
        .def("GetColumnDepthInCGS", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorPosition const &) const
                    )(&DetectorModel::GetColumnDepthInCGS))
        .def("GetColumnDepthInCGS", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorPosition const &) const
                    )(&DetectorModel::GetColumnDepthInCGS))
        .def("DistanceForColumnDepthFromPoint", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthFromPoint))
        .def("DistanceForColumnDepthFromPoint", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthFromPoint))
        .def("DistanceForColumnDepthToPoint", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthToPoint))
        .def("DistanceForColumnDepthToPoint", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double) const
                    )(&DetectorModel::DistanceForColumnDepthToPoint))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        std::set<SI::dataclasses::ParticleType>) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetMassDensity", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        std::set<SI::dataclasses::ParticleType>) const
                    )(&DetectorModel::GetMassDensity))
        .def("GetParticleDensity", (
                    std::vector<double> (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        std::set<SI::dataclasses::ParticleType>) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetParticleDensity", (
                    std::vector<double> (DetectorModel::*)(
                        DetectorPosition const &,
                        std::set<SI::dataclasses::ParticleType>) const
                    )(&DetectorModel::GetParticleDensity))
        .def("GetInteractionDepthInCGS", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorPosition const &,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDepthInCGS))
        .def("GetInteractionDepthInCGS", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorPosition const &,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::GetInteractionDepthInCGS))
        .def("DistanceForInteractionDepthFromPoint", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthFromPoint))
        .def("DistanceForInteractionDepthFromPoint", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthFromPoint))
        .def("DistanceForInteractionDepthToPoint", (
                    double (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthToPoint))
        .def("DistanceForInteractionDepthToPoint", (
                    double (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorDirection const &,
                        double,
                        std::vector<SI::dataclasses::ParticleType> const &,
                        std::vector<double> const &,
                        double const &) const
                    )(&DetectorModel::DistanceForInteractionDepthToPoint))
        .def("GetParticleColumnDepth", (
                    std::vector<double> (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &,
                        DetectorPosition const &,
                        std::vector<SI::dataclasses::ParticleType> const &) const
                    )(&DetectorModel::GetParticleColumnDepth))
        .def("GetContainingSector", (
                    DetectorSector (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const & p0) const
                    )(&DetectorModel::GetContainingSector))
        .def("GetContainingSector", (
                    DetectorSector (DetectorModel::*)(
                        DetectorPosition const & p0) const
                    )(&DetectorModel::GetContainingSector))
        .def_property("Path", &DetectorModel::GetPath, &DetectorModel::SetPath)
        .def_property("Materials", &DetectorModel::GetMaterials, &DetectorModel::SetMaterials)
        .def_property("Sectors", &DetectorModel::GetSectors, &DetectorModel::SetSectors)
        .def_property("DetectorOrigin", &DetectorModel::GetDetectorOrigin, &DetectorModel::SetDetectorOrigin)
        .def("AddSector", &DetectorModel::AddSector)
        .def("GetSector", &DetectorModel::GetSector)
        .def("ClearSectors", &DetectorModel::ClearSectors)
        .def("GetIntersections", (
                    SI::geometry::Geometry::IntersectionList (DetectorModel::*)(
                        DetectorPosition const &,
                        DetectorDirection const &) const
                    )(&DetectorModel::GetIntersections))
        .def_static("SortIntersections", (
                    void (*)(SI::geometry::Geometry::IntersectionList &)
                    )(&DetectorModel::SortIntersections))
        .def_static("SortIntersections", (
                    void (*)(std::vector<SI::geometry::Geometry::Intersection> &)
                    )(&DetectorModel::SortIntersections))
        .def_static("SectorLoop", &DetectorModel::SectorLoop)
        .def_static("GetOuterBoundsFromIntersections", (
                    SI::geometry::Geometry::IntersectionList (*)(
                        SI::geometry::Geometry::IntersectionList const &)
                    )(&DetectorModel::GetOuterBounds))
        .def("GetOuterBounds", [](SI::detector::DetectorModel const & detector, SI::geometry::Geometry::IntersectionList const & intersections){
                return detector.GetOuterBounds(intersections);
                })
        .def("GetOuterBounds", (
                SI::geometry::Geometry::IntersectionList (DetectorModel::*)(
                    DetectorPosition const &,
                    DetectorDirection const &) const
                )(&DetectorModel::GetOuterBounds))
        .def("GetAvailableTargets", (
                    std::set<SI::dataclasses::ParticleType> (DetectorModel::*)(
                        SI::geometry::Geometry::IntersectionList const &,
                        DetectorPosition const &) const
                    )(&DetectorModel::GetAvailableTargets))
        .def("GetAvailableTargets", (
                    std::set<SI::dataclasses::ParticleType> (DetectorModel::*)(
                        DetectorPosition const &) const
                    )(&DetectorModel::GetAvailableTargets))
        .def("GetTargetMass", &DetectorModel::GetTargetMass)
        .def("GetMaterials",&DetectorModel::GetMaterials)
        .def("GetDetectorOrigin",&DetectorModel::GetDetectorOrigin)
        .def("GeoPositionToDetPosition", (
                    DetectorPosition (DetectorModel::*)(GeometryPosition const &) const
                    )(&DetectorModel::ToDet))
        .def("GeoDirectionToDetDirection", (
                    DetectorDirection (DetectorModel::*)(GeometryDirection const &) const
                    )(&DetectorModel::ToDet))
        .def("DetPositionToGeoPosition", (
                    GeometryPosition (DetectorModel::*)(DetectorPosition const &) const
                    )(&DetectorModel::ToGeo))
        .def("DetDirectionToGeoDirection", (
                    GeometryDirection (DetectorModel::*)(DetectorDirection const &) const
                    )(&DetectorModel::ToGeo))
        .def("ParseFiducialVolume", (std::shared_ptr<SI::geometry::Geometry> (*)(std::string, std::string))(&DetectorModel::ParseFiducialVolume))
        .def("ParseFiducialVolume", (std::shared_ptr<SI::geometry::Geometry> (*)(std::string, SI::math::Vector3D, SI::math::Quaternion))(&DetectorModel::ParseFiducialVolume))
        ;
}
