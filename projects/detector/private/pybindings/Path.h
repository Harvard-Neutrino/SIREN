#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Coordinates.h"
#include "../../public/LeptonInjector/detector/Path.h"

void register_Path(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;
    using namespace LI::math;

    class_<Path, std::shared_ptr<Path>>(m, "Path")
        .def(init<std::shared_ptr<const DetectorModel>>())
        .def(init<std::shared_ptr<const DetectorModel>, DetectorPosition const &, DetectorPosition const &>())
        .def(init<std::shared_ptr<const DetectorModel>, DetectorPosition const &, DetectorDirection const &, double>())

        .def("HasDetectorModel", &Path::HasDetectorModel)
        .def("HasPoints", &Path::HasPoints)
        .def("HasIntersections", &Path::HasIntersections)
        .def("HasColumnDepth", &Path::HasColumnDepth)

        .def("GetDetectorModel", &Path::GetDetectorModel)
        .def("GetFirstPoint", [](Path & p)->Vector3D{return p.GetFirstPoint().get(); })
        .def("GetLastPoint", [](Path & p)->Vector3D{return p.GetLastPoint().get(); })
        .def("GetDirection", [](Path & p)->Vector3D{return p.GetDirection().get(); })
        .def("GetGeoFirstPoint", &Path::GetGeoFirstPoint)
        .def("GetGeoLastPoint", &Path::GetGeoLastPoint)
        .def("GetGeoDirection", &Path::GetGeoDirection)
        .def("GetDistance", &Path::GetDistance)
        .def("GetIntersections", &Path::GetIntersections)

        .def("SetDetectorModel", &Path::SetDetectorModel)
        .def("EnsureDetectorModel", &Path::EnsureDetectorModel)

        .def("SetPoints", overload_cast<
                DetectorPosition,
                DetectorPosition
                >(&Path::SetPoints))
        .def("SetPointsWithRay", overload_cast<
                DetectorPosition,
                DetectorDirection,
                double
                >(&Path::SetPointsWithRay))
        .def("EnsurePoints", &Path::EnsurePoints)

        .def("SetIntersections", &Path::SetIntersections)
        .def("ComputeIntersections", &Path::ComputeIntersections)
        .def("EnsureIntersections", &Path::EnsureIntersections)

        .def("ClipToOuterBounds", &Path::ClipToOuterBounds)

        .def("Flip", &Path::Flip)

        // Extend/Shrink By
        .def("ExtendFromEndByDistance", &Path::ExtendFromEndByDistance)
        .def("ExtendFromStartByDistance", &Path::ExtendFromStartByDistance)
        .def("ShrinkFromEndByDistance", &Path::ShrinkFromEndByDistance)
        .def("ShrinkFromStartByDistance", &Path::ShrinkFromStartByDistance)

        .def("ExtendFromEndByColumnDepth", &Path::ExtendFromEndByColumnDepth)
        .def("ExtendFromStartByColumnDepth", &Path::ExtendFromStartByColumnDepth)
        .def("ShrinkFromEndByColumnDepth", &Path::ShrinkFromEndByColumnDepth)
        .def("ShrinkFromStartByColumnDepth", &Path::ShrinkFromStartByColumnDepth)

        .def("ExtendFromEndByInteractionDepth", &Path::ExtendFromEndByInteractionDepth)
        .def("ExtendFromStartByInteractionDepth", &Path::ExtendFromStartByInteractionDepth)
        .def("ShrinkFromEndByInteractionDepth", &Path::ShrinkFromEndByInteractionDepth)
        .def("ShrinkFromStartByInteractionDepth", &Path::ShrinkFromStartByInteractionDepth)

        // Extend/Shrink To
        .def("ExtendFromEndToDistance", &Path::ExtendFromEndToDistance)
        .def("ExtendFromStartToDistance", &Path::ExtendFromStartToDistance)
        .def("ShrinkFromEndToDistance", &Path::ShrinkFromEndToDistance)
        .def("ShrinkFromStartToDistance", &Path::ShrinkFromStartToDistance)

        .def("ExtendFromEndToColumnDepth", &Path::ExtendFromEndToColumnDepth)
        .def("ExtendFromStartToColumnDepth", &Path::ExtendFromStartToColumnDepth)
        .def("ShrinkFromEndToColumnDepth", &Path::ShrinkFromEndToColumnDepth)
        .def("ShrinkFromStartToColumnDepth", &Path::ShrinkFromStartToColumnDepth)

        .def("ExtendFromEndToInteractionDepth", &Path::ExtendFromEndToInteractionDepth)
        .def("ExtendFromStartToInteractionDepth", &Path::ExtendFromStartToInteractionDepth)
        .def("ShrinkFromEndToInteractionDepth", &Path::ShrinkFromEndToInteractionDepth)
        .def("ShrinkFromStartToInteractionDepth", &Path::ShrinkFromStartToInteractionDepth)
        //

        // Get
        .def("GetColumnDepthInBounds", &Path::GetColumnDepthInBounds)
        .def("GetInteractionDepthInBounds", &Path::GetInteractionDepthInBounds)
        //

        // Get * From
        .def("GetColumnDepthFromStartInBounds", &Path::GetColumnDepthFromStartInBounds)
        .def("GetColumnDepthFromEndInBounds", &Path::GetColumnDepthFromEndInBounds)
        .def("GetColumnDepthFromStartAlongPath", &Path::GetColumnDepthFromStartAlongPath)
        .def("GetColumnDepthFromEndAlongPath", &Path::GetColumnDepthFromEndAlongPath)
        .def("GetColumnDepthFromStartInReverse", &Path::GetColumnDepthFromStartInReverse)
        .def("GetColumnDepthFromEndInReverse", &Path::GetColumnDepthFromEndInReverse)

        .def("GetInteractionDepthFromStartInBounds", &Path::GetInteractionDepthFromStartInBounds)
        .def("GetInteractionDepthFromEndInBounds", &Path::GetInteractionDepthFromEndInBounds)
        .def("GetInteractionDepthFromStartAlongPath", &Path::GetInteractionDepthFromStartAlongPath)
        .def("GetInteractionDepthFromEndAlongPath", &Path::GetInteractionDepthFromEndAlongPath)
        .def("GetInteractionDepthFromStartInReverse", &Path::GetInteractionDepthFromStartInReverse)
        .def("GetInteractionDepthFromEndInReverse", &Path::GetInteractionDepthFromEndInReverse)
        //

        // Get Distance From
        .def("GetDistanceFromStartInBounds",
                overload_cast<double>
                (&Path::GetDistanceFromStartInBounds))
        .def("GetDistanceFromEndInBounds",
                overload_cast<double>
                (&Path::GetDistanceFromEndInBounds))
        .def("GetDistanceFromStartAlongPath",
                overload_cast<double>
                (&Path::GetDistanceFromStartAlongPath))
        .def("GetDistanceFromEndAlongPath",
                overload_cast<double>
                (&Path::GetDistanceFromEndAlongPath))
        .def("GetDistanceFromStartInReverse",
                overload_cast<double>
                (&Path::GetDistanceFromStartInReverse))
        .def("GetDistanceFromEndInReverse",
                overload_cast<double>
                (&Path::GetDistanceFromEndInReverse))

        .def("GetDistanceFromStartInBounds", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromStartInBounds))
        .def("GetDistanceFromEndInBounds", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromEndInBounds))
        .def("GetDistanceFromStartAlongPath", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromStartAlongPath))
        .def("GetDistanceFromEndAlongPath", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromEndAlongPath))
        .def("GetDistanceFromStartInReverse", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromStartInReverse))
        .def("GetDistanceFromEndInReverse", overload_cast<
                double,
                std::vector<LI::dataclasses::Particle::ParticleType> const &,
                std::vector<double> const &,
                double const &
                >(&Path::GetDistanceFromEndInReverse))
        //

        .def("IsWithinBounds", overload_cast<DetectorPosition>(&Path::IsWithinBounds))
        .def("GetDistanceFromStartInBounds", overload_cast<DetectorPosition>(&Path::GetDistanceFromStartInBounds))
        ;
}
