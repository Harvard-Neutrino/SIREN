#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Path.h"

void register_Path(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<Path, std::shared_ptr<Path>>(m, "Path")
        .def(init<std::shared_ptr<const DetectorModel>>())
        .def(init<std::shared_ptr<const DetectorModel>, LI::math::Vector3D const &, LI::math::Vector3D const &>())
        .def(init<std::shared_ptr<const DetectorModel>, LI::math::Vector3D const &, LI::math::Vector3D const &, double>())

        .def("HasDetectorModel",&Path::HasDetectorModel)
        .def("HasPoints",&Path::HasPoints)
        .def("HasIntersections",&Path::HasIntersections)
        .def("HasColumnDepth",&Path::HasColumnDepth)
        .def("GetDetectorModel",&Path::GetDetectorModel)
        .def("HasDetectorModel",&Path::HasDetectorModel)
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
