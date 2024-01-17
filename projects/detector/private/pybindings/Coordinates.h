#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../../math/public/LeptonInjector/math/Vector3D.h"
#include "../../public/LeptonInjector/detector/DetectorModel.h"
#include "../../public/LeptonInjector/detector/Coordinates.h"

void register_Coordinates(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::math;
    using namespace LI::detector;

    // Register the NamedTypes

    class_<DetectorPosition, std::shared_ptr<DetectorPosition>>(m, "DetectorPosition")
        .def(init<Vector3D const &>())
        .def("get", (
                    Vector3D & (DetectorPosition::*)()
                    )(&DetectorPosition::get))
        ;

    class_<DetectorDirection, std::shared_ptr<DetectorDirection>>(m, "DetectorDirection")
        .def(init<Vector3D const &>())
        .def("get", (
                    Vector3D & (DetectorDirection::*)()
                    )(&DetectorDirection::get))
        ;

    class_<GeometryPosition, std::shared_ptr<GeometryPosition>>(m, "GeometryPosition")
        .def(init<Vector3D const &>())
        .def("get", (
                    Vector3D & (GeometryPosition::*)()
                    )(&GeometryPosition::get))
        ;

    class_<GeometryDirection, std::shared_ptr<GeometryDirection>>(m, "GeometryDirection")
        .def(init<Vector3D const &>())
        .def("get", (
                    Vector3D & (GeometryDirection::*)()
                    )(&GeometryDirection::get))
        ;

    // Allow implicit conversion from Vector3D to DetectorPosition and DetectorDirection
    // But only on the python side!
    // This works because the interface exposed to python does not include methods that use GeometryPosition and GeometryDirection
    implicitly_convertible<Vector3D, DetectorPosition>();
    implicitly_convertible<Vector3D, DetectorDirection>();
}
