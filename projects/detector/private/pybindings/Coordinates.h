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

    // Register only `DetectorPosition` and `DetectorDirection`
    // The python user does not need to know about or manipulate the Geometry coordinates

    class_<DetectorPosition, std::shared_ptr<DetectorPosition>>(m, "DetectorPosition")
        .def(init<Vector3D const &>())
        ;

    class_<DetectorDirection, std::shared_ptr<DetectorDirection>>(m, "DetectorDirection")
        .def(init<Vector3D const &>())
        ;

    // Allow implicit conversion from Vector3D to DetectorPosition and DetectorDirection
    // But only on the python side!
    implicitly_convertible<Vector3D, DetectorPosition>();
    implicitly_convertible<Vector3D, DetectorDirection>();
}
