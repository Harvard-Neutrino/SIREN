#include <memory>

#include <pybind11/pybind11.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"

void register_EarthSector(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<EarthSector, std::shared_ptr<EarthSector>>(m, "EarthSector")
        .def(init<>())
        .def_readwrite("name",&EarthSector::name)
        .def_readwrite("material_id",&EarthSector::material_id)
        .def_readwrite("level", &EarthSector::level)
        .def_readwrite("geo",&EarthSector::geo)
        .def_readwrite("density",&EarthSector::density);

}
