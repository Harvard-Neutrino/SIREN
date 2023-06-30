#include <memory>
#include <sstream>
#include <iostream>

#include <pybind11/pybind11.h>

#include "../../public/LeptonInjector/detector/EarthModel.h"

std::string to_str(LI::detector::EarthSector const & sector) {
    std::stringstream ss;
    sector.Print(ss);
    return ss.str();
}

void register_EarthSector(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::detector;

    class_<EarthSector, std::shared_ptr<EarthSector>>(m, "EarthSector")
        .def(init<>())
        .def("__str__", &to_str)
        .def_readwrite("name",&EarthSector::name)
        .def_readwrite("material_id",&EarthSector::material_id)
        .def_readwrite("level", &EarthSector::level)
        .def_readwrite("geo",&EarthSector::geo)
        .def_readwrite("density",&EarthSector::density);
}
