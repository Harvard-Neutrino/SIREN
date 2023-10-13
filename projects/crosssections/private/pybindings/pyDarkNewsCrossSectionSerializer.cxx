#include "./DarkNewsCrossSection.h"
#include "../../public/LeptonInjector/crosssections/DarkNewsCrossSection.h"

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

namespace py = pybind11;

std::string pyDarkNewsCrossSection_dumper(LI::crosssections::pyDarkNewsCrossSection & object) {
    pybind11::object obj;
    if(object.self) {
        obj = object.self;
    } else {
        auto *tinfo = pybind11::detail::get_type_info(typeid(LI::crosssections::DarkNewsCrossSection));
        pybind11::handle self_handle = get_object_handle(static_cast<const LI::crosssections::DarkNewsCrossSection *>(&object), tinfo);
        obj = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
    }
    py::module pkl = py::module::import("pickle");
    py::bytes bytes = pkl.attr("dumps")(obj);
    std::string str = (std::string)(bytes.attr("hex")().cast<std::string>());
    return str;
}

LI::crosssections::pyDarkNewsCrossSection pyDarkNewsCrossSection_loader(std::string & state) {
    LI::crosssections::pyDarkNewsCrossSection object;
    py::module pkl = py::module::import("pickle");

    py::object fromhex = py::globals()["__builtins__"].attr("bytes").attr("fromhex");
    py::object bytes = fromhex(state);

    pkl.attr("loads")(bytes);
    object.self = pkl.attr("loads")(bytes);
    return object;
}

PYBIND11_MODULE(pyDarkNewsCrossSectionSerializer, m) {
    m.def("save", &pyDarkNewsCrossSection_dumper);
    m.def("load", &pyDarkNewsCrossSection_loader);
}
