#include "./DarkNewsCrossSection.h"
#include "../../public/LeptonInjector/interactions/DarkNewsCrossSection.h"
#include "./DarkNewsDecay.h"
#include "../../public/LeptonInjector/interactions/DarkNewsDecay.h"

#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

namespace py = pybind11;

std::string pyDarkNewsCrossSection_dumper(LI::interactions::pyDarkNewsCrossSection & object) {
    pybind11::object obj;
    if(object.self) {
        obj = object.self;
    } else {
        auto *tinfo = pybind11::detail::get_type_info(typeid(LI::interactions::DarkNewsCrossSection));
        pybind11::handle self_handle = get_object_handle(static_cast<const LI::interactions::DarkNewsCrossSection *>(&object), tinfo);
        obj = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
    }
    py::module pkl = py::module::import("pickle");
    py::bytes bytes = pkl.attr("dumps")(obj);
    std::string str = (std::string)(bytes.attr("hex")().cast<std::string>());
    return str;
}

LI::interactions::pyDarkNewsCrossSection pyDarkNewsCrossSection_loader(std::string & state) {
    LI::interactions::pyDarkNewsCrossSection object;
    py::module pkl = py::module::import("pickle");

    py::object fromhex = py::globals()["__builtins__"].attr("bytes").attr("fromhex");
    py::object bytes = fromhex(state);

    pkl.attr("loads")(bytes);
    object.self = pkl.attr("loads")(bytes);
    return object;
}

std::string pyDarkNewsDecay_dumper(LI::interactions::pyDarkNewsDecay & object) {
    pybind11::object obj;
    if(object.self) {
        obj = object.self;
    } else {
        auto *tinfo = pybind11::detail::get_type_info(typeid(LI::interactions::DarkNewsDecay));
        pybind11::handle self_handle = get_object_handle(static_cast<const LI::interactions::DarkNewsDecay *>(&object), tinfo);
        obj = pybind11::reinterpret_borrow<pybind11::object>(self_handle);
    }
    py::module pkl = py::module::import("pickle");
    py::bytes bytes = pkl.attr("dumps")(obj);
    std::string str = (std::string)(bytes.attr("hex")().cast<std::string>());
    return str;
}

LI::interactions::pyDarkNewsDecay pyDarkNewsDecay_loader(std::string & state) {
    LI::interactions::pyDarkNewsDecay object;
    py::module pkl = py::module::import("pickle");

    py::object fromhex = py::globals()["__builtins__"].attr("bytes").attr("fromhex");
    py::object bytes = fromhex(state);

    pkl.attr("loads")(bytes);
    object.self = pkl.attr("loads")(bytes);
    return object;
}

PYBIND11_MODULE(pyDarkNewsSerializer, m) {
    m.def("save_xsec", &pyDarkNewsCrossSection_dumper);
    m.def("load_xsec", &pyDarkNewsCrossSection_loader);
    m.def("save_decay", &pyDarkNewsDecay_dumper);
    m.def("load_decay", &pyDarkNewsDecay_loader);
}
