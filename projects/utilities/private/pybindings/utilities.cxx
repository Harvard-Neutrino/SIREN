
#include <vector>

#include "../../public/LeptonInjector/utilities/Random.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(utilities,m) {
  using namespace LI::utilities;

  class_<LI_random, std::shared_ptr<LI_random>>(m, "LI_random")
    .def(init<>())
    .def(init<unsigned int>())
    .def("Uniform",&LI_random::Uniform)
    .def("set_seed",&LI_random::set_seed);
}
