
#include <vector>

#include "../../public/SIREN/utilities/Random.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(utilities,m) {
  using namespace siren::utilities;

  class_<SIREN_random, std::shared_ptr<SIREN_random>>(m, "SIREN_random")
    .def(init<>())
    .def(init<unsigned int>())
    .def("Uniform",&SIREN_random::Uniform)
    .def("set_seed",&SIREN_random::set_seed);
}
