
#include <vector>

#include "../../public/LeptonInjector/injection/Process.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(Process,m) {
  using namespace LI::injection;

  class_<Process>(m, "Process")
    .def_readwrite("primary_type",&Process::primary_type)
    .def_readwrite("cross_sections",&Process::cross_sections);

  class_<InjectionProcess, Process>(m, "InjectionProcess")
    .def_readwrite("injection_distributions",&InjectionProcess::injection_distributions);
  
  class_<PhysicalProcess, Process>(m, "PhysicalProcess")
    .def_readwrite("physical_distributions",&PhysicalProcess::physical_distributions);
}
