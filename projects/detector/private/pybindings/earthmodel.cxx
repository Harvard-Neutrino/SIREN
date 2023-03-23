
#include <vector>

#include "../../public/LeptonInjector/detector/EarthModel.h"

#include <pybind11/pybind11.h>



using namespace pybind11;

PYBIND11_MODULE(EarthModel,m) {
  using namespace LI::detector;

  class_<EarthModel>(m, "EarthModel")
    .def(init<>())
    .def(init<std::string const &, std::string const &>())
    .def(init<std::string const &, std::string const &, std::string const &>())
    .def("LoadEarthModel",&EarthModel::LoadEarthModel)
    .def("LoadMaterialModel",&EarthModel::LoadMaterialModel)
    .def("GetSectors",&EarthModel::GetSectors);
  
  class_<EarthSector>(m, "EarthSector")
    .def_readwrite("name",&EarthSector::name)
    .def_readwrite("material_id",&EarthSector::material_id)
    .def_readwrite("geo",&EarthSector::geo);

}
