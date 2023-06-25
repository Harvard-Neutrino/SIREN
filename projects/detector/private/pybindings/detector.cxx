#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "./EarthModel.h"
#include "./EarthSector.h"
#include "./Path.h"
#include "./DensityDistribution.h"

using namespace pybind11;

PYBIND11_MODULE(detector,m) {
    register_EarthModel(m);
    register_EarthSector(m);
    register_Path(m);
    register_DensityDistribution(m);
}
