#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "./DetectorModel.h"
#include "./DetectorSector.h"
#include "./Path.h"
#include "./DensityDistribution.h"
#include "./Distribution1D.h"
#include "./ConstantDistribution1D.h"
#include "./ExponentialDistribution1D.h"
#include "./PolynomialDistribution1D.h"
#include "./Axis1D.h"
#include "./RadialAxis1D.h"
#include "./CartesianAxis1D.h"
#include "./CartesianAxisExponentialDensityDistribution.h"
#include "./CartesianAxisPolynomialDensityDistribution.h"
#include "./ConstantDensityDistribution.h"
#include "./MaterialModel.h"

using namespace pybind11;

PYBIND11_MODULE(detector,m) {
    register_DetectorModel(m);
    register_DetectorSector(m);
    register_Path(m);
    register_DensityDistribution(m);
    register_Distribution1D(m);
    register_ConstantDistribution1D(m);
    register_ExponentialDistribution1D(m);
    register_PolynomialDistribution1D(m);
    register_Axis1D(m);
    register_CartesianAxis1D(m);
    register_RadialAxis1D(m);

    register_CartesianAxisExponentialDensityDistribution(m);
    register_CartesianAxisPolynomialDensityDistribution(m);
    register_ConstantDensityDistribution(m);

    register_MaterialModel(m);
}
