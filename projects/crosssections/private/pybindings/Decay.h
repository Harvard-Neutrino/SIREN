#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/LeptonInjector/crosssections/CrossSection.h"
#include "../../public/LeptonInjector/crosssections/Decay.h"
#include "../../../dataclasses/public/LeptonInjector/dataclasses/Particle.h"
#include "../../../geometry/public/LeptonInjector/geometry/Geometry.h"
#include "../../../utilities/public/LeptonInjector/utilities/Random.h"

void register_Decay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace LI::crosssections;

    class_<Decay, std::shared_ptr<Decay>>(m, "Decay")
        .def("TotalDecayLength",&Decay::TotalDecayLength)
        .def("TotalDecayLengthForFinalState",&Decay::TotalDecayLengthForFinalState);

}
