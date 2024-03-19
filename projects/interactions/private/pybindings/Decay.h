#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/Decay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_Decay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<Decay, std::shared_ptr<Decay>>(m, "Decay")
        .def("TotalDecayLength",&Decay::TotalDecayLength)
        .def("TotalDecayLengthForFinalState",&Decay::TotalDecayLengthForFinalState);

}
