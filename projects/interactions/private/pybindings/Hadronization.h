#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/Hadronization.h"

#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

using namespace pybind11;
using namespace siren::interactions;


void register_Hadronization(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<Hadronization, std::shared_ptr<Hadronization>>(m, "Hadronization")
        // .def(init<>())
        .def("__eq__", [](const Hadronization &self, const Hadronization &other){ return self == other; })
        .def("equal", &Hadronization::equal)
        .def("SampleFinalState", &Hadronization::SampleFinalState)
        .def("GetPossibleSignatures", &Hadronization::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent", &Hadronization::GetPossibleSignaturesFromParent)
        .def("FragmentationFraction", &Hadronization::FragmentationFraction) 
        ;
}
