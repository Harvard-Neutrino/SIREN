#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/Hadronization.h"
#include "../../public/SIREN/interactions/CharmHadronization.h"


#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionRecord.h"
#include "../../../dataclasses/public/SIREN/dataclasses/InteractionSignature.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

using namespace pybind11;
using namespace siren::interactions;


void register_CharmHadronization(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<CharmHadronization, std::shared_ptr<CharmHadronization>, Hadronization> charmhadronization(m, "CharmHadronization");

    charmhadronization

        .def(init<>())
        .def(self == self)
        .def("SampleFinalState",&CharmHadronization::SampleFinalState)
        .def("GetPossibleSignatures",&CharmHadronization::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&CharmHadronization::GetPossibleSignaturesFromParent)
        .def("FragmentationFraction",&CharmHadronization::FragmentationFraction);
        ;
}
