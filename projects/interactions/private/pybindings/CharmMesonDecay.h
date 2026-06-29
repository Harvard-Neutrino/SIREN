#include <set>
#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "../../public/SIREN/interactions/CrossSection.h"
#include "../../public/SIREN/interactions/Decay.h"
#include "../../public/SIREN/interactions/CharmMesonDecay.h"

#include "../../public/SIREN/interactions/pyDecay.h"
#include "../../../dataclasses/public/SIREN/dataclasses/Particle.h"
#include "../../../geometry/public/SIREN/geometry/Geometry.h"
#include "../../../utilities/public/SIREN/utilities/Random.h"

void register_CharmMesonDecay(pybind11::module_ & m) {
    using namespace pybind11;
    using namespace siren::interactions;

    class_<CharmMesonDecay, std::shared_ptr<CharmMesonDecay>, Decay> charmmesondecay(m, "CharmMesonDecay");

    charmmesondecay

        .def(init<>())
        .def(init<siren::dataclasses::Particle::ParticleType>(),
                arg("primary_type"))
        .def(self == self)
        .def("SampleFinalState",&CharmMesonDecay::SampleFinalState)
        .def("GetPossibleSignatures",&CharmMesonDecay::GetPossibleSignatures)
        .def("GetPossibleSignaturesFromParent",&CharmMesonDecay::GetPossibleSignaturesFromParent);

}
