#include <set>                                    // for set
#include <memory>                                 // for shared_ptr
#include <string>                                 // for string
#include <vector>                                 // for vector
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include "SIREN/interactions/pyDarkNewsDecay.h"
#include "SIREN/interactions/DarkNewsDecay.h"

#include "SIREN/interactions/Decay.h"   // for Decay
#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/utilities/Pybind11Trampoline.h" // for Pybind11Trampoline

#include "SIREN/dataclasses/InteractionSignature.h"
#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"

namespace siren {
namespace interactions {

// Trampoline class for DarkNewsDecay
    pyDarkNewsDecay::pyDarkNewsDecay(DarkNewsDecay && parent) : DarkNewsDecay(std::move(parent)) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsDecay)))));
    }
    pyDarkNewsDecay::pyDarkNewsDecay(DarkNewsDecay const & parent) : DarkNewsDecay(parent) {
        self = pybind11::reinterpret_borrow<pybind11::object>(pybind11::handle(get_object_handle(&parent, pybind11::detail::get_type_info(typeid(DarkNewsDecay)))));
    }
    //pybind11::object self;

double pyDarkNewsDecay::TotalDecayWidth(dataclasses::InteractionRecord const & interaction) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            std::cref(interaction)
        )
    }

double pyDarkNewsDecay::TotalDecayWidthForFinalState(dataclasses::InteractionRecord const & interaction) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidthForFinalState,
            "TotalDecayWidthForFinalState",
            std::cref(interaction)
        )
    }

double pyDarkNewsDecay::TotalDecayWidth(siren::dataclasses::ParticleType primary) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            TotalDecayWidth,
            "TotalDecayWidth",
            primary
        )
    }

double pyDarkNewsDecay::DifferentialDecayWidth(dataclasses::InteractionRecord const & interaction) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            DifferentialDecayWidth,
            "DifferentialDecayWidth",
            std::cref(interaction)
        )
    }

void pyDarkNewsDecay::SampleRecordFromDarkNews(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleRecordFromDarkNews,
            "SampleRecordFromDarkNews",
            std::ref(record),
            random
        )
    }

void pyDarkNewsDecay::SampleFinalState(dataclasses::CrossSectionDistributionRecord & record, std::shared_ptr<siren::utilities::SIREN_random> random) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            void,
            SampleFinalState,
            "SampleFinalState",
            std::ref(record),
            random
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> pyDarkNewsDecay::GetPossibleSignatures() const {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignatures,
            "GetPossibleSignatures"
        )
    }

    std::vector<siren::dataclasses::InteractionSignature> pyDarkNewsDecay::GetPossibleSignaturesFromParent(siren::dataclasses::ParticleType primary_type) const {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<siren::dataclasses::InteractionSignature>,
            GetPossibleSignaturesFromParent,
            "GetPossibleSignaturesFromParent",
            primary_type
        )
    }

    std::vector<std::string> pyDarkNewsDecay::DensityVariables() const {
        SELF_OVERRIDE_PURE(
            self,
            DarkNewsDecay,
            std::vector<std::string>,
            DensityVariables,
            "DensityVariables"
        )
    }

double pyDarkNewsDecay::FinalStateProbability(dataclasses::InteractionRecord const & record) const {
        SELF_OVERRIDE(
            self,
            DarkNewsDecay,
            double,
            FinalStateProbability,
            "FinalStateProbability",
            std::cref(record)
        )
    }

} // namespace interactions
} // namespace siren

