#pragma once
#include <memory>
#include <array>
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include "SIREN/dataclasses/InteractionTree.h"

namespace siren { namespace injection {

// Declarative secondary expansion, using the same true=stop polarity as the
// legacy callback. Each rule is (parent type, child type,
// child index or -1 wildcard, exclusive maximum parent depth or -1 unlimited, any-child flag).
// Rules are OR'ed. No matching rule means stop. There are no Python callbacks.
class SecondaryExpansion {
public:
    using Rule = std::array<std::int64_t, 5>;
    explicit SecondaryExpansion(std::vector<Rule> rules = {}) : rules_(std::move(rules)) { Validate(); }
    bool ShouldStop(siren::dataclasses::InteractionTree const & tree,
                    std::shared_ptr<siren::dataclasses::InteractionTreeDatum> parent,
                    std::size_t index) const {
        if (!parent || index >= parent->record.signature.secondary_types.size())
            throw std::out_of_range("SecondaryExpansion child index is outside the parent signature");
        for (auto const & rule : rules_) {
            if (rule[0] != static_cast<int>(parent->record.signature.primary_type)) continue;
            if (!rule[4] && rule[1] != static_cast<int>(parent->record.signature.secondary_types[index])) continue;
            if (rule[2] >= 0 && static_cast<std::size_t>(rule[2]) != index) continue;
            if (rule[3] >= 0 && parent->depth(tree) >= rule[3]) continue;
            return false;
        }
        return true;
    }
    bool operator==(SecondaryExpansion const & other) const { return rules_ == other.rules_; }
    std::vector<Rule> const & GetRules() const { return rules_; }
    template<class Archive> void serialize(Archive & archive, std::uint32_t version) {
        if (version != 1) throw std::runtime_error("Legacy SecondaryExpansion wildcard semantics are ambiguous; regenerate the archive");
        archive(cereal::make_nvp("Rules", rules_));
        Validate();
    }
private:
    void Validate() const {
        for (auto const & r : rules_)
            if (r[0] == static_cast<int>(siren::dataclasses::ParticleType::unknown)
                || r[2] < -1 || r[3] < -1 || (r[4]!=0 && r[4]!=1))
                throw std::invalid_argument("Invalid SecondaryExpansion parent, index or depth");
    }
    std::vector<Rule> rules_;
};
}}
CEREAL_CLASS_VERSION(siren::injection::SecondaryExpansion, 1);
