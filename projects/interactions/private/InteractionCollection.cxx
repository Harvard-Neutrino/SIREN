#include "SIREN/interactions/InteractionCollection.h"

#include <cmath>
#include <algorithm>
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"

#include <map>                                                // for map
#include <set>                                                // for operator==
#include <tuple>                                              // for tie
#include <limits>                                             // for numeric...
#include <vector>                                             // for vector
#include <utility>                                            // for pair

#include "SIREN/interactions/Interaction.h"          // for Interaction
#include "SIREN/interactions/CrossSection.h"        // for CrossSe...
#include "SIREN/interactions/Decay.h"               // for Decay
#include "SIREN/dataclasses/InteractionRecord.h"     // for Interac...
#include "SIREN/dataclasses/InteractionSignature.h"  // for Interac...
#include "SIREN/dataclasses/Particle.h"              // for Particle

namespace siren {
namespace interactions {

void InteractionCollection::InitializeTargetTypes() {
    target_types.clear();
    cross_sections_by_target.clear();
    for(unsigned int i=0; i<cross_sections.size(); ++i) {
        // Gather target types
        std::vector<siren::dataclasses::ParticleType> xs_targets = cross_sections[i]->GetPossibleTargets();
        //target_types.reserve(target_types.size() + std::distance(xs_targets.begin(), xs_targets.end()));
        for(auto xs : xs_targets)
            target_types.insert(xs);

        // Track cross sections by their target type
        for(unsigned int j=0; j<xs_targets.size(); ++j) {
            siren::dataclasses::ParticleType target = xs_targets[j];
            std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(target);
            if(it == cross_sections_by_target.end()) {
                cross_sections_by_target.insert(it, std::make_pair(target, std::vector<std::shared_ptr<CrossSection>>{cross_sections[i]}));
            } else {
                cross_sections_by_target[target].push_back(cross_sections[i]);
            }
        }
    }

    // Remove duplicate target types
    // std::set<siren::dataclasses::ParticleType> target_set(target_types.begin(), target_types.end());
    // target_types.resize(target_set.size());
    // std::copy(target_set.begin(), target_set.end(), target_types.begin());
}

const std::vector<std::shared_ptr<CrossSection>> InteractionCollection::empty = {};

InteractionCollection::InteractionCollection() {}

InteractionCollection::InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections) : primary_type(primary_type), cross_sections(cross_sections) {
    InitializeTargetTypes();
}

InteractionCollection::InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<Decay>> decays) : primary_type(primary_type), decays(decays) {
    InitializeTargetTypes();
}

InteractionCollection::InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::vector<std::shared_ptr<Decay>> decays) : primary_type(primary_type), cross_sections(cross_sections), decays(decays) {
    InitializeTargetTypes();
}

InteractionCollection::InteractionCollection(siren::dataclasses::ParticleType primary_type, std::vector<std::shared_ptr<Interaction>> interactions) : primary_type(primary_type) {
    for(auto interaction : interactions) {
        std::shared_ptr<CrossSection> xs = std::dynamic_pointer_cast<CrossSection>(interaction);
        if(xs) {
            cross_sections.push_back(xs);
        } else {
            std::shared_ptr<Decay> dec = std::dynamic_pointer_cast<Decay>(interaction);
            if(dec) {
                decays.push_back(dec);
            } else {
                throw std::runtime_error("InteractionCollection: Interaction is neither a CrossSection nor a Decay");
            }
        }
    }
    InitializeTargetTypes();
}

bool InteractionCollection::operator==(InteractionCollection const & other) const {
    return
        std::tie(primary_type, target_types, cross_sections, decays, decay_channels)
        ==
        std::tie(other.primary_type, other.target_types, other.cross_sections, other.decays, other.decay_channels);
}

std::vector<std::shared_ptr<CrossSection>> const & InteractionCollection::GetCrossSectionsForTarget(siren::dataclasses::ParticleType p) const {
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(p);
    if(it != cross_sections_by_target.end()) {
        return it->second;
    } else {
        return empty;
    }
}

void InteractionCollection::SetDecayChannels(
        std::optional<std::vector<dataclasses::InteractionSignature>> channels) {
    if(channels) std::sort(channels->begin(), channels->end());
    auto previous = decay_channels;
    decay_channels = std::move(channels);
    try {
        ValidateDecayChannels();
    } catch(...) {
        decay_channels = std::move(previous);
        throw;
    }
}

bool InteractionCollection::AllowsDecay(dataclasses::InteractionSignature const & signature) const {
    return !decay_channels || std::find(decay_channels->begin(), decay_channels->end(), signature) != decay_channels->end();
}

void InteractionCollection::ValidateDecayChannels() const {
    if(!decay_channels) return;
    if(decay_channels->empty()) {
        throw siren::utilities::ConfigurationError("decay_channels must be non-empty; use None for all decays");
    }
    std::set<dataclasses::InteractionSignature> available;
    for(auto const & decay : decays) {
        for(auto const & signature : decay->GetPossibleSignaturesFromParent(primary_type)) available.insert(signature);
    }
    std::set<dataclasses::InteractionSignature> seen;
    for(auto const & signature : *decay_channels) {
        if(signature.primary_type != primary_type || signature.target_type != dataclasses::ParticleType::Decay || !available.count(signature)) {
            throw siren::utilities::ConfigurationError("decay_channels contains a signature absent from this parent's decay models");
        }
        if(!seen.insert(signature).second) {
            throw siren::utilities::ConfigurationError("decay_channels contains a duplicate signature");
        }
    }
}

double InteractionCollection::ParentDecayWidth(dataclasses::InteractionRecord const & record) const {
    double total = 0, owned = 0;
    for (auto const & decay : decays) {
        double declared = decay->ParentDecayWidth(record);
        if (!std::isfinite(declared) || declared < 0)
            throw std::invalid_argument("Parent decay width must be finite and nonnegative");
        if (declared > 0) {
            if (total > 0 && std::abs(total-declared) > 1e-12*std::max(total,declared))
                throw std::invalid_argument("Conflicting parent decay widths in one InteractionCollection");
            total = declared;
        }
    }
    // Preserve legacy Fixed workflows that never queried all-state widths.
    // Only an explicit parent declaration enables this validation contract.
    if (total == 0) return 0;
    for (auto const & decay : decays) {
        double width = decay->TotalDecayWidthAllFinalStates(record);
        if (!std::isfinite(width) || width < 0)
            throw std::invalid_argument("Owned decay widths must be finite and nonnegative");
        owned += width;
    }
    if (!std::isfinite(owned) || owned > total*(1+1e-12))
        throw std::invalid_argument("Owned decay widths exceed the declared parent width");
    return total;
}

double InteractionCollection::TotalDecayWidthAllFinalStates(dataclasses::InteractionRecord const & record) const {
  double width = ParentDecayWidth(record);
  if(width > 0) return width;
  if(!HasDecays()) return width;
  for(auto dec : decays) {
    width += dec->TotalDecayWidthAllFinalStates(record);
  }
  return width;
}

double InteractionCollection::TotalDecayLengthAllFinalStates(dataclasses::InteractionRecord const & record) const {
  double parent_width = ParentDecayWidth(record);
  if (parent_width > 0) {
    auto const & p = record.primary_momentum;
    double momentum = std::hypot(p[1],p[2],p[3]);
    if (!(record.primary_mass > 0)) throw std::invalid_argument("Decaying parent needs positive mass");
    return momentum / record.primary_mass * siren::utilities::Constants::hbarc / parent_width;
  }
  double inv_length = 0;
  if(!HasDecays()) return std::numeric_limits<double>::infinity();
  for(auto dec : decays) {
    inv_length += 1./dec->TotalDecayLengthAllFinalStates(record);
  }
  if(inv_length == 0) return std::numeric_limits<double>::infinity();
  return 1./inv_length;
}

bool InteractionCollection::MatchesPrimary(dataclasses::InteractionRecord const & record) const {
    return primary_type == record.signature.primary_type;
}

siren::dataclasses::ParticleType InteractionCollection::GetPrimaryType() const {
    return primary_type;
}

void InteractionCollection::SetPrimaryType(siren::dataclasses::ParticleType primary_type) {
    auto previous = this->primary_type;
    this->primary_type = primary_type;
    try {
        ValidateDecayChannels();
    } catch(...) {
        this->primary_type = previous;
        throw;
    }
}

std::map<siren::dataclasses::ParticleType, double> InteractionCollection::TotalCrossSectionByTarget(siren::dataclasses::InteractionRecord const & record) const {
    std::map<siren::dataclasses::ParticleType, double> result;
    for(siren::dataclasses::ParticleType target : target_types) {
        siren::dataclasses::InteractionRecord fake_record = record;
        fake_record.signature.target_type = target;

        double total_xs = 0;
        for(auto xs : cross_sections_by_target.at(target)) {
            total_xs += xs->TotalCrossSection(fake_record);
        }
        result.insert(std::make_pair(target, total_xs));
    }
    return result;
}

std::map<siren::dataclasses::ParticleType, double> InteractionCollection::TotalCrossSectionByTargetAllFinalStates(siren::dataclasses::InteractionRecord const & record) const {
    std::map<siren::dataclasses::ParticleType, double> result;
    for(siren::dataclasses::ParticleType target : target_types) {
        siren::dataclasses::InteractionRecord fake_record = record;
        fake_record.signature.target_type = target;

        double total_xs = 0;
        for(auto xs : cross_sections_by_target.at(target)) {
            total_xs += xs->TotalCrossSectionAllFinalStates(fake_record);
        }
        result.insert(std::make_pair(target, total_xs));
    }
    return result;
}

} // namespace interactions
} // namespace siren
