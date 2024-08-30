#include "SIREN/interactions/InteractionCollection.h"

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
        std::tie(primary_type, target_types, cross_sections, decays)
        ==
        std::tie(other.primary_type, other.target_types, other.cross_sections, other.decays);
}

std::vector<std::shared_ptr<CrossSection>> const & InteractionCollection::GetCrossSectionsForTarget(siren::dataclasses::ParticleType p) const {
    std::map<siren::dataclasses::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(p);
    if(it != cross_sections_by_target.end()) {
        return it->second;
    } else {
        return empty;
    }
}

double InteractionCollection::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
  double width = 0;
  if(!HasDecays()) return width;
  for(auto dec : decays) {
    width += dec->TotalDecayWidth(record);
  }
  return width;
}

double InteractionCollection::TotalDecayLength(dataclasses::InteractionRecord const & record) const {
  double inv_length = 0;
  if(!HasDecays()) return std::numeric_limits<double>::infinity();
  for(auto dec : decays) {
    inv_length += 1./dec->TotalDecayLength(record);
  }
  return 1./inv_length;
}

bool InteractionCollection::MatchesPrimary(dataclasses::InteractionRecord const & record) const {
    return primary_type == record.signature.primary_type;
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
