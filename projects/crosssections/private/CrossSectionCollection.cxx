#include <map>
#include <set>
#include <vector>
#include <string>
#include <memory>

#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/CrossSectionCollection.h"

namespace LI {
namespace crosssections {

void CrossSectionCollection::InitializeTargetTypes() {
    target_types.clear();
    cross_sections_by_target.clear();
    for(unsigned int i=0; i<cross_sections.size(); ++i) {
        // Gather target types
        std::vector<LI::dataclasses::Particle::ParticleType> xs_targets = cross_sections[i]->GetPossibleTargets();
        //target_types.reserve(target_types.size() + std::distance(xs_targets.begin(), xs_targets.end()));
        for(auto xs : xs_targets)
            target_types.insert(xs);

        // Track cross sections by their target type
        for(unsigned int j=0; j<xs_targets.size(); ++j) {
            LI::dataclasses::Particle::ParticleType target = xs_targets[j];
            std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(target);
            if(it == cross_sections_by_target.end()) {
                cross_sections_by_target.insert(it, std::make_pair(target, std::vector<std::shared_ptr<CrossSection>>{cross_sections[i]}));
            } else {
                cross_sections_by_target[target].push_back(cross_sections[i]);
            }
        }
    }

    // Remove duplicate target types
    // std::set<LI::dataclasses::Particle::ParticleType> target_set(target_types.begin(), target_types.end());
    // target_types.resize(target_set.size());
    // std::copy(target_set.begin(), target_set.end(), target_types.begin());
}

const std::vector<std::shared_ptr<CrossSection>> CrossSectionCollection::empty = {};

CrossSectionCollection::CrossSectionCollection() {}

CrossSectionCollection::CrossSectionCollection(LI::dataclasses::Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections) : primary_type(primary_type), cross_sections(cross_sections) {
    InitializeTargetTypes();
}

CrossSectionCollection::CrossSectionCollection(LI::dataclasses::Particle::ParticleType primary_type, std::vector<std::shared_ptr<Decay>> decays) : primary_type(primary_type), decays(decays) {
    InitializeTargetTypes();
}

CrossSectionCollection::CrossSectionCollection(LI::dataclasses::Particle::ParticleType primary_type, std::vector<std::shared_ptr<CrossSection>> cross_sections, std::vector<std::shared_ptr<Decay>> decays) : primary_type(primary_type), cross_sections(cross_sections), decays(decays) {
    InitializeTargetTypes();
}

bool CrossSectionCollection::operator==(CrossSectionCollection const & other) const {
    return
        std::tie(primary_type, target_types, cross_sections, decays)
        ==
        std::tie(other.primary_type, other.target_types, other.cross_sections, other.decays);
}

std::vector<std::shared_ptr<CrossSection>> const & CrossSectionCollection::GetCrossSectionsForTarget(LI::dataclasses::Particle::ParticleType p) const {
    std::map<LI::dataclasses::Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>>::const_iterator it = cross_sections_by_target.find(p);
    if(it != cross_sections_by_target.end()) {
        return it->second;
    } else {
        return empty;
    }
}

double CrossSectionCollection::TotalDecayWidth(dataclasses::InteractionRecord const & record) const {
  double width = 0;
  if(!HasDecays()) return width;
  for(auto dec : decays) {
    width += dec->TotalDecayWidth(record);
  }
  return width;
}

bool CrossSectionCollection::MatchesPrimary(dataclasses::InteractionRecord const & record) const {
    return primary_type == record.signature.primary_type;
}

} // namespace crosssections
} // namespace LI
