#include "LeptonInjector/injection/WeightingUtils.h"

#include <set>                                                    // for set
#include <array>                                                  // for array
#include <vector>                                                 // for vector

#include "LeptonInjector/interactions/CrossSection.h"            // for Cro...
#include "LeptonInjector/interactions/InteractionCollection.h"  // for Cro...
#include "LeptonInjector/interactions/Decay.h"                   // for Decay
#include "LeptonInjector/dataclasses/InteractionRecord.h"         // for Int...
#include "LeptonInjector/dataclasses/InteractionSignature.h"      // for Int...
#include "LeptonInjector/dataclasses/Particle.h"                  // for Par...
#include "LeptonInjector/detector/EarthModel.h"                   // for Ear...
#include "LeptonInjector/geometry/Geometry.h"                     // for Geo...
#include "LeptonInjector/math/Vector3D.h"                         // for Vec...
#include "LeptonInjector/utilities/Constants.h"                   // for cm

namespace LI {
namespace injection {

double CrossSectionProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::interactions::InteractionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) {
    std::set<LI::dataclasses::Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();
    std::set<LI::dataclasses::Particle::ParticleType> available_targets_list = earth_model->GetAvailableTargets(earth_model->GetEarthCoordPosFromDetCoordPos(record.interaction_vertex));
    std::set<LI::dataclasses::Particle::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    LI::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    LI::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    LI::geometry::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), earth_model->GetEarthCoordDirFromDetCoordDir(primary_direction));

    double total_prob = 0.0;
    // double selected_prob = 0.0;
    double selected_final_state = 0.0;
    LI::dataclasses::InteractionRecord fake_record = record;
    // first let's check decays
    std::vector<std::shared_ptr<LI::interactions::Decay>> decays = cross_sections->GetDecays();
    for(auto const & decay : decays) {
      std::vector<LI::dataclasses::InteractionSignature> signatures = decay->GetPossibleSignaturesFromParent(record.signature.primary_type);
      for(auto const & signature : signatures) {
        fake_record.signature = signature;
        // fake_prob has units of 1/cm to match cross section probabilities
        double decay_prob = 1./(decay->TotalDecayLengthForFinalState(fake_record)/LI::utilities::Constants::cm);
        total_prob += decay_prob;
        if(signature == record.signature) {
            selected_final_state += decay_prob * decay->FinalStateProbability(record);
        }
      }
    }
    // next let's check cross sections
    for(auto const target : available_targets) {
        if(possible_targets.find(target) != possible_targets.end()) {
            // Get target density
            double target_density = earth_model->GetParticleDensity(intersections, earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<LI::interactions::CrossSection>> const & target_cross_sections = cross_sections->GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<LI::dataclasses::InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                for(auto const & signature : signatures) {
                    fake_record.signature = signature;
                    fake_record.target_mass = earth_model->GetTargetMass(target);
                    fake_record.target_momentum = {fake_record.target_mass,0,0,0};
                    // Add total cross section times density to the total prob
                    double target_prob = target_density * cross_section->TotalCrossSection(fake_record);
                    total_prob += target_prob;
                    // Add up total cross section times density times final state prob for matching signatures
                    if(signature == record.signature) {
                        // selected_prob += target_prob;
                        selected_final_state += target_prob * cross_section->FinalStateProbability(record);
                    }
                }
            }
        }
    }
    return selected_final_state / total_prob;
}

} // namespace injection
} // namespace LI

