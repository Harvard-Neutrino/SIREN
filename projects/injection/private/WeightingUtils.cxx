#include "SIREN/injection/WeightingUtils.h"

#include <set>                                                    // for set
#include <array>                                                  // for array
#include <vector>                                                 // for vector

#include "SIREN/interactions/CrossSection.h"            // for Cro...
#include "SIREN/interactions/InteractionCollection.h"  // for Cro...
#include "SIREN/interactions/Decay.h"                   // for Decay
#include "SIREN/dataclasses/InteractionRecord.h"         // for Int...
#include "SIREN/dataclasses/InteractionSignature.h"      // for Int...
#include "SIREN/dataclasses/Particle.h"                  // for Par...
#include "SIREN/detector/DetectorModel.h"                   // for Ear...
#include "SIREN/detector/Coordinates.h"
#include "SIREN/geometry/Geometry.h"                     // for Geo...
#include "SIREN/math/Vector3D.h"                         // for Vec...
#include "SIREN/utilities/Constants.h"                   // for cm

namespace siren {
namespace injection {

using detector::DetectorPosition;
using detector::DetectorDirection;

double CrossSectionProbability(std::shared_ptr<siren::detector::DetectorModel const> detector_model, std::shared_ptr<siren::interactions::InteractionCollection const> interactions, siren::dataclasses::InteractionRecord const & record) {

    siren::math::Vector3D interaction_vertex(
            record.interaction_vertex[0],
            record.interaction_vertex[1],
            record.interaction_vertex[2]);

    siren::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    primary_direction.normalize();

    siren::geometry::Geometry::IntersectionList intersections = detector_model->GetIntersections(DetectorPosition(interaction_vertex), DetectorDirection(primary_direction));
    std::set<siren::dataclasses::ParticleType> const & possible_targets = interactions->TargetTypes();
    std::set<siren::dataclasses::ParticleType> available_targets_list = detector_model->GetAvailableTargets(intersections, DetectorPosition(record.interaction_vertex));
    std::set<siren::dataclasses::ParticleType> available_targets(available_targets_list.begin(), available_targets_list.end());

    double total_prob = 0.0;
    // double selected_prob = 0.0;
    double selected_final_state = 0.0;
    siren::dataclasses::InteractionRecord fake_record = record;
    // first let's check decays
    std::vector<std::shared_ptr<siren::interactions::Decay>> decays = interactions->GetDecays();
    for(auto const & decay : decays) {
      std::vector<siren::dataclasses::InteractionSignature> signatures = decay->GetPossibleSignaturesFromParent(record.signature.primary_type);
      for(auto const & signature : signatures) {
        fake_record.signature = signature;
        // fake_prob has units of 1/cm to match cross section probabilities
        double decay_prob = 1./(decay->TotalDecayLength(fake_record)/siren::utilities::Constants::cm);
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
            double target_density = detector_model->GetParticleDensity(intersections, DetectorPosition(interaction_vertex), target);
            // Loop over cross sections that have this target
            std::vector<std::shared_ptr<siren::interactions::CrossSection>> const & target_cross_sections = interactions->GetCrossSectionsForTarget(target);
            for(auto const & cross_section : target_cross_sections) {
                // Loop over cross section signatures with the same target
                std::vector<siren::dataclasses::InteractionSignature> signatures = cross_section->GetPossibleSignaturesFromParents(record.signature.primary_type, target);
                for(auto const & signature : signatures) {
                    fake_record.signature = signature;
                    fake_record.target_mass = detector_model->GetTargetMass(target);
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
    if (total_prob == 0)
        return 0.0;
    return selected_final_state / total_prob;
}

} // namespace injection
} // namespace siren

