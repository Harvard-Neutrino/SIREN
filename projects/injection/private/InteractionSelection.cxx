#include "InteractionSelection.h"

#include <set>

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/Coordinates.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/interactions/CrossSection.h"
#include "SIREN/interactions/Decay.h"
#include "SIREN/interactions/InteractionCollection.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Constants.h"
#include "SIREN/utilities/Errors.h"

namespace siren {
namespace injection {
namespace detail {

std::vector<InteractionCandidate> EnumerateInteractionCandidates(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    std::shared_ptr<siren::interactions::InteractionCollection const> interactions,
    siren::dataclasses::InteractionRecord const & record)
{
    if (!interactions) {
        throw siren::utilities::ConfigurationError(
            "Interaction selection requires an InteractionCollection "
            "[siren-docs: errors#configuration]");
    }
    if (!detector_model) {
        throw siren::utilities::ConfigurationError(
            "Interaction selection requires a DetectorModel "
            "[siren-docs: errors#configuration]");
    }

    siren::math::Vector3D interaction_vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);
    siren::math::Vector3D primary_direction(
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3]);
    primary_direction.normalize();

    siren::geometry::Geometry::IntersectionList intersections =
        detector_model->GetIntersections(
            siren::detector::DetectorPosition(interaction_vertex),
            siren::detector::DetectorDirection(primary_direction));
    std::set<siren::dataclasses::ParticleType> available_targets =
        detector_model->GetAvailableTargets(
            intersections,
            siren::detector::DetectorPosition(record.interaction_vertex));
    std::set<siren::dataclasses::ParticleType> const & possible_targets =
        interactions->TargetTypes();

    std::vector<InteractionCandidate> candidates;
    siren::dataclasses::InteractionRecord candidate_record = record;

    // Preserve the historical sampling order: cross sections first, decays
    // second. The probability is carried by each entry, so model ordering
    // affects only deterministic RNG interval placement, never precedence.
    for (auto const target : available_targets) {
        if (possible_targets.find(target) == possible_targets.end()) continue;
        double target_density = detector_model->GetParticleDensity(
            intersections,
            siren::detector::DetectorPosition(interaction_vertex), target);
        double target_mass = detector_model->GetTargetMass(target);
        for (auto const & cross_section :
             interactions->GetCrossSectionsForTarget(target)) {
            for (auto const & signature :
                 cross_section->GetPossibleSignaturesFromParents(
                     record.signature.primary_type, target)) {
                candidate_record.signature = signature;
                candidate_record.target_mass = target_mass;
                candidates.push_back(InteractionCandidate{
                    signature,
                    target_mass,
                    target_density
                        * cross_section->TotalCrossSection(candidate_record),
                    cross_section});
            }
        }
    }

    if (interactions->HasDecays()) {
        double decay_target_mass = detector_model->GetTargetMass(
            siren::dataclasses::ParticleType::Decay);
        for (auto const & decay : interactions->GetDecays()) {
            for (auto const & signature :
                 decay->GetPossibleSignaturesFromParent(
                     record.signature.primary_type)) {
                candidate_record.signature = signature;
                candidates.push_back(InteractionCandidate{
                    signature,
                    decay_target_mass,
                    1.0 / (decay->TotalDecayLength(candidate_record)
                           / siren::utilities::Constants::cm),
                    decay});
            }
        }
    }

    return candidates;
}

} // namespace detail
} // namespace injection
} // namespace siren
