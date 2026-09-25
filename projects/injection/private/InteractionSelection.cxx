#include "InteractionSelection.h"

#include <set>
#include <cmath>
#include <algorithm>

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
    siren::dataclasses::InteractionRecord const & record,
    bool for_generation)
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
    siren::geometry::Geometry::IntersectionList intersections;
    std::set<siren::dataclasses::ParticleType> available_targets;
    // Decay rates do not depend on material. In particular, long-lived
    // particles can decay outside the finite detector hierarchy; querying
    // material intersections there is unnecessary and can lose precision.
    if (interactions->HasCrossSections()) {
        siren::math::Vector3D primary_direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
        primary_direction.normalize();
        intersections = detector_model->GetIntersections(
            siren::detector::DetectorPosition(interaction_vertex),
            siren::detector::DetectorDirection(primary_direction));
        available_targets = detector_model->GetAvailableTargets(
            intersections,
            siren::detector::DetectorPosition(record.interaction_vertex));
    }
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

    // Python models can change their advertised signatures after configuration.
    // Revalidate at use so a stale selection cannot silently lose support.
    if (for_generation) interactions->ValidateDecayChannels();
    // Pure decays share the same boost factor in every inverse flight length.
    // At rest that factor is singular: use widths for the discrete competition
    // on BOTH the generation and physical sides instead. Moving-parent rates
    // retain their historical arithmetic and sampling order.
    bool const rest_decay_rates = !interactions->HasCrossSections()
        && record.primary_mass > 0.0
        && record.primary_momentum[1] == 0.0
        && record.primary_momentum[2] == 0.0
        && record.primary_momentum[3] == 0.0;
    bool const validate_widths = interactions->HasDecayChannels() || rest_decay_rates;
    auto check_width = [](double width) {
        if (!std::isfinite(width) || width < 0.0) {
            throw siren::utilities::ConfigurationError(
                "Decay-channel rates require finite nonnegative decay widths");
        }
    };
    if (interactions->HasDecays()) {
        double decay_target_mass = detector_model->GetTargetMass(
            siren::dataclasses::ParticleType::Decay);
        double total_width = 0.0;
        for (auto const & decay : interactions->GetDecays()) {
            if (validate_widths) {
                // All models contribute to propagation, including models whose
                // generated signatures are excluded by the selection.
                double width = decay->TotalDecayWidthAllFinalStates(record);
                check_width(width);
                total_width += width;
                check_width(total_width); // finite terms can overflow in sum
            }
            for (auto const & signature :
                 decay->GetPossibleSignaturesFromParent(
                     record.signature.primary_type)) {
                candidate_record.signature = signature;
                double width = 0.0;
                if (validate_widths) {
                    width = decay->TotalDecayWidth(candidate_record);
                    check_width(width);
                }
                if (for_generation && !interactions->AllowsDecay(signature)) continue;
                candidates.push_back(InteractionCandidate{
                    signature,
                    decay_target_mass,
                    rest_decay_rates ? width
                        : 1.0 / (decay->TotalDecayLength(candidate_record)
                                 / siren::utilities::Constants::cm),
                    decay});
            }
        }
    }

    if (rest_decay_rates) {
        // Only ratios enter this competition. Rescale before multiplying by
        // final-state densities so tiny but finite widths do not underflow.
        double scale = 0.0;
        for (auto const & candidate : candidates) scale = std::max(scale, candidate.rate);
        if (scale > 0.0) {
            for (auto & candidate : candidates) candidate.rate /= scale;
        }
    }
    if (for_generation && (interactions->HasDecayChannels() || rest_decay_rates)) {
        for (auto const & candidate : candidates) {
            if (!std::isfinite(candidate.rate) || candidate.rate < 0.0) {
                throw siren::utilities::ConfigurationError("Selected decay generation requires finite nonnegative interaction rates");
            }
        }
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
            [](InteractionCandidate const & c) { return c.rate == 0.0; }), candidates.end());
    }
    return candidates;
}

} // namespace detail
} // namespace injection
} // namespace siren
