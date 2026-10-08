#pragma once
#ifndef SIREN_VertexWeightingMode_H
#define SIREN_VertexWeightingMode_H

#include <string>
#include <cstdint>                            // for uint32_t
#include <stdexcept>                          // for runtime_error
#include <cereal/cereal.hpp>                  // for make_nvp, CEREAL_CLASS_VERSION

namespace siren {
namespace dataclasses {

// Controls how the Weighter computes probability factors at a vertex.
//
// The event weight at each vertex is a product of factors:
//   InteractionProbability * PositionProbability * FinalStateProbability
//
// Not all factors are meaningful for every vertex type. A dk2nu pion
// decays at a fixed point -- InteractionProbability is meaningless.
// A gamma-ray segment provides external bounds for the interaction
// probability integral. This struct lets the user control which
// factors are computed and where the integration bounds come from.
//
// Note: even a Fixed vertex still charges the channel-selection
// probability when multiple channels compete, because the injector
// rate-selects the channel independently of this mode; only the
// path-interaction and position factors are suppressed here.
struct VertexWeightingMode {

    // Which probability factors to compute at this vertex.
    // When false, the factor is treated as 1.0 (no contribution).
    bool compute_interaction_probability = true;
    bool compute_position_probability = true;

    // Where the integration bounds come from for the interaction
    // probability integral.
    enum class BoundSource {
        Geometry,      // Bounds from detector geometry intersections
        Distribution,  // Bounds from the vertex position distribution
        None           // No bounds (vertex is externally fixed)
        // Mixed bounds (an externally supplied starting point combined with a geometry-derived end
        // point, as dk2nu neutrino propagation would need) are not represented by this enum.
    };
    BoundSource bound_source = BoundSource::Geometry;

    // Include the survival from the record's creation point
    // (primary_initial_position) to the entry of its injection bounds in the
    // physical weight, as exp(-depth). A geometry-bounded vertex is otherwise
    // weighted from that entry: its interaction and position factors give the
    // probability of interacting inside the bounds once the particle has
    // reached them. Requires both factors and BoundSource::Geometry.
    bool survival_from_creation = false;

    bool operator==(VertexWeightingMode const & o) const {
        return compute_interaction_probability == o.compute_interaction_probability
            && compute_position_probability == o.compute_position_probability
            && bound_source == o.bound_source
            && survival_from_creation == o.survival_from_creation;
    }
    bool operator!=(VertexWeightingMode const & o) const {
        return !(*this == o);
    }

    // ---- Named presets ----

    // Standard propagation: particle travels through geometry,
    // interaction probability computed from cross section * density.
    // This is the default and matches all existing SIREN behavior.
    static VertexWeightingMode Propagated() {
        return {true, true, BoundSource::Geometry, false};
    }

    // Propagated, with the survival from the creation point to the entry of
    // the injection bounds included in the weight. For a long-lived particle
    // that must cross many decay or interaction lengths before reaching its
    // bounds, this is the physical rate rather than the rate given arrival.
    static VertexWeightingMode PropagatedFromCreation() {
        return {true, true, BoundSource::Geometry, true};
    }

    // Fixed vertex: position externally determined (input from dk2nu mesons, etc.).
    // No interaction or position probability computed.
    static VertexWeightingMode Fixed() {
        return {false, false, BoundSource::None, false};
    }

    // External bounds: interaction and position probabilities computed
    // within bounds provided by the primary vertex distribution itself, e.g.
    // a PrimaryExternalDistribution in segment mode ("length" column: gamma
    // steps in a production target). Injection bounds always come from the
    // vertex distribution's InjectionBounds; the Injector requires this mode
    // for distributions whose ProvidesExternalBounds() is true and rejects
    // it for all others.
    static VertexWeightingMode ExternalBounds() {
        return {true, true, BoundSource::Distribution, false};
    }

    // Throws if this combination of fields cannot be weighted.
    void Validate() const {
        if(survival_from_creation
           && !(compute_interaction_probability && compute_position_probability
                && bound_source == BoundSource::Geometry)) {
            throw std::runtime_error("VertexWeightingMode: survival_from_creation requires "
                "compute_interaction_probability, compute_position_probability and "
                "BoundSource::Geometry (use PropagatedFromCreation())");
        }
    }

    // ---- Serialization ----
    // bound_source is archived as its underlying int so the enum encoding
    // is stable across compilers. Version 1 adds survival_from_creation;
    // version-0 archives load with it false.
    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if(version <= 1) {
            archive(::cereal::make_nvp("ComputeInteractionProbability", compute_interaction_probability));
            archive(::cereal::make_nvp("ComputePositionProbability", compute_position_probability));
            archive(::cereal::make_nvp("BoundSource", static_cast<int>(bound_source)));
            if(version >= 1) {
                archive(::cereal::make_nvp("SurvivalFromCreation", survival_from_creation));
            }
        } else {
            throw std::runtime_error("VertexWeightingMode only supports version <= 1!");
        }
    }
    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if(version <= 1) {
            int bound_source_int;
            archive(::cereal::make_nvp("ComputeInteractionProbability", compute_interaction_probability));
            archive(::cereal::make_nvp("ComputePositionProbability", compute_position_probability));
            archive(::cereal::make_nvp("BoundSource", bound_source_int));
            if(bound_source_int != static_cast<int>(BoundSource::Geometry)
               && bound_source_int != static_cast<int>(BoundSource::Distribution)
               && bound_source_int != static_cast<int>(BoundSource::None)) {
                throw std::runtime_error("VertexWeightingMode: invalid BoundSource value "
                    + std::to_string(bound_source_int) + " in archive");
            }
            bound_source = static_cast<BoundSource>(bound_source_int);
            survival_from_creation = false;
            if(version >= 1) {
                archive(::cereal::make_nvp("SurvivalFromCreation", survival_from_creation));
            }
            Validate();
        } else {
            throw std::runtime_error("VertexWeightingMode only supports version <= 1!");
        }
    }
};

std::string VertexWeightingModeName(VertexWeightingMode const & mode);

} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::VertexWeightingMode, 1);

#endif // SIREN_VertexWeightingMode_H
