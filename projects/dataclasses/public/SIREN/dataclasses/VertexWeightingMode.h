#pragma once
#ifndef SIREN_VertexWeightingMode_H
#define SIREN_VertexWeightingMode_H

#include <string>

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
    };
    BoundSource bound_source = BoundSource::Geometry;

    bool operator==(VertexWeightingMode const & o) const {
        return compute_interaction_probability == o.compute_interaction_probability
            && compute_position_probability == o.compute_position_probability
            && bound_source == o.bound_source;
    }
    bool operator!=(VertexWeightingMode const & o) const {
        return !(*this == o);
    }

    // ---- Named presets ----

    // Standard propagation: particle travels through geometry,
    // interaction probability computed from cross section * density.
    // This is the default and matches all existing SIREN behavior.
    static VertexWeightingMode Propagated() {
        return {true, true, BoundSource::Geometry};
    }

    // Fixed vertex: position externally determined (dk2nu, GENIE).
    // No interaction or position probability computed.
    static VertexWeightingMode Fixed() {
        return {false, false, BoundSource::None};
    }

    // External bounds: interaction probability computed but within
    // bounds provided by the distribution (e.g., gamma-ray segments).
    static VertexWeightingMode ExternalBounds() {
        return {true, true, BoundSource::Distribution};
    }
};

std::string VertexWeightingModeName(VertexWeightingMode const & mode);

} // namespace dataclasses
} // namespace siren

#endif // SIREN_VertexWeightingMode_H
