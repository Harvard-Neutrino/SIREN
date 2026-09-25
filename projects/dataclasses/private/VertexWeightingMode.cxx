#include "SIREN/dataclasses/VertexWeightingMode.h"

namespace siren {
namespace dataclasses {

std::string VertexWeightingModeName(VertexWeightingMode const & mode) {
    if (!mode.compute_interaction_probability &&
        !mode.compute_position_probability &&
        mode.bound_source == VertexWeightingMode::BoundSource::None) {
        return "Fixed";
    }
    if (mode.compute_interaction_probability &&
        mode.compute_position_probability) {
        if (mode.bound_source == VertexWeightingMode::BoundSource::Geometry) {
            return "Propagated";
        }
        if (mode.bound_source == VertexWeightingMode::BoundSource::Distribution) {
            return "ExternalBounds";
        }
    }
    return "Custom";
}

} // namespace dataclasses
} // namespace siren
