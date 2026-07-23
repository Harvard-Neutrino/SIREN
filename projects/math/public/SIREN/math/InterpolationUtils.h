#pragma once
#ifndef SIREN_InterpolationUtils_H
#define SIREN_InterpolationUtils_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace siren {
namespace math {

// Return the adjacent nodes that bracket x in an ascending grid, clamping to
// the first/last interval. Exact interior nodes can belong to either adjacent
// bin so callers can preserve their interpolation convention.
template<typename T>
inline std::pair<std::size_t, std::size_t> InterpolationBracket(
    std::vector<T> const & grid,
    T const & x,
    bool right_bin_at_node = true)
{
    if (grid.size() < 2) {
        throw std::runtime_error(
            "InterpolationBracket requires at least two grid points");
    }
    if (x <= grid.front()) return {0, 1};
    if (x >= grid.back()) return {grid.size() - 2, grid.size() - 1};
    auto upper = right_bin_at_node
        ? std::upper_bound(grid.begin(), grid.end(), x)
        : std::lower_bound(grid.begin(), grid.end(), x);
    std::size_t hi = static_cast<std::size_t>(
        std::distance(grid.begin(), upper));
    return {hi - 1, hi};
}

} // namespace math
} // namespace siren

#endif // SIREN_InterpolationUtils_H
