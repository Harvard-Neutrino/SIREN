#include "LeptonInjector/utilities/Interpolator.h"

#include <tuple>
#include <vector>

using namespace LI::utilities;

template<>
template<>
bool TableData1D<double>::operator==<double>(TableData1D<double> const & other) const {
    return std::tie(x, f) == std::tie(other.x, other.f);
}

template<>
template<>
bool TableData2D<double>::operator==<double>(TableData2D<double> const & other) const {
    return std::tie(x, y, f) == std::tie(other.x, other.y, other.f);
}

template<>
template<>
bool Interpolator1D<double>::operator==<double>(Interpolator1D<double> const & other) const {
    return original_table == other.original_table;
}

template<>
template<>
bool Interpolator2D<double>::operator==<double>(Interpolator2D<double> const & other) const {
    return original_table == other.original_table;
}

