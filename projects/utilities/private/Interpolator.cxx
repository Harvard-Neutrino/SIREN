#include "SIREN/utilities/Interpolator.h"

#include <tuple>
#include <vector>

// Explicit template instantiation
template struct siren::utilities::TableData1D<double>;
template struct siren::utilities::TableData2D<double>;
template struct siren::utilities::IndexFinderRegular<double>;
template struct siren::utilities::IndexFinderIrregular<double>;
template struct siren::utilities::Indexer1D<double>;
template struct siren::utilities::Interpolator1D<double>;
template struct siren::utilities::Interpolator2D<double>;

using namespace siren::utilities;

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

CEREAL_CLASS_VERSION(siren::utilities::TableData1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::TableData2D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::IndexFinderRegular<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::IndexFinderIrregular<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Indexer1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Interpolator1D<double>, 0);
CEREAL_CLASS_VERSION(siren::utilities::Interpolator2D<double>, 0);

CEREAL_REGISTER_DYNAMIC_INIT(siren_Interpolator);

