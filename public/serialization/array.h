#pragma once
#ifndef LI_serialization_array_H
#define LI_serialization_array_H

#include "cereal/cereal.hpp"
#include <array>

namespace cereal {

//! Saving for std::array all other types
template <class Archive, size_t N> inline
typename std::enable_if<!traits::is_output_serializable<BinaryData<double>, Archive>::value
	|| !std::is_arithmetic<double>::value, void>::type
CEREAL_SAVE_FUNCTION_NAME( Archive & ar, std::array<double, N> const & array )
{
    ar( make_size_tag( static_cast<size_type>(N) ) ); // number of elements
	for(auto const & i : array)
		ar( i );
}

template <class Archive, size_t N> inline
typename std::enable_if<!traits::is_output_serializable<BinaryData<int>, Archive>::value
	|| !std::is_arithmetic<int>::value, void>::type
CEREAL_SAVE_FUNCTION_NAME( Archive & ar, std::array<int, N> const & array )
{
	size_type size;
    ar( make_size_tag( size ) );
	for(auto const & i : array)
		ar( i );
}

//! Loading for std::array all other types
template <class Archive, size_t N> inline
typename std::enable_if<!traits::is_input_serializable<BinaryData<double>, Archive>::value
	|| !std::is_arithmetic<double>::value, void>::type
CEREAL_LOAD_FUNCTION_NAME( Archive & ar, std::array<double, N> & array )
{
    ar( make_size_tag( static_cast<size_type>(N) ) ); // number of elements
	for(auto & i : array)
		ar(i);
}

template <class Archive, size_t N> inline
typename std::enable_if<!traits::is_input_serializable<BinaryData<int>, Archive>::value
	|| !std::is_arithmetic<int>::value, void>::type
CEREAL_LOAD_FUNCTION_NAME( Archive & ar, std::array<int, N> & array )
{
	size_type size;
    ar( make_size_tag( size ) );
	for(auto & i : array)
		ar(i);
}

} // namespace cereal

#endif // LI_serialization_array_H

