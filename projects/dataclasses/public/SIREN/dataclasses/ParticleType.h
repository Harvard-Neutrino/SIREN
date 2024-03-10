#pragma once
#ifndef SIREN_ParticleType_H
#define SIREN_ParticleType_H

#include <map>
#include <string>
#include <iostream>
#include <inttypes.h>

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>

#include <cereal/types/set.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>

namespace siren {
namespace dataclasses {

enum class ParticleType : int32_t {
#define X(a, b) a = b,
#include "SIREN/dataclasses/ParticleTypes.def"
#undef X
};

static const std::map<ParticleType, std::string> ParticleTypeNames = {
#define X(a, b) {ParticleType:: a , #a },
#include "SIREN/dataclasses/ParticleTypes.def"
#undef X
};

static const std::map<std::string, ParticleType> ParticleTypeMap = {
#define X(a, b) { #a , ParticleType:: a },
#include "SIREN/dataclasses/ParticleTypes.def"
#undef X
};

static const std::map<int32_t, std::string> ParticleTypeIntNames = {
#define X(a, b) { b , #a },
#include "SIREN/dataclasses/ParticleTypes.def"
#undef X
};

} // namespace siren
} // namespace dataclasses

std::ostream & operator<<(std::ostream & os, siren::dataclasses::ParticleType && p);
std::ostream & operator<<(std::ostream & os, siren::dataclasses::ParticleType const & p);

namespace cereal {
    template <class Archive>
    int save_minimal(Archive const &, siren::dataclasses::ParticleType const & p) {
        return static_cast<int>(p);
    }

    template <class Archive>
    void load_minimal(Archive const &, siren::dataclasses::ParticleType & p, int const & value )
    {
        p = static_cast<siren::dataclasses::ParticleType>(value);
    }

    template <class Archive, class C, class A> inline
    void save( Archive & ar, std::set<siren::dataclasses::ParticleType, C, A> const & set )
    {
      ar( make_size_tag( static_cast<size_type>(set.size()) ) );

      for( const auto & i : set )
        ar( i );
    }

    template <class Archive, class C, class A> inline
    void load( Archive & ar, std::set<siren::dataclasses::ParticleType, C, A> & set )
    {
      size_type size;
      ar( make_size_tag( size ) );

      set.clear();

      auto hint = set.begin();
      for( size_type i = 0; i < size; ++i )
      {
        typename std::set<siren::dataclasses::ParticleType>::key_type key;
        ar( key );
        hint = set.emplace_hint( hint, std::move( key ) );
      }
    }
}

#endif // SIREN_ParticleType_H
