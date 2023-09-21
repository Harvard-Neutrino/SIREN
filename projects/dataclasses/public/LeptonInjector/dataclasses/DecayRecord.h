#ifndef LI_DecayRecord_H
#define LI_DecayRecord_H

#include <array>                                        // for array
#include <vector>                                       // for vector
#include <cstdint>                                      // for uint32_t
#include <ostream>                                      // for ostream
#include <stdexcept>                                    // for runtime_error

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include "LeptonInjector/dataclasses/DecaySignature.h"

namespace LI {
namespace dataclasses {

struct DecayRecord {
    DecaySignature signature;
    double primary_mass;
    std::array<double, 4> primary_momentum;
    double primary_helicity;
    std::array<double, 3> decay_vertex = {0, 0, 0};
    std::vector<double> secondary_masses;
    std::vector<std::array<double, 4>> secondary_momenta;
    std::vector<double> secondary_helicity;
    std::vector<double> decay_parameters;
    bool operator==(DecayRecord const & other) const;
    friend std::ostream& operator<<(std::ostream& os, DecayRecord const& record);
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(::cereal::make_nvp("DecaySignature", signature));
            archive(::cereal::make_nvp("PrimaryMass", primary_mass));
            archive(::cereal::make_nvp("PrimaryMomentum", primary_momentum));
            archive(::cereal::make_nvp("PrimaryHelicity", primary_helicity));
            archive(::cereal::make_nvp("DecayVertex", decay_vertex));
            archive(::cereal::make_nvp("SecondaryMasses", secondary_masses));
            archive(::cereal::make_nvp("SecondaryMomenta", secondary_momenta));
            archive(::cereal::make_nvp("SecondaryHelicity", secondary_helicity));
            archive(::cereal::make_nvp("DecayParameters", decay_parameters));
        } else {
            throw std::runtime_error("DecayRecord only supports version <= 0!");
        }
    };
};

} // namespace dataclasses
} // namespace LI

CEREAL_CLASS_VERSION(LI::dataclasses::DecayRecord, 0);

#endif // LI_DecayRecord_H
