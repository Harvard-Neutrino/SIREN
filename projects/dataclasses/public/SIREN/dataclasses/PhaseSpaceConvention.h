#pragma once
#ifndef SIREN_PhaseSpaceConvention_H
#define SIREN_PhaseSpaceConvention_H

#include <cstdint>
#include <stdexcept>
#include <string>

#include <cereal/cereal.hpp>

namespace siren {
namespace dataclasses {

// ------------------------------------------------------------------ //
//  Topology: structural shape of the final state                      //
// ------------------------------------------------------------------ //

enum class PhaseSpaceTopology {
    Decay2Body,       // P -> A + B                (2 orientational DOF)
    Decay3Body,       // P -> A + B + C            (5 DOF)
    DecayNBody,       // P -> n particles, n >= 4  (3n-4 DOF)
    Scatter2to2,      // A + B -> C + D            (2 DOF, assuming fixed final-state masses)
    Scatter2to3,      // A + B -> C + D + E        (5 DOF, assuming fixed final-state masses)
    Unspecified       // opaque to the framework
};

std::string PhaseSpaceTopologyName(PhaseSpaceTopology topology);

// ------------------------------------------------------------------ //
//  Measure: parameterization the density is differential in           //
// ------------------------------------------------------------------ //

struct PhaseSpaceMeasure {
    // The three-body measures (Recursive2Body, DalitzPair, HelicityAngles)
    // are joint densities over the full three-body chart: conversions
    // between them exchange coordinates pointwise and never integrate an
    // orientation. A model whose orientation is physically uniform includes
    // the uniform factors (for example 1/(4*pi) per isotropic dOmega) in
    // its density values.
    //
    // The 2->2 scattering measures come in azimuth-integrated and
    // explicit-azimuth pairs. An integrated measure declares its omitted
    // beam-axis azimuth uniform, so it lifts to the explicit form by
    // 1/(2*pi); an omitted azimuth that is opaque or nonuniform must use
    // Unspecified instead.
    enum class Type {
        SolidAngleRest,   // dOmega_rest
        SolidAngleLab,    // dOmega_lab
        Recursive2Body,   // ds_pair dOmega_pair dOmega_sub
        DalitzPair,       // ds_12 ds_23 (orientation coordinates shared)
        HelicityAngles,   // ds_pair dOmega_hel dOmega_sub
        MandelstamQ2,     // dQ^2, azimuth integrated with a uniform conditional
        MandelstamQ2Phi,  // dQ^2 dphi
        FixedMassY,       // dy, azimuth integrated with a uniform conditional
        FixedMassYPhi,    // dy dphi for fixed-mass 2->2 scattering
        BjorkenXY,        // dx dy, azimuth integrated with a uniform conditional
        BjorkenXYPhi,     // dx dy dphi
        MandelstamQ2Y,    // dQ^2 dy, azimuth integrated with a uniform conditional
        MandelstamQ2YPhi, // dQ^2 dy dphi
        Unspecified       // model-specific, no auto-conversion
    };

    Type type = Type::Unspecified;

    // Secondary-momenta/masses indices. For SolidAngleLab, `spectator`
    // identifies the daughter whose lab angle defines the measure. All three
    // fields describe the factorization for Recursive2Body, DalitzPair, and
    // HelicityAngles. Ignored for other types.
    int spectator = 0;
    int pair_first = 1;
    int pair_second = 2;

    bool operator==(PhaseSpaceMeasure const & o) const;
    bool operator!=(PhaseSpaceMeasure const & o) const { return !(*this == o); }

    template<class Archive>
    void save(Archive & archive, std::uint32_t const version) const {
        if (version == 0) {
            archive(::cereal::make_nvp("Type", static_cast<int>(type)));
            archive(::cereal::make_nvp("Spectator", spectator));
            archive(::cereal::make_nvp("PairFirst", pair_first));
            archive(::cereal::make_nvp("PairSecond", pair_second));
        } else {
            throw std::runtime_error(
                "PhaseSpaceMeasure only supports version <= 0!");
        }
    }

    template<class Archive>
    void load(Archive & archive, std::uint32_t const version) {
        if (version == 0) {
            int type_int;
            int loaded_spectator;
            int loaded_pair_first;
            int loaded_pair_second;
            archive(::cereal::make_nvp("Type", type_int));
            archive(::cereal::make_nvp("Spectator", loaded_spectator));
            archive(::cereal::make_nvp("PairFirst", loaded_pair_first));
            archive(::cereal::make_nvp("PairSecond", loaded_pair_second));

            switch (static_cast<Type>(type_int)) {
                case Type::SolidAngleRest:
                case Type::SolidAngleLab:
                case Type::Recursive2Body:
                case Type::DalitzPair:
                case Type::HelicityAngles:
                case Type::MandelstamQ2:
                case Type::MandelstamQ2Phi:
                case Type::FixedMassY:
                case Type::FixedMassYPhi:
                case Type::BjorkenXY:
                case Type::BjorkenXYPhi:
                case Type::MandelstamQ2Y:
                case Type::MandelstamQ2YPhi:
                case Type::Unspecified:
                    break;
                default:
                    throw std::runtime_error(
                        "PhaseSpaceMeasure: invalid Type value "
                        + std::to_string(type_int) + " in archive");
            }

            type = static_cast<Type>(type_int);
            spectator = loaded_spectator;
            pair_first = loaded_pair_first;
            pair_second = loaded_pair_second;
        } else {
            throw std::runtime_error(
                "PhaseSpaceMeasure only supports version <= 0!");
        }
    }

    // Convenience factories
    static PhaseSpaceMeasure SolidAngleRest();
    static PhaseSpaceMeasure SolidAngleLab(int daughter_index = 0);
    static PhaseSpaceMeasure Recursive2Body(int spectator = 0,
                                            int pair_first = 1,
                                            int pair_second = 2);
    static PhaseSpaceMeasure DalitzPair(int spectator = 0,
                                        int pair_first = 1,
                                        int pair_second = 2);
    static PhaseSpaceMeasure HelicityAngles(int spectator = 0,
                                            int pair_first = 1,
                                            int pair_second = 2);
    static PhaseSpaceMeasure MandelstamQ2();
    static PhaseSpaceMeasure MandelstamQ2Phi();
    static PhaseSpaceMeasure FixedMassY();
    static PhaseSpaceMeasure FixedMassYPhi();
    static PhaseSpaceMeasure BjorkenXY();
    static PhaseSpaceMeasure BjorkenXYPhi();
    static PhaseSpaceMeasure MandelstamQ2Y();
    static PhaseSpaceMeasure MandelstamQ2YPhi();
    static PhaseSpaceMeasure Unspecified();
};

std::string PhaseSpaceMeasureName(PhaseSpaceMeasure const & measure);

// ------------------------------------------------------------------ //
//  Convention: a topology plus a measure                              //
// ------------------------------------------------------------------ //

// A complete phase-space convention. A density value is only meaningful
// relative to one of these; the weighter resolves a single convention per
// vertex and evaluates both sides of the weight ratio in it.
struct PhaseSpaceConvention {
    PhaseSpaceTopology topology = PhaseSpaceTopology::Unspecified;
    PhaseSpaceMeasure measure = PhaseSpaceMeasure::Unspecified();

    bool operator==(PhaseSpaceConvention const & o) const {
        return topology == o.topology && measure == o.measure;
    }
    bool operator!=(PhaseSpaceConvention const & o) const {
        return !(*this == o);
    }
};

// ------------------------------------------------------------------ //
//  Azimuth taxonomy for the 2->2 scattering measures                  //
// ------------------------------------------------------------------ //

// True when the measure carries the beam-axis azimuth as an explicit
// coordinate (SolidAngleRest and the *Phi scattering measures).
bool MeasureHasExplicitAzimuth(PhaseSpaceMeasure const & measure);

// True when the measure integrates the beam-axis azimuth with a declared
// uniform conditional (MandelstamQ2, FixedMassY, BjorkenXY, MandelstamQ2Y).
bool MeasureIntegratesAzimuth(PhaseSpaceMeasure const & measure);

// The declared explicit-azimuth completion of an azimuth-integrated
// measure (MandelstamQ2 -> MandelstamQ2Phi and so on). Identity for every
// other measure.
PhaseSpaceMeasure MeasureWithExplicitAzimuth(PhaseSpaceMeasure const & measure);

// ------------------------------------------------------------------ //
//  Convertibility                                                     //
// ------------------------------------------------------------------ //

int MeasureConvertibilityGroup(PhaseSpaceTopology topology,
                               PhaseSpaceMeasure const & measure);

// Return whether a density in `from` can be evaluated pointwise in `to`.
// This is intentionally directional: an azimuth-integrated density has a
// declared uniform lift to an explicit-azimuth density, while a general joint
// density cannot be marginalized by a pointwise conversion.
bool PhaseSpaceDensityConvertible(PhaseSpaceTopology topology,
                                  PhaseSpaceMeasure const & from,
                                  PhaseSpaceMeasure const & to);

} // namespace dataclasses
} // namespace siren

CEREAL_CLASS_VERSION(siren::dataclasses::PhaseSpaceMeasure, 0);

#endif // SIREN_PhaseSpaceConvention_H
