#pragma once
#ifndef SIREN_PhaseSpaceConvention_H
#define SIREN_PhaseSpaceConvention_H

#include <string>

namespace siren {
namespace dataclasses {

// ------------------------------------------------------------------ //
//  Topology: structural shape of the final state                      //
// ------------------------------------------------------------------ //

enum class PhaseSpaceTopology {
    Decay2Body,       // P -> A + B               (1 orientational DOF)
    Decay3Body,       // P -> A + B + C            (5 DOF)
    DecayNBody,       // P -> n particles, n >= 4  (3n-7 DOF)
    Scatter2to2,      // A + B -> C + D            (2 DOF)
    Scatter2to3,      // A + B -> C + D + E        (5 DOF)
    Unspecified       // opaque to the framework
};

std::string PhaseSpaceTopologyName(PhaseSpaceTopology topology);

// ------------------------------------------------------------------ //
//  Measure: parameterization the density is differential in           //
// ------------------------------------------------------------------ //

struct PhaseSpaceMeasure {
    enum class Type {
        SolidAngleRest,   // dOmega_rest
        SolidAngleLab,    // dOmega_lab
        Recursive2Body,   // ds_pair dOmega_pair dOmega_sub
        DalitzPair,       // ds_12 ds_23
        HelicityAngles,   // ds_pair dOmega_hel dOmega_sub
        MandelstamQ2,     // dQ^2
        BjorkenXY,        // dx dy
        Unspecified       // model-specific, no auto-conversion
    };

    Type type = Type::Unspecified;

    // 3-body factorization indices (secondary_momenta/masses indices).
    // Meaningful for Recursive2Body, DalitzPair, HelicityAngles.
    // Ignored for other types.
    int spectator = 0;
    int pair_first = 1;
    int pair_second = 2;

    bool operator==(PhaseSpaceMeasure const & o) const;
    bool operator!=(PhaseSpaceMeasure const & o) const { return !(*this == o); }

    // Convenience factories
    static PhaseSpaceMeasure SolidAngleRest();
    static PhaseSpaceMeasure SolidAngleLab();
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
    static PhaseSpaceMeasure BjorkenXY();
    static PhaseSpaceMeasure Unspecified();
};

std::string PhaseSpaceMeasureName(PhaseSpaceMeasure const & measure);

// ------------------------------------------------------------------ //
//  Convertibility                                                     //
// ------------------------------------------------------------------ //

int MeasureConvertibilityGroup(PhaseSpaceTopology topology,
                               PhaseSpaceMeasure const & measure);

bool PhaseSpaceCompatible(PhaseSpaceTopology topo_a, PhaseSpaceMeasure const & meas_a,
                          PhaseSpaceTopology topo_b, PhaseSpaceMeasure const & meas_b);

// ------------------------------------------------------------------ //
//  Legacy PhaseSpaceConvention (deprecated, kept for transition)       //
// ------------------------------------------------------------------ //

enum class PhaseSpaceConvention {
    RestFrameSolidAngle,
    LabFrameSolidAngle,
    Recursive2Body,
    Dalitz,
    HelicityAngles,
    BjorkenXY,
    MandelstamST,
    Custom
};

std::string PhaseSpaceConventionName(PhaseSpaceConvention convention);

PhaseSpaceTopology TopologyFromConvention(PhaseSpaceConvention convention,
                                          int n_secondaries);
PhaseSpaceMeasure MeasureFromConvention(PhaseSpaceConvention convention);

} // namespace dataclasses
} // namespace siren

#endif // SIREN_PhaseSpaceConvention_H
