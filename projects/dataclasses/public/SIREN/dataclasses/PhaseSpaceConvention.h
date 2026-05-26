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

enum class PhaseSpaceMeasure {
    SolidAngleRest,   // dOmega_rest
    SolidAngleLab,    // dOmega_lab
    Recursive2Body,   // ds_pair dOmega_pair dOmega_sub
    DalitzPair,       // ds_12 ds_23
    HelicityAngles,   // ds_pair dOmega_hel dOmega_sub
    MandelstamQ2,     // dQ^2
    BjorkenXY,        // dx dy
    Unspecified        // model-specific, no auto-conversion
};

std::string PhaseSpaceMeasureName(PhaseSpaceMeasure measure);

// ------------------------------------------------------------------ //
//  Convertibility                                                     //
// ------------------------------------------------------------------ //

// Returns the convertibility group index for a (topology, measure) pair.
// Measures in the same group within the same topology can be converted
// via analytic Jacobians.  Returns -1 for Unspecified or incompatible
// combinations.
int MeasureConvertibilityGroup(PhaseSpaceTopology topology,
                               PhaseSpaceMeasure measure);

// Returns true if two (topology, measure) pairs can be combined in a
// MultiChannelPhaseSpace (same topology, same convertibility group).
bool PhaseSpaceCompatible(PhaseSpaceTopology topo_a, PhaseSpaceMeasure meas_a,
                          PhaseSpaceTopology topo_b, PhaseSpaceMeasure meas_b);

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

// Convert legacy convention to the new (topology, measure) pair.
// Topology is inferred from convention + n_secondaries.
PhaseSpaceTopology TopologyFromConvention(PhaseSpaceConvention convention,
                                          int n_secondaries);
PhaseSpaceMeasure MeasureFromConvention(PhaseSpaceConvention convention);

} // namespace dataclasses
} // namespace siren

#endif // SIREN_PhaseSpaceConvention_H
