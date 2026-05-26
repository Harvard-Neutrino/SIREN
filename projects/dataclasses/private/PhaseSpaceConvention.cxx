#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren {
namespace dataclasses {

// ------------------------------------------------------------------ //
//  Topology names                                                     //
// ------------------------------------------------------------------ //

std::string PhaseSpaceTopologyName(PhaseSpaceTopology topology) {
    switch (topology) {
        case PhaseSpaceTopology::Decay2Body:   return "Decay2Body";
        case PhaseSpaceTopology::Decay3Body:   return "Decay3Body";
        case PhaseSpaceTopology::DecayNBody:   return "DecayNBody";
        case PhaseSpaceTopology::Scatter2to2:  return "Scatter2to2";
        case PhaseSpaceTopology::Scatter2to3:  return "Scatter2to3";
        case PhaseSpaceTopology::Unspecified:  return "Unspecified";
    }
    return "Unknown";
}

// ------------------------------------------------------------------ //
//  Measure names                                                      //
// ------------------------------------------------------------------ //

std::string PhaseSpaceMeasureName(PhaseSpaceMeasure measure) {
    switch (measure) {
        case PhaseSpaceMeasure::SolidAngleRest:  return "SolidAngleRest";
        case PhaseSpaceMeasure::SolidAngleLab:   return "SolidAngleLab";
        case PhaseSpaceMeasure::Recursive2Body:  return "Recursive2Body";
        case PhaseSpaceMeasure::DalitzPair:      return "DalitzPair";
        case PhaseSpaceMeasure::HelicityAngles:  return "HelicityAngles";
        case PhaseSpaceMeasure::MandelstamQ2:    return "MandelstamQ2";
        case PhaseSpaceMeasure::BjorkenXY:       return "BjorkenXY";
        case PhaseSpaceMeasure::Unspecified:     return "Unspecified";
    }
    return "Unknown";
}

// ------------------------------------------------------------------ //
//  Convertibility groups                                              //
// ------------------------------------------------------------------ //
//
// Within each topology, measures in the same group can be converted
// via analytic Jacobians.  The group numbering is per-topology.
//
//   Decay2Body:
//     Group 0: SolidAngleRest, SolidAngleLab
//
//   Decay3Body:
//     Group 0: Recursive2Body, DalitzPair, HelicityAngles
//     Group 1: SolidAngleRest  (orientation-separated, not convertible to Group 0)
//
//   DecayNBody:
//     Group 0: SolidAngleRest
//
//   Scatter2to2:
//     Group 0: SolidAngleRest, SolidAngleLab, MandelstamQ2, BjorkenXY
//              (all fully interconvertible)
//
//   Scatter2to3:
//     Group 0: Recursive2Body, DalitzPair
//     Group 1: SolidAngleRest
//

int MeasureConvertibilityGroup(PhaseSpaceTopology topology,
                               PhaseSpaceMeasure measure)
{
    if (measure == PhaseSpaceMeasure::Unspecified) return -1;

    switch (topology) {

    case PhaseSpaceTopology::Decay2Body:
        switch (measure) {
            case PhaseSpaceMeasure::SolidAngleRest:
            case PhaseSpaceMeasure::SolidAngleLab:
                return 0;
            default: return -1;
        }

    case PhaseSpaceTopology::Decay3Body:
        switch (measure) {
            case PhaseSpaceMeasure::Recursive2Body:
            case PhaseSpaceMeasure::DalitzPair:
            case PhaseSpaceMeasure::HelicityAngles:
                return 0;
            case PhaseSpaceMeasure::SolidAngleRest:
                return 1;
            default: return -1;
        }

    case PhaseSpaceTopology::DecayNBody:
        switch (measure) {
            case PhaseSpaceMeasure::SolidAngleRest:
                return 0;
            default: return -1;
        }

    case PhaseSpaceTopology::Scatter2to2:
        switch (measure) {
            case PhaseSpaceMeasure::SolidAngleRest:
            case PhaseSpaceMeasure::SolidAngleLab:
            case PhaseSpaceMeasure::MandelstamQ2:
            case PhaseSpaceMeasure::BjorkenXY:
                return 0;
            default: return -1;
        }

    case PhaseSpaceTopology::Scatter2to3:
        switch (measure) {
            case PhaseSpaceMeasure::Recursive2Body:
            case PhaseSpaceMeasure::DalitzPair:
                return 0;
            case PhaseSpaceMeasure::SolidAngleRest:
                return 1;
            default: return -1;
        }

    case PhaseSpaceTopology::Unspecified:
        return -1;
    }
    return -1;
}

bool PhaseSpaceCompatible(PhaseSpaceTopology topo_a, PhaseSpaceMeasure meas_a,
                          PhaseSpaceTopology topo_b, PhaseSpaceMeasure meas_b)
{
    if (topo_a != topo_b) return false;
    if (topo_a == PhaseSpaceTopology::Unspecified) {
        // Unspecified topology: only compatible if measures also match
        return meas_a == meas_b;
    }
    if (meas_a == meas_b) return true;
    int ga = MeasureConvertibilityGroup(topo_a, meas_a);
    int gb = MeasureConvertibilityGroup(topo_b, meas_b);
    if (ga < 0 || gb < 0) return false;
    return ga == gb;
}

// ------------------------------------------------------------------ //
//  Legacy PhaseSpaceConvention                                        //
// ------------------------------------------------------------------ //

std::string PhaseSpaceConventionName(PhaseSpaceConvention convention) {
    switch (convention) {
        case PhaseSpaceConvention::RestFrameSolidAngle:
            return "RestFrameSolidAngle";
        case PhaseSpaceConvention::LabFrameSolidAngle:
            return "LabFrameSolidAngle";
        case PhaseSpaceConvention::Recursive2Body:
            return "Recursive2Body";
        case PhaseSpaceConvention::Dalitz:
            return "Dalitz";
        case PhaseSpaceConvention::HelicityAngles:
            return "HelicityAngles";
        case PhaseSpaceConvention::BjorkenXY:
            return "BjorkenXY";
        case PhaseSpaceConvention::MandelstamST:
            return "MandelstamST";
        case PhaseSpaceConvention::Custom:
            return "Custom";
    }
    return "Unknown";
}

PhaseSpaceTopology TopologyFromConvention(PhaseSpaceConvention convention,
                                          int n_secondaries)
{
    switch (convention) {
        case PhaseSpaceConvention::RestFrameSolidAngle:
        case PhaseSpaceConvention::LabFrameSolidAngle:
            if (n_secondaries == 2) return PhaseSpaceTopology::Decay2Body;
            if (n_secondaries == 3) return PhaseSpaceTopology::Decay3Body;
            if (n_secondaries > 3)  return PhaseSpaceTopology::DecayNBody;
            return PhaseSpaceTopology::Unspecified;
        case PhaseSpaceConvention::Recursive2Body:
        case PhaseSpaceConvention::Dalitz:
        case PhaseSpaceConvention::HelicityAngles:
            return PhaseSpaceTopology::Decay3Body;
        case PhaseSpaceConvention::BjorkenXY:
        case PhaseSpaceConvention::MandelstamST:
            return PhaseSpaceTopology::Scatter2to2;
        case PhaseSpaceConvention::Custom:
            return PhaseSpaceTopology::Unspecified;
    }
    return PhaseSpaceTopology::Unspecified;
}

PhaseSpaceMeasure MeasureFromConvention(PhaseSpaceConvention convention)
{
    switch (convention) {
        case PhaseSpaceConvention::RestFrameSolidAngle:
            return PhaseSpaceMeasure::SolidAngleRest;
        case PhaseSpaceConvention::LabFrameSolidAngle:
            return PhaseSpaceMeasure::SolidAngleLab;
        case PhaseSpaceConvention::Recursive2Body:
            return PhaseSpaceMeasure::Recursive2Body;
        case PhaseSpaceConvention::Dalitz:
            return PhaseSpaceMeasure::DalitzPair;
        case PhaseSpaceConvention::HelicityAngles:
            return PhaseSpaceMeasure::HelicityAngles;
        case PhaseSpaceConvention::BjorkenXY:
            return PhaseSpaceMeasure::BjorkenXY;
        case PhaseSpaceConvention::MandelstamST:
            return PhaseSpaceMeasure::MandelstamQ2;
        case PhaseSpaceConvention::Custom:
            return PhaseSpaceMeasure::Unspecified;
    }
    return PhaseSpaceMeasure::Unspecified;
}

} // namespace dataclasses
} // namespace siren
