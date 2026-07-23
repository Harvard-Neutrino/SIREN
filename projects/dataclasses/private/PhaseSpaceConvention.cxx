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
//  PhaseSpaceMeasure                                                  //
// ------------------------------------------------------------------ //

namespace {
bool IndicesRelevant(PhaseSpaceMeasure::Type t) {
    return t == PhaseSpaceMeasure::Type::Recursive2Body
        || t == PhaseSpaceMeasure::Type::DalitzPair
        || t == PhaseSpaceMeasure::Type::HelicityAngles;
}
} // anonymous namespace

bool PhaseSpaceMeasure::operator==(PhaseSpaceMeasure const & o) const {
    if (type != o.type) return false;
    if (type == Type::SolidAngleLab) return spectator == o.spectator;
    if (!IndicesRelevant(type)) return true;
    return spectator == o.spectator
        && pair_first == o.pair_first
        && pair_second == o.pair_second;
}

PhaseSpaceMeasure PhaseSpaceMeasure::SolidAngleRest() {
    return {Type::SolidAngleRest, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::SolidAngleLab(int daughter_index) {
    return {Type::SolidAngleLab, daughter_index, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::Recursive2Body(int s, int f, int sec) {
    return {Type::Recursive2Body, s, f, sec};
}
PhaseSpaceMeasure PhaseSpaceMeasure::DalitzPair(int s, int f, int sec) {
    return {Type::DalitzPair, s, f, sec};
}
PhaseSpaceMeasure PhaseSpaceMeasure::HelicityAngles(int s, int f, int sec) {
    return {Type::HelicityAngles, s, f, sec};
}
PhaseSpaceMeasure PhaseSpaceMeasure::MandelstamQ2() {
    return {Type::MandelstamQ2, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::MandelstamQ2Phi() {
    return {Type::MandelstamQ2Phi, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::FixedMassY() {
    return {Type::FixedMassY, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::FixedMassYPhi() {
    return {Type::FixedMassYPhi, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::BjorkenXY() {
    return {Type::BjorkenXY, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::BjorkenXYPhi() {
    return {Type::BjorkenXYPhi, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::MandelstamQ2Y() {
    return {Type::MandelstamQ2Y, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::MandelstamQ2YPhi() {
    return {Type::MandelstamQ2YPhi, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::Unspecified() {
    return {Type::Unspecified, 0, 1, 2};
}

std::string PhaseSpaceMeasureName(PhaseSpaceMeasure const & measure) {
    switch (measure.type) {
        case PhaseSpaceMeasure::Type::SolidAngleRest:  return "SolidAngleRest";
        case PhaseSpaceMeasure::Type::SolidAngleLab:   return "SolidAngleLab";
        case PhaseSpaceMeasure::Type::Recursive2Body:  return "Recursive2Body";
        case PhaseSpaceMeasure::Type::DalitzPair:      return "DalitzPair";
        case PhaseSpaceMeasure::Type::HelicityAngles:  return "HelicityAngles";
        case PhaseSpaceMeasure::Type::MandelstamQ2:    return "MandelstamQ2";
        case PhaseSpaceMeasure::Type::MandelstamQ2Phi: return "MandelstamQ2Phi";
        case PhaseSpaceMeasure::Type::FixedMassY:      return "FixedMassY";
        case PhaseSpaceMeasure::Type::FixedMassYPhi:   return "FixedMassYPhi";
        case PhaseSpaceMeasure::Type::BjorkenXY:       return "BjorkenXY";
        case PhaseSpaceMeasure::Type::BjorkenXYPhi:    return "BjorkenXYPhi";
        case PhaseSpaceMeasure::Type::MandelstamQ2Y:   return "MandelstamQ2Y";
        case PhaseSpaceMeasure::Type::MandelstamQ2YPhi: return "MandelstamQ2YPhi";
        case PhaseSpaceMeasure::Type::Unspecified:     return "Unspecified";
    }
    return "Unknown";
}

// ------------------------------------------------------------------ //
//  Convertibility groups                                              //
// ------------------------------------------------------------------ //

int MeasureConvertibilityGroup(PhaseSpaceTopology topology,
                               PhaseSpaceMeasure const & measure)
{
    using T = PhaseSpaceMeasure::Type;
    if (measure.type == T::Unspecified) return -1;

    switch (topology) {

    case PhaseSpaceTopology::Decay2Body:
        switch (measure.type) {
            case T::SolidAngleRest:
            case T::SolidAngleLab:
                return 0;
            default: return -1;
        }

    case PhaseSpaceTopology::Decay3Body:
        switch (measure.type) {
            case T::Recursive2Body:
            case T::DalitzPair:
            case T::HelicityAngles:
                return 0;
            case T::SolidAngleRest:
                return 1;
            default: return -1;
        }

    case PhaseSpaceTopology::DecayNBody:
        switch (measure.type) {
            case T::SolidAngleRest:
                return 0;
            default: return -1;
        }

    case PhaseSpaceTopology::Scatter2to2:
        // SolidAngleLab is excluded: ConvertDensity implements the
        // rest<->lab boost only for Decay2Body, so listing it here would
        // promise an auto-conversion that throws at evaluation time.
        switch (measure.type) {
            case T::SolidAngleRest:
            case T::MandelstamQ2:
            case T::MandelstamQ2Phi:
            case T::FixedMassY:
            case T::FixedMassYPhi:
                return 0;
            case T::BjorkenXY:
            case T::BjorkenXYPhi:
            case T::MandelstamQ2Y:
            case T::MandelstamQ2YPhi:
                return 1;
            default: return -1;
        }

    case PhaseSpaceTopology::Scatter2to3:
        switch (measure.type) {
            case T::Recursive2Body:
            case T::DalitzPair:
                return 0;
            case T::SolidAngleRest:
                return 1;
            default: return -1;
        }

    case PhaseSpaceTopology::Unspecified:
        return -1;
    }
    return -1;
}

bool MeasureHasExplicitAzimuth(PhaseSpaceMeasure const & measure) {
    using T = PhaseSpaceMeasure::Type;
    switch (measure.type) {
        case T::SolidAngleRest:
        case T::MandelstamQ2Phi:
        case T::FixedMassYPhi:
        case T::BjorkenXYPhi:
        case T::MandelstamQ2YPhi:
            return true;
        default:
            return false;
    }
}

bool MeasureIntegratesAzimuth(PhaseSpaceMeasure const & measure) {
    using T = PhaseSpaceMeasure::Type;
    switch (measure.type) {
        case T::MandelstamQ2:
        case T::FixedMassY:
        case T::BjorkenXY:
        case T::MandelstamQ2Y:
            return true;
        default:
            return false;
    }
}

PhaseSpaceMeasure MeasureWithExplicitAzimuth(PhaseSpaceMeasure const & measure) {
    using T = PhaseSpaceMeasure::Type;
    switch (measure.type) {
        case T::MandelstamQ2:  return PhaseSpaceMeasure::MandelstamQ2Phi();
        case T::FixedMassY:    return PhaseSpaceMeasure::FixedMassYPhi();
        case T::BjorkenXY:     return PhaseSpaceMeasure::BjorkenXYPhi();
        case T::MandelstamQ2Y: return PhaseSpaceMeasure::MandelstamQ2YPhi();
        default:               return measure;
    }
}

bool PhaseSpaceDensityConvertible(PhaseSpaceTopology topology,
                                  PhaseSpaceMeasure const & from,
                                  PhaseSpaceMeasure const & to)
{
    if (from == to) return true;

    int from_family = MeasureConvertibilityGroup(topology, from);
    int to_family = MeasureConvertibilityGroup(topology, to);
    if (from_family < 0 || from_family != to_family) return false;

    // A pointwise conversion can lift a declared-uniform azimuth but never
    // integrate an explicit one. Only the 2->2 families mix the two forms.
    if (MeasureHasExplicitAzimuth(from) && MeasureIntegratesAzimuth(to)) {
        return false;
    }
    return true;
}

} // namespace dataclasses
} // namespace siren
