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
    if (!IndicesRelevant(type)) return true;
    return spectator == o.spectator
        && pair_first == o.pair_first
        && pair_second == o.pair_second;
}

PhaseSpaceMeasure PhaseSpaceMeasure::SolidAngleRest() {
    return {Type::SolidAngleRest, 0, 1, 2};
}
PhaseSpaceMeasure PhaseSpaceMeasure::SolidAngleLab() {
    return {Type::SolidAngleLab, 0, 1, 2};
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
PhaseSpaceMeasure PhaseSpaceMeasure::BjorkenXY() {
    return {Type::BjorkenXY, 0, 1, 2};
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
        case PhaseSpaceMeasure::Type::BjorkenXY:       return "BjorkenXY";
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
        switch (measure.type) {
            case T::SolidAngleRest:
            case T::SolidAngleLab:
            case T::MandelstamQ2:
            case T::BjorkenXY:
                return 0;
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

bool PhaseSpaceCompatible(PhaseSpaceTopology topo_a, PhaseSpaceMeasure const & meas_a,
                          PhaseSpaceTopology topo_b, PhaseSpaceMeasure const & meas_b)
{
    if (topo_a != topo_b) return false;
    if (topo_a == PhaseSpaceTopology::Unspecified) {
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
            return PhaseSpaceMeasure::SolidAngleRest();
        case PhaseSpaceConvention::LabFrameSolidAngle:
            return PhaseSpaceMeasure::SolidAngleLab();
        case PhaseSpaceConvention::Recursive2Body:
            return PhaseSpaceMeasure::Recursive2Body();
        case PhaseSpaceConvention::Dalitz:
            return PhaseSpaceMeasure::DalitzPair();
        case PhaseSpaceConvention::HelicityAngles:
            return PhaseSpaceMeasure::HelicityAngles();
        case PhaseSpaceConvention::BjorkenXY:
            return PhaseSpaceMeasure::BjorkenXY();
        case PhaseSpaceConvention::MandelstamST:
            return PhaseSpaceMeasure::MandelstamQ2();
        case PhaseSpaceConvention::Custom:
            return PhaseSpaceMeasure::Unspecified();
    }
    return PhaseSpaceMeasure::Unspecified();
}

} // namespace dataclasses
} // namespace siren
