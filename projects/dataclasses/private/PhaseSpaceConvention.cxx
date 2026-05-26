#include "SIREN/dataclasses/PhaseSpaceConvention.h"

namespace siren {
namespace dataclasses {

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

} // namespace dataclasses
} // namespace siren
