#pragma once
#ifndef SIREN_PhaseSpaceConvention_H
#define SIREN_PhaseSpaceConvention_H

#include <string>

namespace siren {
namespace dataclasses {

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

} // namespace dataclasses
} // namespace siren

#endif // SIREN_PhaseSpaceConvention_H
