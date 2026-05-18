#include "GDMLParserPrivate.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

#include "SIREN/utilities/Constants.h"

using namespace siren::math;

namespace siren {
namespace detector {
namespace gdml {
// ---- Built-in constants (CLHEP/Geant4 unit system) ----
// Single source of truth for all constants seeded before user-defined <define>.
// Values are in CLHEP base units: mm=1 (length), MeV=1 (energy), radian=1
// (angle), nanosecond=1 (time), kelvin=1 (temp), mole=1 (amount).
// ParseLength/ParseAngle handle conversion from CLHEP units to SIREN units.
//
// Unit values derived from CLHEP SystemOfUnits.h and PhysicalConstants.h:
//   https://gitlab.cern.ch/CLHEP/CLHEP/-/blob/master/Units/Units/SystemOfUnits.h
// Math constants from CLHEP Evaluator (setStdMath):
//   https://gitlab.cern.ch/CLHEP/CLHEP/-/blob/master/Evaluator/src/Evaluator.cc
//
// GDML schema and element semantics per GDML User's Guide v2.9:
//   https://gdml.web.cern.ch/doc/GDMLmanual.pdf

struct BuiltInEntry {
    const char* name;
    double value;
};

static const BuiltInEntry kBuiltInConstants[] = {
    // Math (from CLHEP Evaluator::setStdMath)
    {"pi", M_PI},
    {"twopi", 2.0 * M_PI},
    {"TWOPI", 2.0 * M_PI},
    {"halfpi", M_PI / 2.0},
    {"e", 2.7182818284590452354},
    {"gamma", 0.5772156649015328606},

    // Angle (radian = 1)
    {"rad", 1.0},
    {"radian", 1.0},
    {"mrad", 1e-3},
    {"milliradian", 1e-3},
    {"deg", M_PI / 180.0},
    {"degree", M_PI / 180.0},

    // Length (mm = 1)
    {"mm", 1.0},
    {"millimeter", 1.0},
    {"cm", 10.0},
    {"centimeter", 10.0},
    {"m", 1000.0},
    {"meter", 1000.0},
    {"km", 1e6},
    {"kilometer", 1e6},
    {"micrometer", 1e-3},
    {"nanometer", 1e-6},
    {"angstrom", 1e-7},
    {"fermi", 1e-12},
    {"inch", 25.4},
    {"in", 25.4},

    // Energy (MeV = 1)
    {"eV", 1e-6},
    {"keV", 1e-3},
    {"MeV", 1.0},
    {"GeV", 1e3},
    {"TeV", 1e6},
    {"PeV", 1e9},
    {"joule", 6.241509074460763e12},

    // Time (nanosecond = 1)
    {"ns", 1.0},
    {"nanosecond", 1.0},
    {"ps", 1e-3},
    {"picosecond", 1e-3},
    {"us", 1e3},
    {"microsecond", 1e3},
    {"ms", 1e6},
    {"millisecond", 1e6},
    {"second", 1e9},

    // Mass (derived: kg = joule * s^2 / m^2)
    {"milligram", 6.241509074460763e18},
    {"gram", 6.241509074460763e21},
    {"kg", 6.241509074460763e24},
    {"kilogram", 6.241509074460763e24},

    // Volumic mass (CLHEP mass / mm^3)
    {"g/cm3", 6.241509074460763e18},
    {"mg/cm3", 6.241509074460763e15},
    {"kg/m3", 6.241509074460763e15},

    // Amount (mole = 1)
    {"mole", 1.0},
    {"mol", 1.0},

    // Temperature (kelvin = 1)
    {"kelvin", 1.0},
};

static const std::map<std::string, double>& GetBuiltInMap() {
    static const std::map<std::string, double> m = []() {
        std::map<std::string, double> m;
        for(auto const & e : kBuiltInConstants) m[e.name] = e.value;
        return m;
    }();
    return m;
}

static const std::unordered_set<std::string>& GetBuiltInNames() {
    static const std::unordered_set<std::string> names = []() {
        std::unordered_set<std::string> s;
        for(auto const & e : kBuiltInConstants) s.insert(e.name);
        return s;
    }();
    return names;
}

bool IsBuiltInConstant(std::string const & name) {
    return GetBuiltInNames().count(name) != 0;
}

void SeedBuiltInConstants(GDMLData & data) {
    for(auto const & e : kBuiltInConstants) {
        data.constants[e.name] = e.value;
    }
}

// Get the length scale factor for a given unit string
double LengthScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return siren::utilities::Constants::mm; // GDML default is mm
    }
    std::string u(unit);
    if(u == "mm") return siren::utilities::Constants::mm;
    if(u == "cm") return siren::utilities::Constants::cm;
    if(u == "m")  return siren::utilities::Constants::m;
    if(u == "um") return siren::utilities::Constants::um;
    if(u == "nm") return siren::utilities::Constants::nm;
    if(u == "km") return siren::utilities::Constants::km;
    throw std::runtime_error("GDML error: unrecognized length unit '" + u + "'");
}

// Convert a GDML length value+unit to meters (SIREN base unit)
// GDML default length unit is mm
double ParseLength(const char* value, const char* unit,
                   std::map<std::string, double> const & constants) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

    return val * LengthScale(unit);
}

// Get the angle scale factor for a given unit string
double AngleScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        // GDML default angle unit is radians
        return siren::utilities::Constants::radian;
    }
    std::string u(unit);
    if(u == "deg" || u == "degree") return siren::utilities::Constants::degrees;
    if(u == "rad" || u == "radian") return siren::utilities::Constants::radian;
    if(u == "mrad") return siren::utilities::Constants::mrad;
    throw std::runtime_error("GDML error: unrecognized angle unit '" + u + "'");
}

// Convert a GDML angle value+unit to radians
// GDML default angle unit is radians (matches AngleScale)
double ParseAngle(const char* value, const char* unit,
                  std::map<std::string, double> const & constants) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);
    return val * AngleScale(unit);
}

// Build a quaternion from GDML rotation convention:
// Extrinsic rotations: rotate around X by rx, then Y by ry, then Z by rz
// (each about the static world axes), so the composite quaternion is
// qz * qy * qx applied right-to-left.
//
// Built explicitly from the half-angle factors for clarity and to keep
// this convention self-contained; it equals siren::math::QFromXYZs and
// QuaternionFromEulerAngles(EulerAngles(XYZs, rx, ry, rz)).
Quaternion QuatFromGDMLRotation(double rx, double ry, double rz) {
    // Quaternion(x, y, z, w)
    Quaternion qx(std::sin(rx / 2.0), 0, 0, std::cos(rx / 2.0));
    Quaternion qy(0, std::sin(ry / 2.0), 0, std::cos(ry / 2.0));
    Quaternion qz(0, 0, std::sin(rz / 2.0), std::cos(rz / 2.0));
    return qz * qy * qx;
}

// Convert a value from the given unit to the GDML default for that unit type.
// For length: convert to mm (GDML default length)
// For angle: convert to rad (GDML default angle)
// For energy: convert to MeV (common default)
// Returns the scale factor to multiply the raw value by.
// Infer quantity type from a unit string.
static GDMLQuantityType TypeFromUnit(std::string const & u) {
    if(u.empty()) return GDMLQuantityType::NONE;
    if(u == "mm" || u == "cm" || u == "m" || u == "km" ||
       u == "um" || u == "nm" || u == "in" ||
       u == "millimeter" || u == "centimeter" || u == "meter" ||
       u == "kilometer" || u == "micrometer" || u == "nanometer" ||
       u == "angstrom" || u == "fermi" || u == "inch")
        return GDMLQuantityType::LENGTH;
    if(u == "rad" || u == "deg" || u == "mrad" ||
       u == "radian" || u == "degree" || u == "milliradian")
        return GDMLQuantityType::ANGLE;
    if(u == "eV" || u == "keV" || u == "MeV" || u == "GeV" ||
       u == "TeV" || u == "PeV" || u == "joule")
        return GDMLQuantityType::ENERGY;
    if(u == "ns" || u == "ps" || u == "us" || u == "ms" ||
       u == "nanosecond" || u == "picosecond" || u == "microsecond" ||
       u == "millisecond" || u == "second")
        return GDMLQuantityType::TIME;
    if(u == "g/cm3" || u == "mg/cm3" || u == "kg/m3")
        return GDMLQuantityType::DENSITY;
    if(u == "kelvin")
        return GDMLQuantityType::TEMPERATURE;
    if(u == "mole" || u == "mol")
        return GDMLQuantityType::AMOUNT;
    return GDMLQuantityType::NONE;
}

// Infer quantity type from the type attribute string.
static GDMLQuantityType TypeFromTypeAttr(std::string const & type_str) {
    if(type_str.empty()) return GDMLQuantityType::NONE;
    std::string t = type_str;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    if(t == "length" || t == "lenght" || t == "coordinate")
        return GDMLQuantityType::LENGTH;
    if(t == "angle") return GDMLQuantityType::ANGLE;
    if(t == "energy" || t == "threshold") return GDMLQuantityType::ENERGY;
    if(t == "time") return GDMLQuantityType::TIME;
    if(t == "density") return GDMLQuantityType::DENSITY;
    if(t == "temperature") return GDMLQuantityType::TEMPERATURE;
    return GDMLQuantityType::NONE;
}

// Resolve the quantity type from both the type attribute and the unit string.
// If both are present and disagree, warn and prefer the unit-inferred type.
GDMLQuantityType ResolveQuantityType(
        const char* type_attr, const char* unit,
        std::string const & name,
        GDMLData & data, GDMLParseOptions const & options) {
    GDMLQuantityType from_attr = TypeFromTypeAttr(type_attr ? type_attr : "");
    GDMLQuantityType from_unit = TypeFromUnit(unit ? unit : "");

    if(from_attr != GDMLQuantityType::NONE && from_unit != GDMLQuantityType::NONE) {
        if(from_attr != from_unit) {
            EmitWarning(data, options, "quantity '" + name + "': type '"
                + std::string(type_attr) + "' conflicts with unit '"
                + std::string(unit) + "'; using unit-inferred type");
        }
        return from_unit;
    }
    if(from_unit != GDMLQuantityType::NONE) return from_unit;
    return from_attr;
}

// Convert a quantity value from its declared unit to the canonical storage
// representation for its type.
//
// LENGTH/ANGLE/ENERGY/TIME: stored in CLHEP base units (mm, rad, MeV, ns).
//   These equal the GDML default units, so downstream consumption via
//   lscale/ascale works correctly.
// DENSITY: stored in g/cm3. The <D> parser expects g/cm3, so we must NOT
//   use CLHEP density units here.
double QuantityUnitScale(const char* unit, GDMLQuantityType resolved_type) {
    if(!unit || unit[0] == '\0') return 1.0;

    std::string u(unit);

    if(resolved_type == GDMLQuantityType::DENSITY) {
        if(u == "g/cm3") return 1.0;
        if(u == "mg/cm3") return 1.0e-3;
        if(u == "kg/m3") return 1.0e-3;
        throw std::runtime_error("GDML error: unrecognized density unit '" + u + "'");
    }

    auto it = GetBuiltInMap().find(u);
    if(it != GetBuiltInMap().end()) return it->second;
    throw std::runtime_error("GDML error: unrecognized unit '" + u + "'");
}

} // namespace gdml
} // namespace detector
} // namespace siren
