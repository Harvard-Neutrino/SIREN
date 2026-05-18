#pragma once
#ifndef SIREN_GDMLParser_H
#define SIREN_GDMLParser_H

#include <map>
#include <string>
#include <vector>
#include <memory>

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Quaternion.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace detector {

// Parsed GDML data, ready for conversion to DetectorModel sectors
struct GDMLMaterial {
    std::string name;
    double density; // g/cm^3
    double Z;       // atomic number (for simple materials)
    double A;       // atomic mass (g/mol, for simple materials)
    // For composite materials: map of component name -> mass fraction
    std::map<std::string, double> composition;

    bool operator==(GDMLMaterial const & o) const {
        return name == o.name && density == o.density
            && Z == o.Z && A == o.A && composition == o.composition;
    }
};

struct GDMLPhysVol {
    std::string volume_ref;
    math::Vector3D position;     // in meters (converted from GDML units)
    math::Quaternion rotation;   // converted from GDML Euler angles
};

struct GDMLVolume {
    std::string name;
    std::string solid_ref;
    std::string material_ref;
    bool is_assembly = false;
    std::vector<GDMLPhysVol> children;
};

// Physical dimension type for a GDML quantity.
// Used to check that typed quantities are consumed in compatible contexts
// and to select the correct unit-conversion path (e.g. density -> g/cm3).
enum class GDMLQuantityType {
    NONE, LENGTH, ANGLE, ENERGY, TIME, DENSITY, TEMPERATURE, AMOUNT
};

// Options controlling GDML parse behavior
struct GDMLParseOptions {
    // If true, any warning condition throws std::runtime_error instead of
    // collecting a warning and continuing. Useful for validating GDML files.
    bool strict = false;
};

struct GDMLData {
    // Defined positions and rotations (resolved to SIREN units)
    std::map<std::string, math::Vector3D> positions;
    std::map<std::string, math::Quaternion> rotations;
    std::map<std::string, double> constants;

    // Type information for named quantities. Keys are a subset of constants.
    // Only quantities with explicit type/unit have entries here.
    std::map<std::string, GDMLQuantityType> quantity_types;

    // Materials
    std::map<std::string, GDMLMaterial> materials;

    // Matrix definitions: name -> {coldim, flat_values}
    std::map<std::string, std::pair<int, std::vector<double>>> matrices;

    // Solids (as Geometry objects, in SIREN units)
    std::map<std::string, std::shared_ptr<geometry::Geometry>> solids;

    // Volume hierarchy
    std::map<std::string, GDMLVolume> volumes;

    // World volume name (from <setup>)
    std::string world_volume;

    // Instance counts: how many definitions of each name were seen during parsing.
    // Names with count > 1 are ambiguous for direct volume references.
    // Boolean solids resolve operand references using instance tracking at
    // definition time, so they are not affected by ambiguity.
    std::map<std::string, int> solid_instance_counts;
    std::map<std::string, int> material_instance_counts;

    // Warnings collected during parsing (unsupported features, skipped solids, etc.)
    // In strict mode this stays empty since warnings throw instead.
    std::vector<std::string> warnings;
};

// Parse a GDML file and return the parsed data.
// Accepts any root element (e.g. <gdml>, <gdml_simple_extension>) as long as
// it contains the standard GDML child sections (define, materials, solids, etc.).
// Multiple <solids>/<structure>/<define>/<materials> sections are merged.
GDMLData ParseGDML(std::string const & filename, GDMLParseOptions const & options = {});

} // namespace detector
} // namespace siren

#endif // SIREN_GDMLParser_H
