#pragma once
#ifndef SIREN_GDMLParser_H
#define SIREN_GDMLParser_H

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>

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

    // If true (default), strip Geant4 pointer-address suffixes from names
    // (e.g. "Si0x7f3a4c002a" -> "Si") and dedup identical definitions.
    // Disable if the stripping interferes with intentional naming.
    bool strip_pointer_suffixes = true;
};

struct GDMLData {
    // Defined positions and rotations (resolved to SIREN units)
    std::unordered_map<std::string, math::Vector3D> positions;
    std::unordered_map<std::string, math::Quaternion> rotations;
    std::unordered_map<std::string, double> constants;

    // Type information for named quantities. Keys are a subset of constants.
    // Only quantities with explicit type/unit have entries here.
    std::unordered_map<std::string, GDMLQuantityType> quantity_types;

    // GDML material types, kept in separate namespaces to avoid collisions
    // between e.g. <element name="Si"> and <material name="Si">.
    std::unordered_map<std::string, GDMLMaterial> isotopes;
    std::unordered_map<std::string, GDMLMaterial> elements;
    std::unordered_map<std::string, GDMLMaterial> materials;

    // Matrix definitions: name -> {coldim, flat_values}
    std::unordered_map<std::string, std::pair<int, std::vector<double>>> matrices;

    // Solids (as Geometry objects, in SIREN units)
    std::unordered_map<std::string, std::shared_ptr<geometry::Geometry>> solids;

    // Volume hierarchy
    std::unordered_map<std::string, GDMLVolume> volumes;

    // World volume name (from <setup>)
    std::string world_volume;

    // Base directory of the GDML file, for resolving relative <file> paths.
    std::string base_dir;

    // Name aliases: maps original names (possibly with Geant4 pointer suffixes
    // like "Si0x7f3a") to their canonical stripped/deduped names ("Si").
    // Used to resolve references that use the original exported names.
    std::unordered_map<std::string, std::string> name_aliases;

    // Instance counts: how many definitions of each name were seen during parsing.
    // Names with count > 1 are ambiguous for direct volume references.
    // Boolean solids resolve operand references using instance tracking at
    // definition time, so they are not affected by ambiguity.
    std::unordered_map<std::string, int> solid_instance_counts;
    std::unordered_map<std::string, int> material_instance_counts;

    // Warnings collected during parsing (unsupported features, skipped solids, etc.)
    // In strict mode this stays empty since warnings throw instead.
    std::vector<std::string> warnings;

    // Fatal errors collected during parsing. Parsing continues to gather as
    // many diagnostics as possible, then throws with the full list at the end
    // of each parse phase. Check with ThrowIfErrors() after each phase.
    std::vector<std::string> errors;

    // Throw a std::runtime_error if any errors have been collected, with a
    // message listing all of them. Clears the error list after throwing.
    void ThrowIfErrors();

    // Merge another GDMLData (typically from a <file> sub-parse) into this
    // one, handling name collisions automatically.
    //
    // Identical definitions (same material properties, same solid geometry)
    // are silently reused. Differing definitions with the same name are
    // auto-renamed with a numeric suffix in the sub-data before merging.
    //
    // Returns a rename map (original name -> new name) for any names that
    // were changed, so the caller can update volume references.
    std::unordered_map<std::string, std::string> MergeFrom(GDMLData & other);
};

/// Parse a GDML file and return the parsed data.
/// Accepts any root element (e.g. <gdml>, <gdml_simple_extension>) as long as
/// it contains the standard GDML child sections (define, materials, solids, etc.).
/// Multiple <solids>/<structure>/<define>/<materials> sections are merged.
///
/// @warning ParseGDML resolves ENTITY SYSTEM and xi:include directives
/// by reading files relative to the GDML file's directory, including
/// absolute paths. Only call this on trusted GDML input.
GDMLData ParseGDML(std::string const & filename, GDMLParseOptions const & options = {});

} // namespace detector
} // namespace siren

#endif // SIREN_GDMLParser_H
