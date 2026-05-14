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
    std::vector<GDMLPhysVol> children;
};

struct GDMLData {
    // Defined positions and rotations (resolved to SIREN units)
    std::map<std::string, math::Vector3D> positions;
    std::map<std::string, math::Quaternion> rotations;
    std::map<std::string, double> constants;

    // Materials
    std::map<std::string, GDMLMaterial> materials;

    // Solids (as Geometry objects, in SIREN units)
    std::map<std::string, std::shared_ptr<geometry::Geometry>> solids;

    // Volume hierarchy
    std::map<std::string, GDMLVolume> volumes;

    // World volume name (from <setup>)
    std::string world_volume;
};

// Parse a GDML file and return the parsed data
GDMLData ParseGDML(std::string const & filename);

} // namespace detector
} // namespace siren

#endif // SIREN_GDMLParser_H
