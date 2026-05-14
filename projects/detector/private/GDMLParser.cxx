#include "SIREN/detector/GDMLParser.h"

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>

#include <cereal/external/rapidxml/rapidxml.hpp>

namespace rapidxml = cereal::rapidxml;

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/Quaternion.h"

#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/ExtrPoly.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/BooleanGeometry.h"

using namespace siren::math;
using namespace siren::geometry;

namespace siren {
namespace detector {

namespace {

static const double PI = 3.141592653589793238462643383279502884197;

// Safe attribute value accessor: returns empty string if attribute is missing
const char* SafeAttrVal(rapidxml::xml_node<>* node, const char* attr_name) {
    if(!node) return "";
    rapidxml::xml_attribute<>* attr = node->first_attribute(attr_name);
    if(!attr) return "";
    return attr->value();
}

// Parse a double from a string, returning 0 if empty
double SafeParseDouble(const char* val) {
    if(!val || val[0] == '\0') return 0.0;
    return std::stod(std::string(val));
}

// Convert a GDML length value+unit to meters (SIREN base unit)
// GDML default length unit is mm
double ParseLength(const char* value, const char* unit) {
    if(!value || value[0] == '\0') return 0.0;
    double val = std::stod(std::string(value));

    if(!unit || unit[0] == '\0') {
        // GDML default is mm
        return val * 0.001;
    }

    std::string u(unit);
    if(u == "mm") return val * 0.001;
    if(u == "cm") return val * 0.01;
    if(u == "m")  return val * 1.0;

    // Unknown unit, assume mm
    return val * 0.001;
}

// Get the length scale factor for a given unit string
double LengthScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return 0.001; // GDML default is mm
    }
    std::string u(unit);
    if(u == "mm") return 0.001;
    if(u == "cm") return 0.01;
    if(u == "m")  return 1.0;
    return 0.001; // fallback to mm
}

// Convert a GDML angle value+unit to radians
// GDML default angle unit is radians
double ParseAngle(const char* value, const char* unit) {
    if(!value || value[0] == '\0') return 0.0;
    double val = std::stod(std::string(value));

    if(!unit || unit[0] == '\0') {
        // GDML default is radians
        return val;
    }

    std::string u(unit);
    if(u == "deg") return val * PI / 180.0;
    if(u == "rad") return val;

    // Unknown unit, assume radians
    return val;
}

// Get the angle scale factor for a given unit string
double AngleScale(const char* unit) {
    if(!unit || unit[0] == '\0') {
        return 1.0; // GDML default is radians
    }
    std::string u(unit);
    if(u == "deg") return PI / 180.0;
    if(u == "rad") return 1.0;
    return 1.0; // fallback to radians
}

// Build a quaternion from GDML rotation convention:
// Extrinsic rotations: rotate around X by rx, then Y by ry, then Z by rz
// Equivalent to intrinsic Z-Y-X, so quaternion = Qz * Qy * Qx
Quaternion QuatFromGDMLRotation(double rx, double ry, double rz) {
    // Individual axis quaternions: Quaternion(x, y, z, w)
    Quaternion qx(std::sin(rx / 2.0), 0, 0, std::cos(rx / 2.0));
    Quaternion qy(0, std::sin(ry / 2.0), 0, std::cos(ry / 2.0));
    Quaternion qz(0, 0, std::sin(rz / 2.0), std::cos(rz / 2.0));
    // Applied right-to-left: first X, then Y, then Z
    return qz * qy * qx;
}

} // anonymous namespace


// ---- Section parsers ----

// Parse the <define> section: constants, positions, rotations
static void ParseDefine(rapidxml::xml_node<>* define_node, GDMLData & data) {
    if(!define_node) return;

    for(auto* node = define_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());

        if(tag == "constant") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"));
            if(!name.empty()) {
                data.constants[name] = value;
            }
        }
        else if(tag == "quantity") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"));
            if(!name.empty()) {
                data.constants[name] = value;
            }
        }
        else if(tag == "position") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double x = ParseLength(SafeAttrVal(node, "x"), unit);
            double y = ParseLength(SafeAttrVal(node, "y"), unit);
            double z = ParseLength(SafeAttrVal(node, "z"), unit);
            if(!name.empty()) {
                data.positions[name] = Vector3D(x, y, z);
            }
        }
        else if(tag == "rotation") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double rx = ParseAngle(SafeAttrVal(node, "x"), unit);
            double ry = ParseAngle(SafeAttrVal(node, "y"), unit);
            double rz = ParseAngle(SafeAttrVal(node, "z"), unit);
            if(!name.empty()) {
                data.rotations[name] = QuatFromGDMLRotation(rx, ry, rz);
            }
        }
        else if(tag == "scale") {
            // Scale definitions are not directly used; skip
        }
        else if(tag == "variable") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"));
            if(!name.empty()) {
                data.constants[name] = value;
            }
        }
    }
}


// Parse the <materials> section: elements and materials
static void ParseMaterials(rapidxml::xml_node<>* materials_node, GDMLData & data) {
    if(!materials_node) return;

    for(auto* node = materials_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());

        if(tag == "isotope") {
            // Isotope definitions: store as simple materials
            GDMLMaterial mat;
            mat.name = SafeAttrVal(node, "name");
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"));
            mat.A = SafeParseDouble(SafeAttrVal(node, "N"));
            mat.density = 0.0;

            // Check for <atom> child
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"));
            }

            if(!mat.name.empty()) {
                data.materials[mat.name] = mat;
            }
        }
        else if(tag == "element") {
            GDMLMaterial mat;
            mat.name = SafeAttrVal(node, "name");
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"));
            mat.density = 0.0;

            // Check for <atom> child with mass
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                // atom value might have units like "g/mol"
                const char* a_val = SafeAttrVal(atom_node, "value");
                mat.A = SafeParseDouble(a_val);
            }

            // Check for <fraction> children (element defined as isotope mix)
            for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                std::string ref = SafeAttrVal(frac, "ref");
                double n = SafeParseDouble(SafeAttrVal(frac, "n"));
                if(!ref.empty()) {
                    mat.composition[ref] = n;
                }
            }

            if(!mat.name.empty()) {
                data.materials[mat.name] = mat;
            }
        }
        else if(tag == "material") {
            GDMLMaterial mat;
            mat.name = SafeAttrVal(node, "name");
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"));
            mat.density = 0.0;
            mat.A = 0.0;

            // Parse <D> (density) child
            auto* d_node = node->first_node("D");
            if(d_node) {
                mat.density = SafeParseDouble(SafeAttrVal(d_node, "value"));
                // Density unit: typically g/cm3, which is the SIREN convention
            }

            // Parse <atom> child for simple materials
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"));
            }

            // Parse <fraction> children for composite materials
            for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                std::string ref = SafeAttrVal(frac, "ref");
                double n = SafeParseDouble(SafeAttrVal(frac, "n"));
                if(!ref.empty()) {
                    mat.composition[ref] = n;
                }
            }

            // Parse <composite> children (alternative to fraction)
            for(auto* comp = node->first_node("composite"); comp; comp = comp->next_sibling("composite")) {
                std::string ref = SafeAttrVal(comp, "ref");
                double n = SafeParseDouble(SafeAttrVal(comp, "n"));
                if(!ref.empty()) {
                    mat.composition[ref] = n;
                }
            }

            if(!mat.name.empty()) {
                data.materials[mat.name] = mat;
            }
        }
    }
}


// Parse the <solids> section: convert GDML solids to SIREN Geometry objects
static void ParseSolids(rapidxml::xml_node<>* solids_node, GDMLData & data) {
    if(!solids_node) return;

    for(auto* node = solids_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());
        std::string name = SafeAttrVal(node, "name");
        if(name.empty()) continue;

        std::shared_ptr<Geometry> geo;

        const char* lunit = SafeAttrVal(node, "lunit");
        const char* aunit = SafeAttrVal(node, "aunit");
        double lscale = LengthScale(lunit);
        double ascale = AngleScale(aunit);

        if(tag == "box") {
            // GDML <box> x,y,z are HALF-widths; SIREN Box takes FULL widths
            double hx = SafeParseDouble(SafeAttrVal(node, "x")) * lscale;
            double hy = SafeParseDouble(SafeAttrVal(node, "y")) * lscale;
            double hz = SafeParseDouble(SafeAttrVal(node, "z")) * lscale;
            geo = Box(hx * 2.0, hy * 2.0, hz * 2.0).create();
        }
        else if(tag == "sphere") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin")) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax")) * lscale;
            geo = Sphere(rmax, rmin).create();
        }
        else if(tag == "orb") {
            double r = SafeParseDouble(SafeAttrVal(node, "r")) * lscale;
            geo = Sphere(r, 0).create();
        }
        else if(tag == "tube" || tag == "tubs") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin")) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax")) * lscale;
            // GDML <tube> z is HALF-height; SIREN Cylinder takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z")) * lscale;
            geo = Cylinder(rmax, rmin, hz * 2.0).create();
        }
        else if(tag == "cone") {
            double rmin1 = SafeParseDouble(SafeAttrVal(node, "rmin1")) * lscale;
            double rmax1 = SafeParseDouble(SafeAttrVal(node, "rmax1")) * lscale;
            double rmin2 = SafeParseDouble(SafeAttrVal(node, "rmin2")) * lscale;
            double rmax2 = SafeParseDouble(SafeAttrVal(node, "rmax2")) * lscale;
            // GDML <cone> z is HALF-height; SIREN Cone takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z")) * lscale;
            geo = Cone(rmin1, rmax1, rmin2, rmax2, hz * 2.0).create();
        }
        else if(tag == "trd") {
            // GDML <trd> uses half-widths; SIREN Trd also uses half-widths
            double dx1 = SafeParseDouble(SafeAttrVal(node, "x1")) * lscale;
            double dx2 = SafeParseDouble(SafeAttrVal(node, "x2")) * lscale;
            double dy1 = SafeParseDouble(SafeAttrVal(node, "y1")) * lscale;
            double dy2 = SafeParseDouble(SafeAttrVal(node, "y2")) * lscale;
            double dz  = SafeParseDouble(SafeAttrVal(node, "z"))  * lscale;
            geo = Trd(dx1, dx2, dy1, dy2, dz).create();
        }
        else if(tag == "polycone") {
            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi")) * ascale;
            (void)startphi; // Polycone in SIREN is axially symmetric

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z")) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin")) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax")) * lscale;
                z_planes.push_back(z);
                rmin_vec.push_back(rmin);
                rmax_vec.push_back(rmax);
            }

            if(z_planes.size() >= 2) {
                geo = Polycone(z_planes, rmin_vec, rmax_vec).create();
            }
        }
        else if(tag == "polyhedra") {
            int numSide = std::stoi(std::string(SafeAttrVal(node, "numsides")));
            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi")) * ascale;

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z")) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin")) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax")) * lscale;
                z_planes.push_back(z);
                rmin_vec.push_back(rmin);
                rmax_vec.push_back(rmax);
            }

            if(z_planes.size() >= 2 && numSide >= 3) {
                geo = Polyhedra(numSide, startphi, z_planes, rmin_vec, rmax_vec).create();
            }
        }
        else if(tag == "xtru") {
            std::vector<std::vector<double>> polygon;
            std::vector<ExtrPoly::ZSection> zsections;

            for(auto* vtx = node->first_node("twoDimVertex"); vtx; vtx = vtx->next_sibling("twoDimVertex")) {
                double vx = SafeParseDouble(SafeAttrVal(vtx, "x")) * lscale;
                double vy = SafeParseDouble(SafeAttrVal(vtx, "y")) * lscale;
                std::vector<double> vert;
                vert.push_back(vx);
                vert.push_back(vy);
                polygon.push_back(vert);
            }

            for(auto* sec = node->first_node("section"); sec; sec = sec->next_sibling("section")) {
                double zpos = SafeParseDouble(SafeAttrVal(sec, "zPosition")) * lscale;
                double xoff = SafeParseDouble(SafeAttrVal(sec, "xOffset")) * lscale;
                double yoff = SafeParseDouble(SafeAttrVal(sec, "yOffset")) * lscale;
                double scale = SafeParseDouble(SafeAttrVal(sec, "scalingFactor"));
                if(scale == 0.0) scale = 1.0;
                double offset[2] = {xoff, yoff};
                zsections.push_back(ExtrPoly::ZSection(zpos, offset, scale));
            }

            if(!polygon.empty() && !zsections.empty()) {
                geo = ExtrPoly(polygon, zsections).create();
            }
        }
        else if(tag == "subtraction" || tag == "union" || tag == "intersection") {
            // Boolean solids: reference first and second solids
            BooleanOperation op;
            if(tag == "subtraction")   op = BooleanOperation::SUBTRACTION;
            else if(tag == "union")    op = BooleanOperation::UNION;
            else                       op = BooleanOperation::INTERSECTION;

            auto* first_node = node->first_node("first");
            auto* second_node = node->first_node("second");

            if(first_node && second_node) {
                std::string first_ref = SafeAttrVal(first_node, "ref");
                std::string second_ref = SafeAttrVal(second_node, "ref");

                auto first_it = data.solids.find(first_ref);
                auto second_it = data.solids.find(second_ref);

                if(first_it != data.solids.end() && second_it != data.solids.end()) {
                    auto left = first_it->second;

                    // The second solid may have a position and rotation
                    // relative to the first
                    Vector3D rel_pos(0, 0, 0);
                    Quaternion rel_rot;

                    // Check for inline <position> child
                    auto* pos_node = node->first_node("position");
                    if(pos_node) {
                        const char* punit = SafeAttrVal(pos_node, "unit");
                        double px = ParseLength(SafeAttrVal(pos_node, "x"), punit);
                        double py = ParseLength(SafeAttrVal(pos_node, "y"), punit);
                        double pz = ParseLength(SafeAttrVal(pos_node, "z"), punit);
                        rel_pos = Vector3D(px, py, pz);
                    }

                    // Check for <positionref>
                    auto* posref_node = node->first_node("positionref");
                    if(posref_node) {
                        std::string ref = SafeAttrVal(posref_node, "ref");
                        auto it = data.positions.find(ref);
                        if(it != data.positions.end()) {
                            rel_pos = it->second;
                        }
                    }

                    // Check for inline <rotation> child
                    auto* rot_node = node->first_node("rotation");
                    if(rot_node) {
                        const char* runit = SafeAttrVal(rot_node, "unit");
                        double rx = ParseAngle(SafeAttrVal(rot_node, "x"), runit);
                        double ry = ParseAngle(SafeAttrVal(rot_node, "y"), runit);
                        double rz = ParseAngle(SafeAttrVal(rot_node, "z"), runit);
                        rel_rot = QuatFromGDMLRotation(rx, ry, rz);
                    }

                    // Check for <rotationref>
                    auto* rotref_node = node->first_node("rotationref");
                    if(rotref_node) {
                        std::string ref = SafeAttrVal(rotref_node, "ref");
                        auto it = data.rotations.find(ref);
                        if(it != data.rotations.end()) {
                            rel_rot = it->second;
                        }
                    }

                    // Create the second solid with relative placement
                    auto right = second_it->second->create();
                    Placement rel_placement(rel_pos, rel_rot);
                    right->SetPlacement(rel_placement);

                    geo = std::make_shared<BooleanGeometry>(op,
                        std::const_pointer_cast<const Geometry>(left),
                        std::const_pointer_cast<const Geometry>(right));
                }
            }
        }
        // Unrecognized solid types are silently skipped

        if(geo) {
            data.solids[name] = geo;
        }
    }
}


// Parse the <structure> section: volumes and physical volumes
static void ParseStructure(rapidxml::xml_node<>* structure_node, GDMLData & data) {
    if(!structure_node) return;

    for(auto* vol_node = structure_node->first_node("volume"); vol_node;
        vol_node = vol_node->next_sibling("volume")) {

        GDMLVolume vol;
        vol.name = SafeAttrVal(vol_node, "name");

        // Get <materialref>
        auto* matref = vol_node->first_node("materialref");
        if(matref) {
            vol.material_ref = SafeAttrVal(matref, "ref");
        }

        // Get <solidref>
        auto* solidref = vol_node->first_node("solidref");
        if(solidref) {
            vol.solid_ref = SafeAttrVal(solidref, "ref");
        }

        // Parse <physvol> children
        for(auto* pv = vol_node->first_node("physvol"); pv; pv = pv->next_sibling("physvol")) {
            GDMLPhysVol physvol;
            physvol.position = Vector3D(0, 0, 0);
            physvol.rotation = Quaternion(); // identity

            // Volume reference
            auto* volref = pv->first_node("volumeref");
            if(volref) {
                physvol.volume_ref = SafeAttrVal(volref, "ref");
            }

            // Position: inline or reference
            auto* pos = pv->first_node("position");
            if(pos) {
                const char* punit = SafeAttrVal(pos, "unit");
                double px = ParseLength(SafeAttrVal(pos, "x"), punit);
                double py = ParseLength(SafeAttrVal(pos, "y"), punit);
                double pz = ParseLength(SafeAttrVal(pos, "z"), punit);
                physvol.position = Vector3D(px, py, pz);
            }
            auto* posref = pv->first_node("positionref");
            if(posref) {
                std::string ref = SafeAttrVal(posref, "ref");
                auto it = data.positions.find(ref);
                if(it != data.positions.end()) {
                    physvol.position = it->second;
                }
            }

            // Rotation: inline or reference
            auto* rot = pv->first_node("rotation");
            if(rot) {
                const char* runit = SafeAttrVal(rot, "unit");
                double rx = ParseAngle(SafeAttrVal(rot, "x"), runit);
                double ry = ParseAngle(SafeAttrVal(rot, "y"), runit);
                double rz = ParseAngle(SafeAttrVal(rot, "z"), runit);
                physvol.rotation = QuatFromGDMLRotation(rx, ry, rz);
            }
            auto* rotref = pv->first_node("rotationref");
            if(rotref) {
                std::string ref = SafeAttrVal(rotref, "ref");
                auto it = data.rotations.find(ref);
                if(it != data.rotations.end()) {
                    physvol.rotation = it->second;
                }
            }

            vol.children.push_back(physvol);
        }

        if(!vol.name.empty()) {
            data.volumes[vol.name] = vol;
        }
    }
}


// Parse the <setup> section to find the world volume
static void ParseSetup(rapidxml::xml_node<>* setup_node, GDMLData & data) {
    if(!setup_node) return;

    auto* world_node = setup_node->first_node("world");
    if(world_node) {
        data.world_volume = SafeAttrVal(world_node, "ref");
    }
}


// Main parser function
GDMLData ParseGDML(std::string const & filename) {
    // Read the file into a mutable buffer (rapidxml modifies in-place)
    std::ifstream file(filename.c_str(), std::ios::binary);
    if(!file.is_open()) {
        throw std::runtime_error("Cannot open GDML file: " + filename);
    }

    // Read entire file into string
    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    file.close();

    // rapidxml requires a mutable, null-terminated buffer
    std::vector<char> buffer(content.begin(), content.end());
    buffer.push_back('\0');

    // Parse the XML
    rapidxml::xml_document<> doc;
    try {
        doc.parse<0>(buffer.data());
    } catch(rapidxml::parse_error & e) {
        throw std::runtime_error(std::string("GDML XML parse error: ") + e.what());
    }

    // Find the root <gdml> node
    auto* gdml_node = doc.first_node("gdml");
    if(!gdml_node) {
        throw std::runtime_error("No <gdml> root element found in: " + filename);
    }

    GDMLData data;

    // Parse sections in order
    ParseDefine(gdml_node->first_node("define"), data);
    ParseMaterials(gdml_node->first_node("materials"), data);
    ParseSolids(gdml_node->first_node("solids"), data);
    ParseStructure(gdml_node->first_node("structure"), data);
    ParseSetup(gdml_node->first_node("setup"), data);

    return data;
}

} // namespace detector
} // namespace siren
