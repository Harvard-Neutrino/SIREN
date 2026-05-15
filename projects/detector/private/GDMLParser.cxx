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
#include <iostream>
#include <map>
#include <cctype>

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

// Trim leading and trailing whitespace from a string
std::string TrimWhitespace(std::string const & s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if(start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

// Evaluate a GDML expression string that may contain:
//   - numeric literals: "3.14"
//   - constant references: "det_length"
//   - simple binary ops: "3.14/2", "det_length+10"
//   - negation: "-det_length"
//   - parenthesized sub-expressions: "(a+b)*2"
// Operators handled: +, -, *, /
// Precedence: * and / bind tighter than + and -
double EvalExpression(std::string const & expr, std::map<std::string, double> const & constants) {
    std::string s = TrimWhitespace(expr);
    if(s.empty()) return 0.0;

    // Try simple numeric literal first (fast path)
    {
        size_t pos = 0;
        try {
            double val = std::stod(s, &pos);
            if(pos == s.size()) return val;
        } catch(std::invalid_argument &) {
            // Not a simple number, continue to expression parsing
        } catch(std::out_of_range &) {
            return 0.0;
        }
    }

    // Try constant lookup for bare names (no operators)
    {
        bool is_name = true;
        for(size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if(!std::isalnum(static_cast<unsigned char>(c)) && c != '_') {
                is_name = false;
                break;
            }
        }
        if(is_name) {
            auto it = constants.find(s);
            if(it != constants.end()) return it->second;
            return 0.0; // Unknown constant
        }
    }

    // Handle unary negation at the start: "-expr"
    if(s[0] == '-' && s.size() > 1) {
        // Check if this is just a negative number (already handled above via stod)
        // so this must be negation of an expression like "-det_length"
        return -EvalExpression(s.substr(1), constants);
    }
    if(s[0] == '+' && s.size() > 1) {
        return EvalExpression(s.substr(1), constants);
    }

    // Strip outer parentheses if they match
    if(s.front() == '(' && s.back() == ')') {
        int depth = 0;
        bool matched = true;
        for(size_t i = 0; i < s.size(); ++i) {
            if(s[i] == '(') ++depth;
            else if(s[i] == ')') --depth;
            if(depth == 0 && i < s.size() - 1) {
                matched = false;
                break;
            }
        }
        if(matched) {
            return EvalExpression(s.substr(1, s.size() - 2), constants);
        }
    }

    // Find the lowest-precedence operator not inside parentheses.
    // Scan right-to-left for + or - (lowest precedence, left-associative),
    // then for * or / if no + or - found.
    auto findSplitOp = [&](std::string const & ops) -> int {
        int depth = 0;
        for(int i = static_cast<int>(s.size()) - 1; i >= 0; --i) {
            if(s[i] == ')') ++depth;
            else if(s[i] == '(') --depth;
            if(depth != 0) continue;

            for(char op : ops) {
                if(s[i] == op) {
                    // Don't split on a unary +/- at position 0
                    if(i == 0) continue;
                    // Don't split on +/- after another operator (e.g. "3*-2")
                    char prev = s[i - 1];
                    if(prev == '+' || prev == '-' || prev == '*' || prev == '/') continue;
                    return i;
                }
            }
        }
        return -1;
    };

    // Try + and - first (lowest precedence)
    int splitPos = findSplitOp("+-");
    if(splitPos > 0) {
        double left = EvalExpression(s.substr(0, splitPos), constants);
        double right = EvalExpression(s.substr(splitPos + 1), constants);
        if(s[splitPos] == '+') return left + right;
        else return left - right;
    }

    // Try * and /
    splitPos = findSplitOp("*/");
    if(splitPos > 0) {
        double left = EvalExpression(s.substr(0, splitPos), constants);
        double right = EvalExpression(s.substr(splitPos + 1), constants);
        if(s[splitPos] == '*') return left * right;
        else return (right != 0.0) ? left / right : 0.0;
    }

    // Could not parse -- return 0
    return 0.0;
}

// Parse a double from a string, returning 0 if empty.
// Optionally evaluates expressions using the provided constants map.
double SafeParseDouble(const char* val, std::map<std::string, double> const & constants = {}) {
    if(!val || val[0] == '\0') return 0.0;
    // Fast path: try simple numeric literal
    try {
        size_t pos = 0;
        double result = std::stod(std::string(val), &pos);
        std::string s(val);
        if(pos == s.size()) return result;
    } catch(std::invalid_argument &) {
        // Not a simple number
    } catch(std::out_of_range &) {
        return 0.0;
    }
    // Fall back to expression evaluation
    return EvalExpression(std::string(val), constants);
}

// Convert a GDML length value+unit to meters (SIREN base unit)
// GDML default length unit is mm
double ParseLength(const char* value, const char* unit,
                   std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

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
double ParseAngle(const char* value, const char* unit,
                  std::map<std::string, double> const & constants = {}) {
    if(!value || value[0] == '\0') return 0.0;
    double val = EvalExpression(std::string(value), constants);

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
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
            if(!name.empty()) {
                data.constants[name] = value;
            }
        }
        else if(tag == "quantity") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
            if(!name.empty()) {
                data.constants[name] = value;
            }
        }
        else if(tag == "position") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double x = ParseLength(SafeAttrVal(node, "x"), unit, data.constants);
            double y = ParseLength(SafeAttrVal(node, "y"), unit, data.constants);
            double z = ParseLength(SafeAttrVal(node, "z"), unit, data.constants);
            if(!name.empty()) {
                data.positions[name] = Vector3D(x, y, z);
            }
        }
        else if(tag == "rotation") {
            std::string name = SafeAttrVal(node, "name");
            const char* unit = SafeAttrVal(node, "unit");
            double rx = ParseAngle(SafeAttrVal(node, "x"), unit, data.constants);
            double ry = ParseAngle(SafeAttrVal(node, "y"), unit, data.constants);
            double rz = ParseAngle(SafeAttrVal(node, "z"), unit, data.constants);
            if(!name.empty()) {
                data.rotations[name] = QuatFromGDMLRotation(rx, ry, rz);
            }
        }
        else if(tag == "scale") {
            // Scale definitions are not directly used; skip
        }
        else if(tag == "variable") {
            std::string name = SafeAttrVal(node, "name");
            double value = SafeParseDouble(SafeAttrVal(node, "value"), data.constants);
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
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
            mat.A = SafeParseDouble(SafeAttrVal(node, "N"), data.constants);
            mat.density = 0.0;

            // Check for <atom> child
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"), data.constants);
            }

            if(!mat.name.empty()) {
                data.materials[mat.name] = mat;
            }
        }
        else if(tag == "element") {
            GDMLMaterial mat;
            mat.name = SafeAttrVal(node, "name");
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
            mat.density = 0.0;

            // Check for <atom> child with mass
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                // atom value might have units like "g/mol"
                const char* a_val = SafeAttrVal(atom_node, "value");
                mat.A = SafeParseDouble(a_val, data.constants);
            }

            // Check for <fraction> children (element defined as isotope mix)
            for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                std::string ref = SafeAttrVal(frac, "ref");
                double n = SafeParseDouble(SafeAttrVal(frac, "n"), data.constants);
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
            mat.Z = SafeParseDouble(SafeAttrVal(node, "Z"), data.constants);
            mat.density = 0.0;
            mat.A = 0.0;

            // Parse <D> (density) child
            auto* d_node = node->first_node("D");
            if(d_node) {
                mat.density = SafeParseDouble(SafeAttrVal(d_node, "value"), data.constants);
                // Density unit scaling (SIREN convention is g/cm3)
                const char* d_unit = SafeAttrVal(d_node, "unit");
                if(d_unit && d_unit[0] != '\0') {
                    std::string du(d_unit);
                    if(du == "kg/m3") {
                        mat.density /= 1000.0;
                    } else if(du == "mg/cm3") {
                        mat.density /= 1000.0;
                    } else if(du == "g/cm3") {
                        // No conversion needed
                    } else {
                        std::cerr << "GDML warning: unrecognized density unit '"
                                  << du << "' for material '" << mat.name
                                  << "', assuming g/cm3" << std::endl;
                    }
                }
                // If unit is missing, assume g/cm3 (no conversion)
            }

            // Parse <atom> child for simple materials
            auto* atom_node = node->first_node("atom");
            if(atom_node) {
                mat.A = SafeParseDouble(SafeAttrVal(atom_node, "value"), data.constants);
            }

            // Parse <fraction> children for composite materials
            for(auto* frac = node->first_node("fraction"); frac; frac = frac->next_sibling("fraction")) {
                std::string ref = SafeAttrVal(frac, "ref");
                double n = SafeParseDouble(SafeAttrVal(frac, "n"), data.constants);
                if(!ref.empty()) {
                    mat.composition[ref] = n;
                }
            }

            // Parse <composite> children (alternative to fraction)
            for(auto* comp = node->first_node("composite"); comp; comp = comp->next_sibling("composite")) {
                std::string ref = SafeAttrVal(comp, "ref");
                double n = SafeParseDouble(SafeAttrVal(comp, "n"), data.constants);
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


// Parse the <solids> section: convert GDML solids to SIREN Geometry objects.
// Uses two passes: first non-boolean solids, then boolean solids so that
// forward references to operands are resolved correctly.
static void ParseSolids(rapidxml::xml_node<>* solids_node, GDMLData & data) {
    if(!solids_node) return;

    // Tolerance for comparing angular values to full-rotation defaults
    static const double ANG_TOL = 1e-6;

    // Helper lambda: parse a single non-boolean solid node and return geometry
    auto parsePrimitive = [&](rapidxml::xml_node<>* node) -> std::shared_ptr<Geometry> {
        std::string tag(node->name());
        std::shared_ptr<Geometry> geo;

        const char* lunit = SafeAttrVal(node, "lunit");
        const char* aunit = SafeAttrVal(node, "aunit");
        double lscale = LengthScale(lunit);
        double ascale = AngleScale(aunit);
        std::string name = SafeAttrVal(node, "name");

        if(tag == "box") {
            // GDML <box> x,y,z are HALF-widths; SIREN Box takes FULL widths
            double hx = SafeParseDouble(SafeAttrVal(node, "x"), data.constants) * lscale;
            double hy = SafeParseDouble(SafeAttrVal(node, "y"), data.constants) * lscale;
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            geo = Box(hx * 2.0, hy * 2.0, hz * 2.0).create();
        }
        else if(tag == "sphere") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;

            // Check angular parameters for partial sphere
            double starttheta = SafeParseDouble(SafeAttrVal(node, "starttheta"), data.constants) * ascale;
            double deltatheta = SafeParseDouble(SafeAttrVal(node, "deltatheta"), data.constants) * ascale;
            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi"), data.constants) * ascale;
            double deltaphi = SafeParseDouble(SafeAttrVal(node, "deltaphi"), data.constants) * ascale;
            // Default GDML sphere: starttheta=0, deltatheta=pi, startphi=0, deltaphi=2*pi
            // Use defaults if the attribute was not specified (value is 0 from SafeParseDouble)
            const char* st_val = SafeAttrVal(node, "starttheta");
            const char* dt_val = SafeAttrVal(node, "deltatheta");
            const char* sp_val = SafeAttrVal(node, "startphi");
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if((st_val[0] != '\0' && std::fabs(starttheta) > ANG_TOL) ||
               (dt_val[0] != '\0' && std::fabs(deltatheta - PI) > ANG_TOL) ||
               (sp_val[0] != '\0' && std::fabs(startphi) > ANG_TOL) ||
               (dp_val[0] != '\0' && std::fabs(deltaphi - 2.0 * PI) > ANG_TOL)) {
                std::cerr << "GDML warning: sphere '" << name
                          << "' has non-full-rotation angular parameters"
                          << " (starttheta=" << starttheta
                          << " deltatheta=" << deltatheta
                          << " startphi=" << startphi
                          << " deltaphi=" << deltaphi
                          << "); SIREN does not support partial spheres,"
                          << " creating full sphere" << std::endl;
            }

            geo = Sphere(rmax, rmin).create();
        }
        else if(tag == "orb") {
            double r = SafeParseDouble(SafeAttrVal(node, "r"), data.constants) * lscale;
            geo = Sphere(r, 0).create();
        }
        else if(tag == "tube" || tag == "tubs") {
            double rmin = SafeParseDouble(SafeAttrVal(node, "rmin"), data.constants) * lscale;
            double rmax = SafeParseDouble(SafeAttrVal(node, "rmax"), data.constants) * lscale;
            // GDML <tube> z is HALF-height; SIREN Cylinder takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            // Check deltaphi for partial tube
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    std::cerr << "GDML warning: " << tag << " '" << name
                              << "' has deltaphi=" << deltaphi
                              << " (not 2*pi); SIREN does not support partial tubes,"
                              << " creating full cylinder" << std::endl;
                }
            }

            geo = Cylinder(rmax, rmin, hz * 2.0).create();
        }
        else if(tag == "cone") {
            double rmin1 = SafeParseDouble(SafeAttrVal(node, "rmin1"), data.constants) * lscale;
            double rmax1 = SafeParseDouble(SafeAttrVal(node, "rmax1"), data.constants) * lscale;
            double rmin2 = SafeParseDouble(SafeAttrVal(node, "rmin2"), data.constants) * lscale;
            double rmax2 = SafeParseDouble(SafeAttrVal(node, "rmax2"), data.constants) * lscale;
            // GDML <cone> z is HALF-height; SIREN Cone takes full height
            double hz = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;

            // Check deltaphi for partial cone
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    std::cerr << "GDML warning: cone '" << name
                              << "' has deltaphi=" << deltaphi
                              << " (not 2*pi); SIREN does not support partial cones,"
                              << " creating full cone" << std::endl;
                }
            }

            geo = Cone(rmin1, rmax1, rmin2, rmax2, hz * 2.0).create();
        }
        else if(tag == "trd") {
            // GDML <trd> uses half-widths; SIREN Trd also uses half-widths
            double dx1 = SafeParseDouble(SafeAttrVal(node, "x1"), data.constants) * lscale;
            double dx2 = SafeParseDouble(SafeAttrVal(node, "x2"), data.constants) * lscale;
            double dy1 = SafeParseDouble(SafeAttrVal(node, "y1"), data.constants) * lscale;
            double dy2 = SafeParseDouble(SafeAttrVal(node, "y2"), data.constants) * lscale;
            double dz  = SafeParseDouble(SafeAttrVal(node, "z"), data.constants) * lscale;
            geo = Trd(dx1, dx2, dy1, dy2, dz).create();
        }
        else if(tag == "polycone") {
            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi"), data.constants) * ascale;
            (void)startphi;

            // Check deltaphi for partial polycone
            const char* dp_val = SafeAttrVal(node, "deltaphi");
            if(dp_val[0] != '\0') {
                double deltaphi = SafeParseDouble(dp_val, data.constants) * ascale;
                if(std::fabs(deltaphi - 2.0 * PI) > ANG_TOL) {
                    std::cerr << "GDML warning: polycone '" << name
                              << "' has deltaphi=" << deltaphi
                              << " (not 2*pi); SIREN does not support partial polycones,"
                              << " creating full polycone" << std::endl;
                }
            }

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z"), data.constants) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin"), data.constants) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax"), data.constants) * lscale;
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
            double startphi = SafeParseDouble(SafeAttrVal(node, "startphi"), data.constants) * ascale;

            std::vector<double> z_planes;
            std::vector<double> rmin_vec;
            std::vector<double> rmax_vec;

            for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
                double z = SafeParseDouble(SafeAttrVal(zp, "z"), data.constants) * lscale;
                double rmin = SafeParseDouble(SafeAttrVal(zp, "rmin"), data.constants) * lscale;
                double rmax = SafeParseDouble(SafeAttrVal(zp, "rmax"), data.constants) * lscale;
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
                double vx = SafeParseDouble(SafeAttrVal(vtx, "x"), data.constants) * lscale;
                double vy = SafeParseDouble(SafeAttrVal(vtx, "y"), data.constants) * lscale;
                std::vector<double> vert;
                vert.push_back(vx);
                vert.push_back(vy);
                polygon.push_back(vert);
            }

            for(auto* sec = node->first_node("section"); sec; sec = sec->next_sibling("section")) {
                double zpos = SafeParseDouble(SafeAttrVal(sec, "zPosition"), data.constants) * lscale;
                double xoff = SafeParseDouble(SafeAttrVal(sec, "xOffset"), data.constants) * lscale;
                double yoff = SafeParseDouble(SafeAttrVal(sec, "yOffset"), data.constants) * lscale;
                double scale = SafeParseDouble(SafeAttrVal(sec, "scalingFactor"), data.constants);
                if(scale == 0.0) scale = 1.0;
                double offset[2] = {xoff, yoff};
                zsections.push_back(ExtrPoly::ZSection(zpos, offset, scale));
            }

            if(!polygon.empty() && !zsections.empty()) {
                geo = ExtrPoly(polygon, zsections).create();
            }
        }

        return geo;
    };

    // Helper lambda: parse a boolean solid node
    auto parseBoolean = [&](rapidxml::xml_node<>* node) -> std::shared_ptr<Geometry> {
        std::string tag(node->name());
        std::string name = SafeAttrVal(node, "name");

        BooleanOperation op;
        if(tag == "subtraction")   op = BooleanOperation::SUBTRACTION;
        else if(tag == "union")    op = BooleanOperation::UNION;
        else                       op = BooleanOperation::INTERSECTION;

        auto* first_node = node->first_node("first");
        auto* second_node = node->first_node("second");

        if(!first_node || !second_node) return nullptr;

        std::string first_ref = SafeAttrVal(first_node, "ref");
        std::string second_ref = SafeAttrVal(second_node, "ref");

        auto first_it = data.solids.find(first_ref);
        auto second_it = data.solids.find(second_ref);

        if(first_it == data.solids.end() || second_it == data.solids.end()) return nullptr;

        auto left = first_it->second;

        // The second solid may have a position and rotation
        // relative to the first
        Vector3D rel_pos(0, 0, 0);
        Quaternion rel_rot;

        // Check for inline <position> child
        auto* pos_node = node->first_node("position");
        if(pos_node) {
            const char* punit = SafeAttrVal(pos_node, "unit");
            double px = ParseLength(SafeAttrVal(pos_node, "x"), punit, data.constants);
            double py = ParseLength(SafeAttrVal(pos_node, "y"), punit, data.constants);
            double pz = ParseLength(SafeAttrVal(pos_node, "z"), punit, data.constants);
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
            double rx = ParseAngle(SafeAttrVal(rot_node, "x"), runit, data.constants);
            double ry = ParseAngle(SafeAttrVal(rot_node, "y"), runit, data.constants);
            double rz = ParseAngle(SafeAttrVal(rot_node, "z"), runit, data.constants);
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

        return std::make_shared<BooleanGeometry>(op,
            std::const_pointer_cast<const Geometry>(left),
            std::const_pointer_cast<const Geometry>(right));
    };

    // First pass: parse all non-boolean solids
    for(auto* node = solids_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());
        if(tag == "subtraction" || tag == "union" || tag == "intersection") continue;

        std::string name = SafeAttrVal(node, "name");
        if(name.empty()) continue;

        auto geo = parsePrimitive(node);
        if(geo) {
            data.solids[name] = geo;
        }
    }

    // Second pass: parse boolean solids (operands now all available)
    for(auto* node = solids_node->first_node(); node; node = node->next_sibling()) {
        std::string tag(node->name());
        if(tag != "subtraction" && tag != "union" && tag != "intersection") continue;

        std::string name = SafeAttrVal(node, "name");
        if(name.empty()) continue;

        auto geo = parseBoolean(node);
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
