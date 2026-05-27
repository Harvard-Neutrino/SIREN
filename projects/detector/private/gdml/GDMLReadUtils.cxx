#include "GDMLParserPrivate.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace siren::math;

namespace siren {
namespace detector {
namespace gdml {

double ReadDoubleAttr(rapidxml::xml_node<>* node, const char* attr,
                      std::unordered_map<std::string, double> const & constants,
                      double scale) {
    return SafeParseDouble(SafeAttrVal(node, attr), constants) * scale;
}

double ReadLengthAttr(rapidxml::xml_node<>* node, const char* attr,
                      double lscale,
                      std::unordered_map<std::string, double> const & constants) {
    return ReadDoubleAttr(node, attr, constants, lscale);
}

double ReadAngleAttr(rapidxml::xml_node<>* node, const char* attr,
                     double ascale,
                     std::unordered_map<std::string, double> const & constants) {
    return ReadDoubleAttr(node, attr, constants, ascale);
}

GDMLPhiRange ReadPhiRange(rapidxml::xml_node<>* node,
                          double ascale,
                          std::unordered_map<std::string, double> const & constants) {
    GDMLPhiRange range;
    range.start = 0.0;
    range.delta = 2.0 * M_PI;

    const char* start = SafeAttrVal(node, "startphi");
    if(start[0] != '\0') {
        range.start = SafeParseDouble(start, constants) * ascale;
    }

    const char* delta = SafeAttrVal(node, "deltaphi");
    if(delta[0] != '\0') {
        range.delta = SafeParseDouble(delta, constants) * ascale;
    }

    return range;
}

Vector3D ReadPositionXYZ(rapidxml::xml_node<>* node, GDMLData const & data) {
    if(!node) return Vector3D(0, 0, 0);

    const char* unit = SafeAttrVal(node, "unit");
    double x = ParseLength(SafeAttrVal(node, "x"), unit, data.constants);
    double y = ParseLength(SafeAttrVal(node, "y"), unit, data.constants);
    double z = ParseLength(SafeAttrVal(node, "z"), unit, data.constants);
    return Vector3D(x, y, z);
}

Quaternion ReadRotationXYZ(rapidxml::xml_node<>* node, GDMLData const & data) {
    if(!node) return Quaternion();

    const char* unit = SafeAttrVal(node, "unit");
    double x = ParseAngle(SafeAttrVal(node, "x"), unit, data.constants);
    double y = ParseAngle(SafeAttrVal(node, "y"), unit, data.constants);
    double z = ParseAngle(SafeAttrVal(node, "z"), unit, data.constants);
    return QuatFromGDMLRotation(x, y, z);
}

Vector3D ReadPositionReferenceAttr(rapidxml::xml_node<>* node,
                                  const char* attr,
                                  GDMLData const & data) {
    std::string ref = SafeAttrVal(node, attr);
    if(ref.empty()) {
        return Vector3D(0, 0, 0);
    }
    auto it = data.positions.find(ref);
    if(it != data.positions.end()) {
        return it->second;
    }
    throw std::runtime_error("GDML error: missing position reference '" + ref + "'");
}

GDMLPlacement ReadPlacement(rapidxml::xml_node<>* node,
                            const char* position_tag,
                            const char* position_ref_tag,
                            const char* rotation_tag,
                            const char* rotation_ref_tag,
                            GDMLData const & data) {
    GDMLPlacement placement;

    auto* pos_node = node->first_node(position_tag);
    if(pos_node) {
        placement.position = ReadPositionXYZ(pos_node, data);
        placement.specified = true;
    }

    auto* posref_node = node->first_node(position_ref_tag);
    if(posref_node) {
        std::string ref = SafeAttrVal(posref_node, "ref");
        auto it = data.positions.find(ref);
        if(it == data.positions.end()) {
            throw std::runtime_error("GDML error: unresolved position reference '" + ref + "'");
        }
        placement.position = it->second;
        placement.specified = true;
    }

    auto* rot_node = node->first_node(rotation_tag);
    if(rot_node) {
        placement.rotation = ReadRotationXYZ(rot_node, data);
        placement.specified = true;
    }

    auto* rotref_node = node->first_node(rotation_ref_tag);
    if(rotref_node) {
        std::string ref = SafeAttrVal(rotref_node, "ref");
        auto it = data.rotations.find(ref);
        if(it == data.rotations.end()) {
            throw std::runtime_error("GDML error: unresolved rotation reference '" + ref + "'");
        }
        placement.rotation = it->second;
        placement.specified = true;
    }

    return placement;
}

GDMLZPlanes ReadZPlanes(rapidxml::xml_node<>* node,
                        double lscale,
                        GDMLData const & data) {
    GDMLZPlanes planes;

    for(auto* zp = node->first_node("zplane"); zp; zp = zp->next_sibling("zplane")) {
        planes.z.push_back(ReadLengthAttr(zp, "z", lscale, data.constants));
        planes.rmin.push_back(ReadLengthAttr(zp, "rmin", lscale, data.constants));
        planes.rmax.push_back(ReadLengthAttr(zp, "rmax", lscale, data.constants));
    }

    return planes;
}

bool ValidateZPlaneMonotonicity(GDMLZPlanes const & planes) {
    if(planes.z.size() < 2) return false;
    // Z-planes must be monotonic (ascending or descending, duplicates allowed).
    // The geometry constructors handle reversal to ascending order internally.
    bool ascending = true, descending = true;
    for(size_t i = 1; i < planes.z.size(); ++i) {
        if(planes.z[i] < planes.z[i - 1]) ascending = false;
        if(planes.z[i] > planes.z[i - 1]) descending = false;
    }
    return ascending || descending;
}

} // namespace gdml
} // namespace detector
} // namespace siren
