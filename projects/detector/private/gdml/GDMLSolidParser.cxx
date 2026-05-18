#include "GDMLParserPrivate.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "SIREN/geometry/Placement.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/geometry/ExtrPoly.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Polyhedra.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/EllipticalTube.h"
#include "SIREN/geometry/CutTube.h"
#include "SIREN/geometry/Trap.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/Para.h"
#include "SIREN/geometry/GenericPolycone.h"
#include "SIREN/geometry/GeometryMesh.h"

using namespace siren::math;
using namespace siren::geometry;

namespace siren {
namespace detector {
namespace gdml {

namespace {

static constexpr double ANG_TOL = 1e-6;

struct SolidContext {
    rapidxml::xml_node<>* node;
    GDMLData & data;
    GDMLParseOptions const & options;
    std::string tag;
    std::string name;
    double lscale;
    double ascale;
};

SolidContext MakeSolidContext(rapidxml::xml_node<>* node,
                              GDMLData & data,
                              GDMLParseOptions const & options) {
    return SolidContext{
        node,
        data,
        options,
        std::string(node->name()),
        std::string(SafeAttrVal(node, "name")),
        LengthScale(SafeAttrVal(node, "lunit")),
        AngleScale(SafeAttrVal(node, "aunit"))
    };
}

double L(SolidContext const & ctx, const char* attr) {
    return ReadLengthAttr(ctx.node, attr, ctx.lscale, ctx.data.constants);
}

double A(SolidContext const & ctx, const char* attr) {
    return ReadAngleAttr(ctx.node, attr, ctx.ascale, ctx.data.constants);
}

double D(SolidContext const & ctx, const char* attr) {
    return ReadDoubleAttr(ctx.node, attr, ctx.data.constants);
}

bool HasPhiCut(GDMLPhiRange const & phi) {
    return std::fabs(phi.start) > ANG_TOL || std::fabs(phi.delta - 2.0 * M_PI) > ANG_TOL;
}

BooleanOperation BooleanOperationFromTag(std::string const & tag) {
    if(tag == "subtraction") return BooleanOperation::SUBTRACTION;
    if(tag == "union") return BooleanOperation::UNION;
    return BooleanOperation::INTERSECTION;
}

std::shared_ptr<Geometry> BuildBooleanGeometry(BooleanOperation op,
                                               std::shared_ptr<Geometry> const & left,
                                               std::shared_ptr<Geometry> const & right_base,
                                               GDMLPlacement const & first,
                                               GDMLPlacement const & second) {
    std::shared_ptr<Geometry> left_placed;
    if(first.specified) {
        left_placed = left->create();
        left_placed->SetPlacement(Placement(first.position, first.rotation));
    }

    auto right = right_base->create();
    right->SetPlacement(Placement(second.position, second.rotation));

    auto left_final = first.specified
        ? std::const_pointer_cast<const Geometry>(left_placed)
        : std::const_pointer_cast<const Geometry>(left);

    return std::make_shared<BooleanGeometry>(op, left_final,
        std::const_pointer_cast<const Geometry>(right));
}

std::shared_ptr<Geometry> ParseBox(SolidContext const & ctx) {
    // GDML <box> x,y,z are half-widths; SIREN Box takes full widths.
    return Box(L(ctx, "x") * 2.0, L(ctx, "y") * 2.0, L(ctx, "z") * 2.0).create();
}

std::shared_ptr<Geometry> ParseSphere(SolidContext const & ctx) {
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    double starttheta = A(ctx, "starttheta");
    double deltatheta = M_PI;
    const char* dt_val = SafeAttrVal(ctx.node, "deltatheta");
    if(dt_val[0] != '\0') {
        deltatheta = SafeParseDouble(dt_val, ctx.data.constants) * ctx.ascale;
    }

    return Sphere(L(ctx, "rmax"), L(ctx, "rmin"),
                  phi.start, phi.delta, starttheta, deltatheta).create();
}

std::shared_ptr<Geometry> ParseOrb(SolidContext const & ctx) {
    return Sphere(L(ctx, "r"), 0).create();
}

std::shared_ptr<Geometry> ParseTube(SolidContext const & ctx) {
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    // GDML <tube> z is half-height; SIREN Cylinder takes full height.
    return Cylinder(L(ctx, "rmax"), L(ctx, "rmin"), L(ctx, "z") * 2.0,
                    phi.start, phi.delta).create();
}

std::shared_ptr<Geometry> ParseCone(SolidContext const & ctx) {
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    // GDML <cone> z is half-height; SIREN Cone takes full height.
    return Cone(L(ctx, "rmin1"), L(ctx, "rmax1"),
                L(ctx, "rmin2"), L(ctx, "rmax2"),
                L(ctx, "z") * 2.0, phi.start, phi.delta).create();
}

std::shared_ptr<Geometry> ParseTrd(SolidContext const & ctx) {
    // GDML <trd> uses half-widths; SIREN Trd also uses half-widths.
    return Trd(L(ctx, "x1"), L(ctx, "x2"), L(ctx, "y1"), L(ctx, "y2"), L(ctx, "z")).create();
}

std::shared_ptr<Geometry> ParsePolycone(SolidContext const & ctx) {
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    GDMLZPlanes planes = ReadZPlanes(ctx.node, ctx.lscale, ctx.data);
    if(planes.z.size() < 2) return nullptr;

    if(!ValidateZPlaneMonotonicity(planes)) {
        throw std::runtime_error("polycone '" + ctx.name + "' has non-monotonic z-planes");
    }

    return Polycone(planes.z, planes.rmin, planes.rmax, phi.start, phi.delta).create();
}

std::shared_ptr<Geometry> ParsePolyhedra(SolidContext const & ctx) {
    int numSide = 0;
    std::string ns_str(SafeAttrVal(ctx.node, "numsides"));
    if(!ns_str.empty()) {
        try {
            numSide = std::stoi(ns_str);
        } catch(std::invalid_argument &) {
            numSide = 0;
        } catch(std::out_of_range &) {
            numSide = 0;
        }
    }
    if(numSide <= 0) return nullptr;

    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    GDMLZPlanes planes = ReadZPlanes(ctx.node, ctx.lscale, ctx.data);
    if(planes.z.size() < 2 || numSide < 3) return nullptr;

    if(!ValidateZPlaneMonotonicity(planes)) {
        throw std::runtime_error("polyhedra '" + ctx.name + "' has non-monotonic z-planes");
    }

    return Polyhedra(numSide, phi.start, planes.z, planes.rmin, planes.rmax, phi.delta).create();
}

std::shared_ptr<Geometry> ParseXtru(SolidContext const & ctx) {
    std::vector<std::vector<double>> polygon;
    std::vector<ExtrPoly::ZSection> zsections;

    for(auto* vtx = ctx.node->first_node("twoDimVertex"); vtx; vtx = vtx->next_sibling("twoDimVertex")) {
        polygon.push_back({
            ReadLengthAttr(vtx, "x", ctx.lscale, ctx.data.constants),
            ReadLengthAttr(vtx, "y", ctx.lscale, ctx.data.constants)
        });
    }

    for(auto* sec = ctx.node->first_node("section"); sec; sec = sec->next_sibling("section")) {
        double zpos = ReadLengthAttr(sec, "zPosition", ctx.lscale, ctx.data.constants);
        double xoff = ReadLengthAttr(sec, "xOffset", ctx.lscale, ctx.data.constants);
        double yoff = ReadLengthAttr(sec, "yOffset", ctx.lscale, ctx.data.constants);
        double scale = ReadDoubleAttr(sec, "scalingFactor", ctx.data.constants);
        if(scale == 0.0) scale = 1.0;
        double offset[2] = {xoff, yoff};
        zsections.push_back(ExtrPoly::ZSection(zpos, offset, scale));
    }

    if(polygon.empty() || zsections.empty()) return nullptr;
    return ExtrPoly(polygon, zsections).create();
}

std::shared_ptr<Geometry> ParseTorus(SolidContext const & ctx) {
    double rmin_t = L(ctx, "rmin");
    double rmax_t = L(ctx, "rmax");
    double rtor = L(ctx, "rtor");
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);

    if(rtor <= 0 || rmax_t <= 0) return nullptr;
    return Torus(rtor, rmax_t, rmin_t, phi.start, phi.delta).create();
}

std::shared_ptr<Geometry> ParseEllipticalTube(SolidContext const & ctx) {
    double dx = L(ctx, "dx");
    double dy = L(ctx, "dy");
    double dz = L(ctx, "dz");
    if(dx <= 0 || dy <= 0 || dz <= 0) return nullptr;
    return EllipticalTube(dx, dy, dz).create();
}

std::shared_ptr<Geometry> ParseCutTube(SolidContext const & ctx) {
    double rmin = L(ctx, "rmin");
    double rmax = L(ctx, "rmax");
    double hz = L(ctx, "z");
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);

    if(rmax <= 0 || hz <= 0) return nullptr;

    Vector3D low_norm(D(ctx, "lowX"), D(ctx, "lowY"), D(ctx, "lowZ"));
    Vector3D high_norm(D(ctx, "highX"), D(ctx, "highY"), D(ctx, "highZ"));
    return CutTube(rmin, rmax, hz, low_norm, high_norm, phi.start, phi.delta).create();
}

std::shared_ptr<Geometry> ParseTrap(SolidContext const & ctx) {
    return Trap(L(ctx, "z"), A(ctx, "theta"), A(ctx, "phi"),
                L(ctx, "y1"), L(ctx, "x1"), L(ctx, "x2"), A(ctx, "alpha1"),
                L(ctx, "y2"), L(ctx, "x3"), L(ctx, "x4"), A(ctx, "alpha2")).create();
}

std::shared_ptr<Geometry> ParseEllipsoid(SolidContext const & ctx) {
    double ax = L(ctx, "ax");
    double by = L(ctx, "by");
    double cz = L(ctx, "cz");
    if(ax <= 0 || by <= 0 || cz <= 0) return nullptr;

    double zcut1 = -cz;
    double zcut2 = cz;
    const char* z1_val = SafeAttrVal(ctx.node, "zcut1");
    if(z1_val[0] != '\0') zcut1 = SafeParseDouble(z1_val, ctx.data.constants) * ctx.lscale;
    const char* z2_val = SafeAttrVal(ctx.node, "zcut2");
    if(z2_val[0] != '\0') zcut2 = SafeParseDouble(z2_val, ctx.data.constants) * ctx.lscale;

    return Ellipsoid(ax, by, cz, zcut1, zcut2).create();
}

std::shared_ptr<Geometry> ParsePara(SolidContext const & ctx) {
    double dx = L(ctx, "x");
    double dy = L(ctx, "y");
    double dz = L(ctx, "z");
    if(dx <= 0 || dy <= 0 || dz <= 0) return nullptr;

    return Para(dx, dy, dz, A(ctx, "alpha"), A(ctx, "theta"), A(ctx, "phi")).create();
}

std::shared_ptr<Geometry> ParseGenericPolycone(SolidContext const & ctx) {
    GDMLPhiRange phi = ReadPhiRange(ctx.node, ctx.ascale, ctx.data.constants);
    std::vector<double> r_vec;
    std::vector<double> z_vec;

    for(auto* rz = ctx.node->first_node("rzpoint"); rz; rz = rz->next_sibling("rzpoint")) {
        r_vec.push_back(ReadLengthAttr(rz, "r", ctx.lscale, ctx.data.constants));
        z_vec.push_back(ReadLengthAttr(rz, "z", ctx.lscale, ctx.data.constants));
    }

    if(r_vec.size() < 3) return nullptr;
    if(HasPhiCut(phi)) return GenericPolycone(r_vec, z_vec, phi.start, phi.delta).create();
    return GenericPolycone(r_vec, z_vec).create();
}

std::shared_ptr<Geometry> ParseTessellated(SolidContext const & ctx) {
    std::vector<std::array<math::Vector3D, 3>> triangles;

    for(auto* facet = ctx.node->first_node(); facet; facet = facet->next_sibling()) {
        std::string ftag(facet->name());
        if(ftag == "triangular") {
            math::Vector3D v1 = ReadPositionReferenceAttr(facet, "vertex1", ctx.data);
            math::Vector3D v2 = ReadPositionReferenceAttr(facet, "vertex2", ctx.data);
            math::Vector3D v3 = ReadPositionReferenceAttr(facet, "vertex3", ctx.data);
            triangles.push_back({{v1, v2, v3}});
        } else if(ftag == "quadrangular") {
            math::Vector3D v1 = ReadPositionReferenceAttr(facet, "vertex1", ctx.data);
            math::Vector3D v2 = ReadPositionReferenceAttr(facet, "vertex2", ctx.data);
            math::Vector3D v3 = ReadPositionReferenceAttr(facet, "vertex3", ctx.data);
            math::Vector3D v4 = ReadPositionReferenceAttr(facet, "vertex4", ctx.data);
            triangles.push_back({{v1, v2, v3}});
            triangles.push_back({{v1, v3, v4}});
        }
    }

    if(triangles.empty()) return nullptr;
    return TriangularMesh(triangles).create();
}

std::shared_ptr<Geometry> ParseArb8(SolidContext const & ctx) {
    double dz = L(ctx, "dz");
    if(dz <= 0) return nullptr;

    static const char* x_attrs[8] = {"v1x", "v2x", "v3x", "v4x", "v5x", "v6x", "v7x", "v8x"};
    static const char* y_attrs[8] = {"v1y", "v2y", "v3y", "v4y", "v5y", "v6y", "v7y", "v8y"};

    math::Vector3D verts[8];
    for(int i = 0; i < 8; ++i) {
        double z = (i < 4) ? -dz : dz;
        verts[i] = math::Vector3D(
            ReadLengthAttr(ctx.node, x_attrs[i], ctx.lscale, ctx.data.constants),
            ReadLengthAttr(ctx.node, y_attrs[i], ctx.lscale, ctx.data.constants),
            z);
    }

    std::vector<std::array<math::Vector3D, 3>> triangles;
    triangles.reserve(12);

    auto addQuad = [&](math::Vector3D const & a, math::Vector3D const & b,
                       math::Vector3D const & c, math::Vector3D const & d) {
        triangles.push_back({{a, b, c}});
        triangles.push_back({{a, c, d}});
    };

    // Bottom face (outward normal in -z): v4, v3, v2, v1.
    addQuad(verts[3], verts[2], verts[1], verts[0]);
    // Top face (outward normal in +z): v5, v6, v7, v8.
    addQuad(verts[4], verts[5], verts[6], verts[7]);
    addQuad(verts[0], verts[1], verts[5], verts[4]);
    addQuad(verts[1], verts[2], verts[6], verts[5]);
    addQuad(verts[2], verts[3], verts[7], verts[6]);
    addQuad(verts[3], verts[0], verts[4], verts[7]);

    return TriangularMesh(triangles).create();
}

std::shared_ptr<Geometry> ParsePrimitiveSolid(SolidContext const & ctx) {
    if(ctx.tag == "box") return ParseBox(ctx);
    if(ctx.tag == "sphere") return ParseSphere(ctx);
    if(ctx.tag == "orb") return ParseOrb(ctx);
    if(ctx.tag == "tube" || ctx.tag == "tubs") return ParseTube(ctx);
    if(ctx.tag == "cone") return ParseCone(ctx);
    if(ctx.tag == "trd") return ParseTrd(ctx);
    if(ctx.tag == "polycone") return ParsePolycone(ctx);
    if(ctx.tag == "polyhedra") return ParsePolyhedra(ctx);
    if(ctx.tag == "xtru") return ParseXtru(ctx);
    if(ctx.tag == "torus") return ParseTorus(ctx);
    if(ctx.tag == "eltube") return ParseEllipticalTube(ctx);
    if(ctx.tag == "cutTube") return ParseCutTube(ctx);
    if(ctx.tag == "trap") return ParseTrap(ctx);
    if(ctx.tag == "ellipsoid") return ParseEllipsoid(ctx);
    if(ctx.tag == "para") return ParsePara(ctx);
    if(ctx.tag == "genericPolycone") return ParseGenericPolycone(ctx);
    if(ctx.tag == "tessellated") return ParseTessellated(ctx);
    if(ctx.tag == "arb8") return ParseArb8(ctx);

    if(!ctx.tag.empty()) {
        EmitWarning(ctx.data, ctx.options,
            "unsupported solid type '" + ctx.tag + "' (name='" + ctx.name + "'), skipping");
    }
    return nullptr;
}

} // namespace

// Parse all <solids> sections: convert GDML solids to SIREN Geometry objects.
// Iterates all <solids> children of root_node, running the primitive parse pass
// for each section first, then performing iterative boolean resolution across
// all sections.
void ParseAllSolids(rapidxml::xml_node<>* root_node, GDMLData & data, GDMLParseOptions const & options) {
    if(!root_node) return;

    std::vector<rapidxml::xml_node<>*> solids_sections;
    for(auto* node = root_node->first_node("solids"); node; node = node->next_sibling("solids")) {
        solids_sections.push_back(node);
    }
    if(solids_sections.empty()) return;

    // Instance-aware storage: each name can have multiple definitions.
    // solids_instances[name] = vector of shared_ptr (one per definition)
    // solid_instance_count[name] = how many definitions seen so far
    std::map<std::string, std::vector<std::shared_ptr<Geometry>>> solids_instances;
    std::map<std::string, int> solid_instance_count;

    auto storeSolid = [&](std::string const & name, std::shared_ptr<Geometry> geo) {
        auto & vec = solids_instances[name];
        if(!vec.empty() && *vec.back() == *geo) {
            return;
        }
        vec.push_back(geo);
        solid_instance_count[name] = (int)vec.size();
    };

    auto lookupSolid = [&](std::string const & name) -> std::shared_ptr<Geometry> {
        auto it = solid_instance_count.find(name);
        if(it != solid_instance_count.end()) {
            return solids_instances[name][it->second - 1];
        }
        return nullptr;
    };

    auto parseBooleanInstanced = [&](rapidxml::xml_node<>* node) -> std::shared_ptr<Geometry> {
        BooleanOperation op = BooleanOperationFromTag(std::string(node->name()));

        auto* first_node = node->first_node("first");
        auto* second_node = node->first_node("second");
        if(!first_node || !second_node) return nullptr;

        auto left = lookupSolid(SafeAttrVal(first_node, "ref"));
        auto right_base = lookupSolid(SafeAttrVal(second_node, "ref"));
        if(!left || !right_base) return nullptr;

        GDMLPlacement first = ReadPlacement(
            node, "firstposition", "firstpositionref", "firstrotation", "firstrotationref", data);
        GDMLPlacement second = ReadPlacement(
            node, "position", "positionref", "rotation", "rotationref", data);
        return BuildBooleanGeometry(op, left, right_base, first, second);
    };

    struct DeferredBoolean {
        rapidxml::xml_node<>* node;
        std::string name;
        std::string first_ref;
        std::string second_ref;
        int first_instance;
        int second_instance;
    };

    std::vector<DeferredBoolean> deferred;

    for(auto* solids_node : solids_sections) {
        for(auto* node = solids_node->first_node(); node; node = node->next_sibling()) {
            std::string tag(node->name());
            std::string name = SafeAttrVal(node, "name");
            if(name.empty()) continue;

            bool is_boolean = (tag == "subtraction" || tag == "union" || tag == "intersection");

            if(!is_boolean) {
                auto geo = ParsePrimitiveSolid(MakeSolidContext(node, data, options));
                if(geo) {
                    storeSolid(name, geo);
                }
                continue;
            }

            auto* first_child = node->first_node("first");
            auto* second_child = node->first_node("second");
            if(!first_child || !second_child) continue;

            std::string first_ref = SafeAttrVal(first_child, "ref");
            std::string second_ref = SafeAttrVal(second_child, "ref");

            auto left = lookupSolid(first_ref);
            auto right = lookupSolid(second_ref);

            if(left && right) {
                auto geo = parseBooleanInstanced(node);
                if(geo) {
                    storeSolid(name, geo);
                }
                continue;
            }

            DeferredBoolean db;
            db.node = node;
            db.name = name;
            db.first_ref = first_ref;
            db.second_ref = second_ref;
            db.first_instance = left ? (solid_instance_count[first_ref] - 1) : 0;
            db.second_instance = right ? (solid_instance_count[second_ref] - 1) : 0;
            deferred.push_back(db);
        }
    }

    bool made_progress = true;
    while(made_progress) {
        made_progress = false;
        for(auto & db : deferred) {
            if(solid_instance_count.find(db.name) != solid_instance_count.end()) {
                if(db.node == nullptr) continue;
            }

            auto first_it = solids_instances.find(db.first_ref);
            auto second_it = solids_instances.find(db.second_ref);

            if(first_it == solids_instances.end() || second_it == solids_instances.end()) continue;
            if(db.first_instance >= (int)first_it->second.size()) continue;
            if(db.second_instance >= (int)second_it->second.size()) continue;

            auto left = first_it->second[db.first_instance];
            auto right_base = second_it->second[db.second_instance];
            if(!left || !right_base) continue;

            BooleanOperation op = BooleanOperationFromTag(std::string(db.node->name()));
            GDMLPlacement first = ReadPlacement(
                db.node, "firstposition", "firstpositionref", "firstrotation", "firstrotationref", data);
            GDMLPlacement second = ReadPlacement(
                db.node, "position", "positionref", "rotation", "rotationref", data);
            auto geo = BuildBooleanGeometry(op, left, right_base, first, second);

            storeSolid(db.name, geo);
            db.node = nullptr;
            made_progress = true;
        }
    }

    for(auto const & db : deferred) {
        if(db.node != nullptr) {
            EmitWarning(data, options, "boolean solid '" + db.name + "' has unresolved operand references; skipping");
        }
    }

    for(auto const & pair : solids_instances) {
        std::string const & name = pair.first;
        std::vector<std::shared_ptr<Geometry>> const & instances = pair.second;
        data.solid_instance_counts[name] = (int)instances.size();
        for(int i = 0; i < (int)instances.size(); ++i) {
            std::string flat_name;
            if(i == 0) {
                flat_name = name;
            } else {
                flat_name = name + "__" + std::to_string(i + 1);
            }
            data.solids[flat_name] = instances[i];
        }
    }
}

} // namespace gdml
} // namespace detector
} // namespace siren
