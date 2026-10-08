#include "SIREN/injection/GeometryVolume.h"

#include "SIREN/geometry/BooleanGeometry.h"
#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Cone.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/Ellipsoid.h"
#include "SIREN/geometry/EllipticalTube.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/GenericPolycone.h"
#include "SIREN/geometry/GeometryMesh.h"
#include "SIREN/geometry/Para.h"
#include "SIREN/geometry/Polycone.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Torus.h"
#include "SIREN/geometry/Trd.h"
#include "SIREN/math/Vector3D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace siren {
namespace injection {

namespace {

constexpr double kTwoPi = 2.0 * M_PI;
constexpr double kMinimumViableFillFraction = 1e-4;
// Check of a supplied volume by chord integration. Rays along kDirections fixed
// directions, spread over a hemisphere (a Fibonacci lattice turned by fixed
// skew angles, so none lies along a coordinate axis), cross the bounding box on
// a jittered kRaysPerSide x kRaysPerSide grid over its projection. Each
// direction's estimate is the grid's area times the mean length of the rays
// inside the solid; a thin shell or slab gives every crossing ray a short
// chord, so it is resolved where counting inside points is not. The volume is
// the plain mean over directions and its standard error their scatter over
// sqrt(kDirections): a direction running parallel to a thin wall misses it (or
// meets it in a rare long chord), which widens the error instead of biasing a
// tight result, and a solid cannot line up its walls with more than a few of
// the directions.
constexpr int kDirections = 48;
constexpr int kRaysPerSide = 48;
constexpr double kSkewAngles[3] = {0.4, 0.7, 0.9};
// At most kSparseDirections directions may have fewer than kResolvedHits rays
// that meet the solid (a direction nearly in the plane of a thin sheet sees it
// edge-on); each direction's estimate is unbiased on its own, so such a
// direction stays in the mean and its noise in the scatter. More sparse
// directions mean the solid is too small or thin for the grid.
constexpr int kResolvedHits = 20;
constexpr int kSparseDirections = kDirections / 10;
// A larger relative standard error makes the accepted range too wide to check.
constexpr double kMaximumRelativeError = 0.05;
// Tolerated fraction of rays whose crossings do not alternate between entering
// and exiting, and of inside-test cross-checks that disagree (numerical edge
// and corner hits); more means the surface is not a well-formed solid.
constexpr double kMalformedFraction = 1e-3;
// Crossings closer than Geometry::GEOMETRY_PRECISION are merged by
// BooleanGeometry, so chords within a few times that may have lost material.
constexpr double kChordResolution = 5.0 * siren::geometry::Geometry::GEOMETRY_PRECISION;
constexpr double kSubResolutionFraction = 0.01;
// Inside test at interval midpoints on the rays of every kInsideCheckStride-th
// grid row and column that meet the solid, for intervals at least
// kInsideCheckLength long.
constexpr int kInsideCheckStride = 4;
constexpr double kInsideCheckLength = 1e3 * siren::geometry::Geometry::GEOMETRY_PRECISION;
// Rays can only weigh the parts of a solid they meet: a small separate part (a
// core inside a thin shell) can hold most of the volume while every direction
// is resolved by the rest. The parts are the primitives of a Boolean tree, the
// connected pieces of a triangle mesh, and the tree's intersection and
// subtraction nodes (which can leave a small solid out of broad operands). A
// primitive with an exact volume inside the solid's bounding box is estimated
// by the same rays on its own: the part of its exact volume the estimate does
// not reach (beyond five standard errors) may have been missed. Any other part
// needs kPartResolvedHits meeting rays, or its bounding box bounds what may
// have been missed (a part met by a few rays shows in the scatter instead).
// Together these must stay below kUnseenFraction of the estimate.
constexpr long long kPartResolvedHits = 10;
constexpr double kUnseenFraction = 0.01;
// The channel samples points with IsInside, which looks for the first
// crossing along world +z beyond GEOMETRY_PRECISION: the last
// GEOMETRY_PRECISION before each exit is sampled as outside. That layer must
// hold less than kInsideLayerFraction of the volume, or the channel's samples
// do not fill the solid its density is normalized to (a solid a few nm thick).
constexpr double kInsideLayerFraction = 0.005;
// A supplied volume may differ from the estimate by 2% or five standard
// errors, whichever is larger.
constexpr double kEstimateRelativeTolerance = 0.02;
constexpr double kEstimateStandardErrors = 5.0;

double AABBVolume(siren::geometry::Geometry const & geometry) {
    auto aabb = geometry.GetWorldBoundingBox();
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    return extent.GetX() * extent.GetY() * extent.GetZ();
}

// Relative rounding of a bounding-box volume computed from world corners: each
// extent hi - lo carries about eps*(|lo| + |hi|).
double AABBRounding(siren::geometry::Geometry const & geometry) {
    auto aabb = geometry.GetWorldBoundingBox();
    double eps = std::numeric_limits<double>::epsilon(), rounding = 1e-12;
    double lo[3] = {aabb.min_corner.GetX(), aabb.min_corner.GetY(), aabb.min_corner.GetZ()};
    double hi[3] = {aabb.max_corner.GetX(), aabb.max_corner.GetY(), aabb.max_corner.GetZ()};
    for (int i = 0; i < 3; ++i)
        if (hi[i] > lo[i]) rounding += 8 * eps * (std::abs(lo[i]) + std::abs(hi[i])) / (hi[i] - lo[i]);
    return rounding;
}

double AzimuthalExtent(double delta_phi) {
    return std::clamp(delta_phi, 0.0, kTwoPi);
}

std::string Format(double value) {
    std::ostringstream stream;
    stream.precision(10);
    stream << value;
    return stream.str();
}

// A part of the solid, tested in the frame its Intersections works in: world
// coordinates taken through its Boolean ancestors' GlobalToLocal transforms.
struct SolidPart {
    std::vector<siren::geometry::Geometry const *> ancestors;
    std::shared_ptr<siren::geometry::Geometry const> owned;  // a mesh piece built here
    siren::geometry::Geometry const * geometry = nullptr;
    std::string name;
    double lo[3] = {0.0, 0.0, 0.0}, hi[3] = {0.0, 0.0, 0.0};  // world bounding box
    double box = 0.0;    // its volume within the solid's bounding box
    double exact = std::numeric_limits<double>::quiet_NaN();  // compared with its own estimate
    long long hits = 0;
    double sum = 0.0;            // chords through it along the current direction
    std::vector<double> volumes; // its estimate along each direction
};

// Connected pieces of a mesh: triangles sharing a vertex belong together.
std::vector<std::vector<std::array<siren::math::Vector3D, 3>>> MeshPieces(
    siren::geometry::TriangularMesh const & mesh)
{
    using siren::math::Vector3D;
    auto const & triangles = mesh.GetTriangles();
    std::vector<std::size_t> root(triangles.size());
    std::iota(root.begin(), root.end(), std::size_t(0));
    auto find = [&root](std::size_t i) {
        while (root[i] != i) i = root[i] = root[root[i]];
        return i;
    };
    std::map<std::array<double, 3>, std::size_t> first;
    for (std::size_t t = 0; t < triangles.size(); ++t) {
        for (auto const & vertex : triangles[t]) {
            auto found = first.emplace(vertex, t);
            if (!found.second) root[find(t)] = find(found.first->second);
        }
    }
    std::map<std::size_t, std::vector<std::array<Vector3D, 3>>> pieces;
    for (std::size_t t = 0; t < triangles.size(); ++t) {
        auto const & tri = triangles[t];
        pieces[find(t)].push_back({Vector3D(tri[0][0], tri[0][1], tri[0][2]),
                                   Vector3D(tri[1][0], tri[1][1], tri[1][2]),
                                   Vector3D(tri[2][0], tri[2][1], tri[2][2])});
    }
    std::vector<std::vector<std::array<Vector3D, 3>>> result;
    for (auto & piece : pieces) result.push_back(std::move(piece.second));
    return result;
}

void CollectParts(siren::geometry::Geometry const & geometry,
                  std::vector<siren::geometry::Geometry const *> ancestors,
                  std::vector<SolidPart> & parts)
{
    if (auto const * boolean = dynamic_cast<siren::geometry::BooleanGeometry const *>(&geometry)) {
        if (!ancestors.empty() && boolean->GetOperation() != siren::geometry::BooleanOperation::UNION) {
            SolidPart node;
            node.ancestors = ancestors;
            node.geometry = &geometry;
            node.name = boolean->GetOperation() == siren::geometry::BooleanOperation::INTERSECTION
                ? "an intersection in its Boolean tree" : "a subtraction in its Boolean tree";
            parts.push_back(std::move(node));
        }
        ancestors.push_back(&geometry);
        for (auto const & child : {boolean->GetLeft(), boolean->GetRight()})
            if (child) CollectParts(*child, ancestors, parts);
        return;
    }
    if (auto const * mesh = dynamic_cast<siren::geometry::TriangularMesh const *>(&geometry)) {
        auto pieces = MeshPieces(*mesh);
        if (pieces.size() > 1) {
            for (auto & piece : pieces) {
                SolidPart part;
                part.ancestors = ancestors;
                part.name = "a " + std::to_string(piece.size()) + "-triangle piece of a mesh";
                part.owned = std::make_shared<siren::geometry::TriangularMesh>(mesh->GetPlacement(), piece);
                part.geometry = part.owned.get();
                parts.push_back(std::move(part));
            }
            return;
        }
    }
    SolidPart part;
    part.ancestors = ancestors;
    part.geometry = &geometry;
    part.name = "a " + geometry.GetName();
    parts.push_back(std::move(part));
}

// Whether a mesh's shells form a solid: split into edge-connected shells,
// every shell wound inward (negative signed volume) must be a cavity inside an
// outward shell, and no outward shell may lie inside another. Otherwise the
// divergence sum and the inside test count different sets (an inward shell
// outside everything makes IsInside count the column below it as inside, and
// nested outward shells are counted twice), wherever the rays happen to go.
// Containment is tested at a vertex-averaged point of the shell.
bool MeshShellsValid(siren::geometry::TriangularMesh const & mesh) {
    using siren::math::Vector3D;
    auto const & triangles = mesh.GetTriangles();
    std::vector<std::size_t> root(triangles.size());
    std::iota(root.begin(), root.end(), std::size_t(0));
    auto find = [&root](std::size_t i) {
        while (root[i] != i) i = root[i] = root[root[i]];
        return i;
    };
    std::map<std::array<double, 6>, std::size_t> edges;
    for (std::size_t t = 0; t < triangles.size(); ++t) {
        for (int k = 0; k < 3; ++k) {
            auto const & a = triangles[t][k];
            auto const & b = triangles[t][(k + 1) % 3];
            std::array<double, 6> key = a < b ? std::array<double, 6>{a[0], a[1], a[2], b[0], b[1], b[2]}
                                              : std::array<double, 6>{b[0], b[1], b[2], a[0], a[1], a[2]};
            auto found = edges.emplace(key, t);
            if (!found.second) root[find(t)] = find(found.first->second);
        }
    }
    struct Shell {
        std::vector<std::array<Vector3D, 3>> triangles;
        double signed_volume = 0.0;
        double point[3] = {0.0, 0.0, 0.0};
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    };
    std::map<std::size_t, Shell> shells;
    for (std::size_t t = 0; t < triangles.size(); ++t) {
        auto const & v = triangles[t];
        Shell & shell = shells[find(t)];
        shell.triangles.push_back({Vector3D(v[0][0], v[0][1], v[0][2]), Vector3D(v[1][0], v[1][1], v[1][2]),
                                   Vector3D(v[2][0], v[2][1], v[2][2])});
        for (int k = 0; k < 3; ++k)
            for (int m = 0; m < 3; ++m) {
                shell.point[m] += v[k][m];
                shell.lo[m] = std::min(shell.lo[m], v[k][m]);
                shell.hi[m] = std::max(shell.hi[m], v[k][m]);
            }
    }
    for (auto & entry : shells) {
        Shell & shell = entry.second;
        double c[3];
        for (int m = 0; m < 3; ++m) { c[m] = 0.5 * (shell.lo[m] + shell.hi[m]); shell.point[m] /= 3.0 * shell.triangles.size(); }
        for (auto const & v : shell.triangles) {
            double a[3] = {v[0].GetX() - c[0], v[0].GetY() - c[1], v[0].GetZ() - c[2]};
            double b[3] = {v[1].GetX() - v[0].GetX(), v[1].GetY() - v[0].GetY(), v[1].GetZ() - v[0].GetZ()};
            double e[3] = {v[2].GetX() - v[0].GetX(), v[2].GetY() - v[0].GetY(), v[2].GetZ() - v[0].GetZ()};
            shell.signed_volume += a[0] * (b[1] * e[2] - b[2] * e[1]) - a[1] * (b[0] * e[2] - b[2] * e[0])
                                 + a[2] * (b[0] * e[1] - b[1] * e[0]);
        }
    }
    if (shells.size() == 1) return shells.begin()->second.signed_volume > 0.0;
    std::vector<Shell const *> outward;
    std::vector<std::unique_ptr<siren::geometry::TriangularMesh>> solids;
    for (auto const & entry : shells) {
        if (entry.second.signed_volume > 0.0) {
            outward.push_back(&entry.second);
            solids.push_back(std::make_unique<siren::geometry::TriangularMesh>(entry.second.triangles));
        } else if (!(entry.second.signed_volume < 0.0)) {
            return false;
        }
    }
    auto contained = [&](Shell const & shell, Shell const * self) {
        for (std::size_t k = 0; k < outward.size(); ++k) {
            Shell const & other = *outward[k];
            if (&other == self) continue;
            bool boxed = true;
            for (int m = 0; m < 3; ++m) boxed = boxed && other.lo[m] <= shell.point[m] && shell.point[m] <= other.hi[m];
            if (boxed && solids[k]->IsInside(Vector3D(shell.point[0], shell.point[1], shell.point[2]))) return true;
        }
        return false;
    };
    for (auto const & entry : shells) {
        Shell const & shell = entry.second;
        bool inside = contained(shell, &shell);
        if (shell.signed_volume < 0.0 ? !inside : inside) return false;
    }
    return true;
}

// Whether the closed polygon (r[i], z[i]) is simple: after dropping repeated
// consecutive vertices, no two edges meet except adjacent ones at their shared
// vertex. The first-moment volume of a revolved profile needs this; a profile
// traced twice, a bow tie or touching loops would be counted wrongly.
bool SimplePolygon(std::vector<double> const & r, std::vector<double> const & z) {
    std::vector<std::array<double, 2>> p;
    for (std::size_t i = 0; i < r.size() && i < z.size(); ++i) {
        std::array<double, 2> v{r[i], z[i]};
        if (p.empty() || p.back() != v) p.push_back(v);
    }
    while (p.size() > 1 && p.front() == p.back()) p.pop_back();
    std::size_t n = p.size();
    if (n < 3) return false;
    // Orientation with its rounding bound: within the bound the three points
    // count as collinear, so a touch hidden by rounding is still found (a
    // conservative answer: the profile then goes to the ray estimate).
    auto orient = [](std::array<double, 2> const & a, std::array<double, 2> const & b, std::array<double, 2> const & c) {
        double left = (b[0] - a[0]) * (c[1] - a[1]), right = (b[1] - a[1]) * (c[0] - a[0]);
        double det = left - right;
        double bound = 8.0 * std::numeric_limits<double>::epsilon() * (std::abs(left) + std::abs(right))
            + 4.0 * std::numeric_limits<double>::epsilon()
              * (std::abs(b[0] - a[0]) + std::abs(b[1] - a[1])) * (std::abs(c[0] - a[0]) + std::abs(c[1] - a[1]));
        return std::abs(det) <= bound ? 0.0 : det;
    };
    auto within = [](std::array<double, 2> const & a, std::array<double, 2> const & b, std::array<double, 2> const & c) {
        return std::min(a[0], b[0]) <= c[0] && c[0] <= std::max(a[0], b[0])
            && std::min(a[1], b[1]) <= c[1] && c[1] <= std::max(a[1], b[1]);
    };
    auto meet = [&](std::array<double, 2> const & a, std::array<double, 2> const & b,
                    std::array<double, 2> const & c, std::array<double, 2> const & d) {
        double d1 = orient(c, d, a), d2 = orient(c, d, b), d3 = orient(a, b, c), d4 = orient(a, b, d);
        if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return true;
        return (d1 == 0 && within(c, d, a)) || (d2 == 0 && within(c, d, b))
            || (d3 == 0 && within(a, b, c)) || (d4 == 0 && within(a, b, d));
    };
    for (std::size_t i = 0; i < n; ++i) {
        auto const & a = p[i];
        auto const & b = p[(i + 1) % n];
        for (std::size_t j = i + 1; j < n; ++j) {
            auto const & c = p[j];
            auto const & d = p[(j + 1) % n];
            // Adjacent edges share a vertex; collinear ones must not turn back
            // along each other (a spike, or a profile retracing itself).
            auto folds = [&](std::array<double, 2> const & shared, std::array<double, 2> const & p,
                             std::array<double, 2> const & q) {
                return orient(p, shared, q) == 0
                    && (p[0] - shared[0]) * (q[0] - shared[0]) + (p[1] - shared[1]) * (q[1] - shared[1]) > 0;
            };
            if (j == i + 1) {
                if (folds(b, a, d)) return false;
                continue;
            }
            if (i == 0 && j == n - 1) {
                if (folds(a, b, c)) return false;
                continue;
            }
            if (meet(a, b, c, d)) return false;
        }
    }
    return true;
}

// Total length inside a solid of the ray start + s d (s > 0), from its sorted
// crossings; zero unless they alternate entering and exiting.
double ChordLength(std::vector<siren::geometry::Geometry::Intersection> crossings) {
    std::sort(crossings.begin(), crossings.end(),
        [](auto const & a, auto const & b) { return a.distance < b.distance; });
    if (crossings.size() % 2 != 0) return 0.0;
    double chord = 0.0;
    for (std::size_t c = 0; c + 1 < crossings.size(); c += 2) {
        if (!crossings[c].entering || crossings[c + 1].entering) return 0.0;
        chord += crossings[c + 1].distance - crossings[c].distance;
    }
    return std::isfinite(chord) ? chord : 0.0;
}

// The triangle meshes of a solid, including those in a Boolean tree.
std::vector<siren::geometry::TriangularMesh const *> MeshesIn(siren::geometry::Geometry const & geometry) {
    if (auto const * mesh = dynamic_cast<siren::geometry::TriangularMesh const *>(&geometry)) return {mesh};
    std::vector<siren::geometry::TriangularMesh const *> result;
    if (auto const * boolean = dynamic_cast<siren::geometry::BooleanGeometry const *>(&geometry)) {
        for (auto const & child : {boolean->GetLeft(), boolean->GetRight()}) {
            if (!child) continue;
            auto more = MeshesIn(*child);
            result.insert(result.end(), more.begin(), more.end());
        }
    }
    return result;
}

// Whether the line start + s d (s > 0) crosses the box [lo, hi].
bool RayMeetsBox(siren::math::Vector3D const & start, siren::math::Vector3D const & d,
                 double const lo[3], double const hi[3])
{
    double s[3] = {start.GetX(), start.GetY(), start.GetZ()}, v[3] = {d.GetX(), d.GetY(), d.GetZ()};
    double enter = 0.0, leave = std::numeric_limits<double>::infinity();
    for (int m = 0; m < 3; ++m) {
        if (v[m] == 0.0) {
            if (s[m] < lo[m] || s[m] > hi[m]) return false;
            continue;
        }
        double a = (lo[m] - s[m]) / v[m], b = (hi[m] - s[m]) / v[m];
        enter = std::max(enter, std::min(a, b));
        leave = std::min(leave, std::max(a, b));
    }
    return enter <= leave;
}

// SplitMix64: a fixed, portable sequence, so the estimate is reproducible.
struct DeterministicUniform {
    std::uint64_t state = 0x5eed5eed5eed5eedULL;
    double operator()() {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
        return (z >> 11) * (1.0 / 9007199254740992.0);
    }
};

// Volume of the layer IsInside samples as outside: GEOMETRY_PRECISION times
// the exits that world +z rays meet over the solid's footprint.
double InsideTestLayer(siren::geometry::Geometry const & geometry) {
    using siren::math::Vector3D;
    auto aabb = geometry.GetWorldBoundingBox();
    double x0 = aabb.min_corner.GetX(), x1 = aabb.max_corner.GetX();
    double y0 = aabb.min_corner.GetY(), y1 = aabb.max_corner.GetY();
    double z0 = aabb.min_corner.GetZ() - 1e-6 * (aabb.max_corner.GetZ() - aabb.min_corner.GetZ()) - 1e-9;
    DeterministicUniform uniform;
    long long exits = 0;
    for (int i = 0; i < kRaysPerSide; ++i) {
        for (int j = 0; j < kRaysPerSide; ++j) {
            double ru = uniform();
            double rv = uniform();
            Vector3D start(x0 + (i + ru) / kRaysPerSide * (x1 - x0), y0 + (j + rv) / kRaysPerSide * (y1 - y0), z0);
            for (auto const & crossing : geometry.Intersections(start, Vector3D(0, 0, 1)))
                if (!crossing.entering && crossing.distance > 0.0) ++exits;
        }
    }
    return siren::geometry::Geometry::GEOMETRY_PRECISION * (x1 - x0) * (y1 - y0)
        * static_cast<double>(exits) / (static_cast<double>(kRaysPerSide) * kRaysPerSide);
}

} // namespace

double ExactGeometryVolume(siren::geometry::Geometry const & geometry) {
    if (auto const * cylinder =
            dynamic_cast<siren::geometry::Cylinder const *>(&geometry)) {
        double outer = cylinder->GetRadius();
        double inner = cylinder->GetInnerRadius();
        double delta_phi = AzimuthalExtent(cylinder->GetDeltaPhi());
        return 0.5 * delta_phi
            * (outer * outer - inner * inner) * cylinder->GetZ();
    }
    if (auto const * sphere =
            dynamic_cast<siren::geometry::Sphere const *>(&geometry)) {
        double outer = sphere->GetRadius();
        double inner = sphere->GetInnerRadius();
        double delta_phi = AzimuthalExtent(sphere->GetDeltaPhi());
        double theta_start = std::clamp(
            sphere->GetStartTheta(), 0.0, M_PI);
        double theta_end = std::clamp(
            theta_start + sphere->GetDeltaTheta(), theta_start, M_PI);
        double angular_integral = delta_phi
            * (std::cos(theta_start) - std::cos(theta_end));
        return (outer * outer * outer - inner * inner * inner)
            * angular_integral / 3.0;
    }
    if (auto const * box =
            dynamic_cast<siren::geometry::Box const *>(&geometry)) {
        return box->GetX() * box->GetY() * box->GetZ();
    }
    if (auto const * ellipsoid =
            dynamic_cast<siren::geometry::Ellipsoid const *>(&geometry)) {
        // Cross-section pi*a*b*(c - z)(c + z)/c^2 between the (clamped) z cuts,
        // integrated in the distance from the nearer pole so a thin cap does
        // not cancel: with u = c - z, the integral is
        // (u1 - u2)*(c*(u1 + u2) - (u1^2 + u1*u2 + u2^2)/3), and the same with
        // w = c + z for a slab in the lower half.
        double a = ellipsoid->GetAx(), b = ellipsoid->GetBy(), c = ellipsoid->GetCz();
        if (!(c > 0.0)) return 0.0;
        double z1 = std::clamp(ellipsoid->GetZcut1(), -c, c);
        double z2 = std::clamp(ellipsoid->GetZcut2(), z1, c);
        double near1 = z1 + z2 >= 0.0 ? c - z1 : c + z2;
        double near2 = z1 + z2 >= 0.0 ? c - z2 : c + z1;
        double bracket = c * (near1 + near2)
            - (near1 * near1 + near1 * near2 + near2 * near2) / 3.0;
        return M_PI * a * b / (c * c) * (z2 - z1) * bracket;
    }
    if (auto const * tube =
            dynamic_cast<siren::geometry::EllipticalTube const *>(&geometry)) {
        // Semi-axes dx, dy and half-length dz.
        return M_PI * tube->GetDx() * tube->GetDy() * 2.0 * tube->GetDz();
    }
    if (auto const * cone =
            dynamic_cast<siren::geometry::Cone const *>(&geometry)) {
        // Full height z; radii vary linearly between the two ends.
        auto frustum = [](double r1, double r2) { return r1 * r1 + r1 * r2 + r2 * r2; };
        return 0.5 * AzimuthalExtent(cone->GetDeltaPhi()) * cone->GetZ() / 3.0
            * (frustum(cone->GetRmax1(), cone->GetRmax2())
               - frustum(cone->GetRmin1(), cone->GetRmin2()));
    }
    if (auto const * torus =
            dynamic_cast<siren::geometry::Torus const *>(&geometry)) {
        // Pappus; a self-intersecting (spindle) torus has no such formula.
        double major = torus->GetMajorRadius();
        double outer = torus->GetMinorRadius(), inner = torus->GetInnerRadius();
        if (outer > major) return std::numeric_limits<double>::quiet_NaN();
        return AzimuthalExtent(torus->GetDeltaPhi()) * major
            * M_PI * (outer * outer - inner * inner);
    }
    if (auto const * trd =
            dynamic_cast<siren::geometry::Trd const *>(&geometry)) {
        // Half-widths vary linearly over the half-height dz.
        double x1 = trd->GetDx1(), x2 = trd->GetDx2();
        double y1 = trd->GetDy1(), y2 = trd->GetDy2();
        return 4.0 * trd->GetDz() / 3.0
            * (2.0 * x1 * y1 + 2.0 * x2 * y2 + x1 * y2 + x2 * y1);
    }
    if (auto const * para =
            dynamic_cast<siren::geometry::Para const *>(&geometry)) {
        // Shears do not change the volume of the half-length box.
        return 8.0 * para->GetDx() * para->GetDy() * para->GetDz();
    }
    if (auto const * polycone =
            dynamic_cast<siren::geometry::GenericPolycone const *>(&geometry)) {
        // Revolution of the R-Z polygon: delta_phi * integral of r dr dz, valid
        // only for a simple polygon (GenericPolycone does not require one).
        auto const & r = polycone->GetR();
        auto const & z = polycone->GetZ();
        if (!SimplePolygon(r, z)) return std::numeric_limits<double>::quiet_NaN();
        // The moment does not change with a shift in z; measuring z from the
        // profile's first vertex keeps a profile far along z from cancelling.
        double z0 = z.empty() ? 0.0 : z[0];
        double moment = 0.0, compensation = 0.0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            std::size_t j = (i + 1) % r.size();
            double term = (r[i] * (z[j] - z0) - r[j] * (z[i] - z0)) * (r[i] + r[j]);
            double next = moment + term;
            compensation += std::abs(moment) >= std::abs(term) ? (moment - next) + term : (term - next) + moment;
            moment = next;
        }
        return AzimuthalExtent(polycone->GetDeltaPhi()) * std::abs(moment + compensation) / 6.0;
    }
    if (auto const * polycone =
            dynamic_cast<siren::geometry::Polycone const *>(&geometry)) {
        // Frustum shells between consecutive planes (radii linear in z).
        auto const & z = polycone->GetZPlanes();
        auto const & rmin = polycone->GetRmin();
        auto const & rmax = polycone->GetRmax();
        double integral = 0.0;
        for (std::size_t k = 0; k + 1 < z.size(); ++k) {
            double outer = rmax[k] * rmax[k] + rmax[k] * rmax[k + 1] + rmax[k + 1] * rmax[k + 1];
            double inner = rmin[k] * rmin[k] + rmin[k] * rmin[k + 1] + rmin[k + 1] * rmin[k + 1];
            integral += (z[k + 1] - z[k]) * (outer - inner) / 3.0;
        }
        return 0.5 * AzimuthalExtent(polycone->GetDeltaPhi()) * integral;
    }
    if (auto const * mesh =
            dynamic_cast<siren::geometry::TriangularMesh const *>(&geometry)) {
        // Divergence theorem: signed tetrahedra over a closed, outward-oriented
        // surface, taken about the centre of its bounding box (about the origin
        // a mesh far away loses its volume to cancellation) and summed with
        // compensation. Only for one connected closed piece; the resolver also
        // requires the rays to resolve the mesh and agree, since a surface
        // that is closed edge by edge need not bound one solid.
        if (!mesh->ValidateClosed().empty() || MeshPieces(*mesh).size() != 1 || !MeshShellsValid(*mesh))
            return std::numeric_limits<double>::quiet_NaN();
        auto local = mesh->GetBoundingBox();
        double c[3] = {0.5 * (local.min_corner.GetX() + local.max_corner.GetX()),
                       0.5 * (local.min_corner.GetY() + local.max_corner.GetY()),
                       0.5 * (local.min_corner.GetZ() + local.max_corner.GetZ())};
        double sum = 0.0, compensation = 0.0;
        for (auto const & t : mesh->GetTriangles()) {
            // det(vk - c, vk+1 - vk, vk+2 - vk), the same value as
            // det(v0-c, v1-c, v2-c), from the vertex opposite the longest edge:
            // its two edges are the shortest, so a thin triangle's cross
            // product does not cancel between two nearly parallel long edges.
            int k = 0;
            double longest = -1.0;
            for (int q = 0; q < 3; ++q) {
                double s = 0.0;
                for (int m = 0; m < 3; ++m) s += (t[(q + 2) % 3][m] - t[(q + 1) % 3][m]) * (t[(q + 2) % 3][m] - t[(q + 1) % 3][m]);
                if (s > longest) { longest = s; k = q; }
            }
            double a[3], b[3], e[3];
            for (int m = 0; m < 3; ++m) {
                a[m] = t[k][m] - c[m];
                b[m] = t[(k + 1) % 3][m] - t[k][m];
                e[m] = t[(k + 2) % 3][m] - t[k][m];
            }
            double term = a[0] * (b[1] * e[2] - b[2] * e[1]) - a[1] * (b[0] * e[2] - b[2] * e[0])
                        + a[2] * (b[0] * e[1] - b[1] * e[0]);
            double next = sum + term;
            compensation += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
            sum = next;
        }
        double six = sum + compensation;
        return six > 0.0 ? six / 6.0 : std::numeric_limits<double>::quiet_NaN();
    }
    return std::numeric_limits<double>::quiet_NaN();
}

GeometryVolumeEstimate EstimateGeometryVolume(
    siren::geometry::Geometry const & geometry)
{
    using siren::math::Vector3D;
    GeometryVolumeEstimate estimate;
    auto aabb = geometry.GetWorldBoundingBox();
    // The ray grid is laid out in the solid's own frame, where its bounding
    // box is tight (a rotated placement inflates the world box, and a slender
    // solid then fills too little of the grid); the rays are mapped to world
    // coordinates. The volume does not change under the rigid placement.
    auto frame_box = geometry.GetBoundingBox();
    Vector3D extent = frame_box.max_corner - frame_box.min_corner;
    double box = extent.GetX() * extent.GetY() * extent.GetZ();
    if (!(box > 0.0) || !std::isfinite(box)) {
        estimate.reason = "its bounding box has no finite volume";
        return estimate;
    }
    double ca = std::cos(kSkewAngles[0]), sa = std::sin(kSkewAngles[0]);
    double cb = std::cos(kSkewAngles[1]), sb = std::sin(kSkewAngles[1]);
    double cg = std::cos(kSkewAngles[2]), sg = std::sin(kSkewAngles[2]);
    double R[3][3] = {
        {cg * cb, cg * sb * sa - sg * ca, cg * sb * ca + sg * sa},
        {sg * cb, sg * sb * sa + cg * ca, sg * sb * ca - cg * sa},
        {-sb, cb * sa, cb * ca}};
    std::vector<Vector3D> corners;
    for (int c = 0; c < 8; ++c)
        corners.emplace_back(c & 1 ? frame_box.max_corner.GetX() : frame_box.min_corner.GetX(),
                             c & 2 ? frame_box.max_corner.GetY() : frame_box.min_corner.GetY(),
                             c & 4 ? frame_box.max_corner.GetZ() : frame_box.min_corner.GetZ());
    // The parts, their world bounding boxes, and the volume each box holds
    // within the solid's box (its most the rays could miss).
    std::vector<SolidPart> parts;
    CollectParts(geometry, {}, parts);
    for (auto const * mesh : MeshesIn(geometry)) {
        if (!mesh->ValidateClosed().empty() || !MeshShellsValid(*mesh)) {
            estimate.reason = "its triangle mesh does not bound a solid (it is not closed, a shell wound"
                " inward is not a cavity inside an outward one, or outward shells are nested)";
            estimate.malformed = true;
            return estimate;
        }
    }
    for (auto & part : parts) {
        auto local = part.geometry->GetWorldBoundingBox();
        for (int m = 0; m < 3; ++m) { part.lo[m] = 1e300; part.hi[m] = -1e300; }
        for (int c = 0; c < 8; ++c) {
            Vector3D corner(c & 1 ? local.max_corner.GetX() : local.min_corner.GetX(),
                            c & 2 ? local.max_corner.GetY() : local.min_corner.GetY(),
                            c & 4 ? local.max_corner.GetZ() : local.min_corner.GetZ());
            for (auto it = part.ancestors.rbegin(); it != part.ancestors.rend(); ++it)
                corner = (*it)->LocalToGlobalPosition(corner);
            double w[3] = {corner.GetX(), corner.GetY(), corner.GetZ()};
            for (int m = 0; m < 3; ++m) { part.lo[m] = std::min(part.lo[m], w[m]); part.hi[m] = std::max(part.hi[m], w[m]); }
        }
        double solid_lo[3] = {aabb.min_corner.GetX(), aabb.min_corner.GetY(), aabb.min_corner.GetZ()};
        double solid_hi[3] = {aabb.max_corner.GetX(), aabb.max_corner.GetY(), aabb.max_corner.GetZ()};
        part.box = 1.0;
        bool inside = true;
        for (int m = 0; m < 3; ++m) {
            part.box *= std::max(0.0, std::min(part.hi[m], solid_hi[m]) - std::max(part.lo[m], solid_lo[m]));
            double slack = 1e-9 * (solid_hi[m] - solid_lo[m]);
            inside = inside && part.lo[m] >= solid_lo[m] - slack && part.hi[m] <= solid_hi[m] + slack;
        }
        // The rays cover the solid's box, so only a primitive inside it can be
        // estimated on its own and compared with its exact volume.
        double exact = ExactGeometryVolume(*part.geometry);
        if (inside && std::isfinite(exact) && exact > 0.0) part.exact = exact;
    }
    DeterministicUniform uniform;
    double const golden = M_PI * (3.0 - std::sqrt(5.0));
    std::vector<double> volumes;
    long long touched = 0, malformed = 0, chords = 0, thin_chords = 0;
    long long inside_tests = 0, inside_mismatches = 0;
    int sparse_directions = 0;
    for (int k = 0; k < kDirections; ++k) {
        double cz = (k + 0.5) / kDirections, rz = std::sqrt(1.0 - cz * cz);
        double base[3] = {rz * std::cos(golden * k), rz * std::sin(golden * k), cz};
        Vector3D frame_d(R[0][0] * base[0] + R[0][1] * base[1] + R[0][2] * base[2],
                         R[1][0] * base[0] + R[1][1] * base[1] + R[1][2] * base[2],
                         R[2][0] * base[0] + R[2][1] * base[1] + R[2][2] * base[2]);
        frame_d.normalize();
        // Orthonormal (u, v) across d giving the smallest rectangle around the
        // box's projection: that rectangle has a side along a projected box
        // edge, so u runs along the projection of x, y or z. A long thin solid
        // lying diagonally would otherwise fill only a sliver of the grid.
        Vector3D u, v;
        double lo[3], hi[3], area = std::numeric_limits<double>::infinity();
        for (Vector3D edge : {Vector3D(1, 0, 0), Vector3D(0, 1, 0), Vector3D(0, 0, 1)}) {
            Vector3D cu = edge - frame_d * siren::math::scalar_product(edge, frame_d);
            if (!(cu.magnitude() > 1e-6)) continue;
            cu.normalize();
            Vector3D cv = siren::math::vector_product(frame_d, cu);
            double clo[3] = {1e300, 1e300, 1e300}, chi[3] = {-1e300, -1e300, -1e300};
            for (auto const & corner : corners) {
                double c[3] = {siren::math::scalar_product(corner, cu),
                               siren::math::scalar_product(corner, cv),
                               siren::math::scalar_product(corner, frame_d)};
                for (int m = 0; m < 3; ++m) { clo[m] = std::min(clo[m], c[m]); chi[m] = std::max(chi[m], c[m]); }
            }
            double carea = (chi[0] - clo[0]) * (chi[1] - clo[1]);
            if (carea < area) {
                area = carea; u = cu; v = cv;
                std::copy(clo, clo + 3, lo); std::copy(chi, chi + 3, hi);
            }
        }
        double margin = 1e-6 * (hi[2] - lo[2]) + 1e-9;
        Vector3D d = geometry.LocalToGlobalDirection(frame_d);
        double sum = 0.0;
        int hits = 0, counted = 0;
        for (int i = 0; i < kRaysPerSide; ++i) {
            for (int j = 0; j < kRaysPerSide; ++j) {
                // Separate statements fix the draw order (argument order is not).
                double ru = uniform();
                double rv = uniform();
                Vector3D start = geometry.LocalToGlobalPosition(
                    u * (lo[0] + (i + ru) / kRaysPerSide * (hi[0] - lo[0]))
                    + v * (lo[1] + (j + rv) / kRaysPerSide * (hi[1] - lo[1])) + frame_d * (lo[2] - margin));
                auto crossings = geometry.Intersections(start, d);
                for (auto & part : parts) {
                    bool compared = std::isfinite(part.exact);
                    if ((!compared && part.hits >= kPartResolvedHits) || !RayMeetsBox(start, d, part.lo, part.hi)) continue;
                    if (part.geometry == &geometry) {
                        part.sum += ChordLength(crossings);
                        ++part.hits;
                        continue;
                    }
                    Vector3D p = start, q = d;
                    for (auto const * ancestor : part.ancestors) {
                        p = ancestor->GlobalToLocalPosition(p);
                        q = ancestor->GlobalToLocalDirection(q);
                    }
                    auto part_crossings = part.geometry->Intersections(p, q);
                    if (part_crossings.empty()) continue;
                    ++part.hits;
                    if (compared) part.sum += ChordLength(std::move(part_crossings));
                }
                std::sort(crossings.begin(), crossings.end(),
                    [](auto const & a, auto const & b) { return a.distance < b.distance; });
                // The ray starts outside the bounding box, so a well-formed
                // solid's crossings alternate entering, exiting, ... and end
                // outside. Anything else (an open, doubled, nested or
                // self-intersecting surface) is a malformed ray; a non-finite
                // distance makes the whole estimate inconclusive.
                bool alternating = crossings.size() % 2 == 0;
                for (std::size_t c = 0; c < crossings.size(); ++c) {
                    double distance = crossings[c].distance;
                    if (!std::isfinite(distance)) {
                        estimate.reason = "its surface crossings include a non-finite distance";
                        estimate.malformed = true;
                        return estimate;
                    }
                    if (!(distance > 0.0) || crossings[c].entering != (c % 2 == 0)) alternating = false;
                }
                if (!crossings.empty()) ++touched;
                if (!alternating) { ++malformed; continue; }
                ++counted;
                double chord = 0.0;
                for (std::size_t c = 0; c + 1 < crossings.size(); c += 2) {
                    double length = crossings[c + 1].distance - crossings[c].distance;
                    chord += length;
                    ++chords;
                    if (length < kChordResolution) ++thin_chords;
                }
                sum += chord;
                if (!(chord > 0.0)) continue;
                ++hits;
                if (i % kInsideCheckStride != 0 || j % kInsideCheckStride != 0) continue;
                // The channel samples points with IsInside, so it must agree
                // with the crossings: inside between an entry and its exit,
                // outside between an exit and the next entry.
                for (std::size_t c = 0; c + 1 < crossings.size(); ++c) {
                    double a = crossings[c].distance, b = crossings[c + 1].distance;
                    if (b - a < kInsideCheckLength) continue;
                    ++inside_tests;
                    if (geometry.IsInside(start + d * (0.5 * (a + b))) != (c % 2 == 0)) ++inside_mismatches;
                }
            }
        }
        estimate.points += static_cast<long long>(kRaysPerSide) * kRaysPerSide;
        if (hits < kResolvedHits) ++sparse_directions;
        volumes.push_back(counted > 0 ? area * sum / counted : 0.0);
        double rays = static_cast<double>(kRaysPerSide) * kRaysPerSide;
        for (auto & part : parts) {
            if (!std::isfinite(part.exact)) continue;
            part.volumes.push_back(area * part.sum / rays);
            part.sum = 0.0;
        }
    }
    auto count = [](long long n, long long of) {
        return std::to_string(n) + " of " + std::to_string(of);
    };
    if (malformed > std::max(2.0, kMalformedFraction * touched)) {
        estimate.reason = "its surface crossings do not alternate between entering and exiting on "
            + count(malformed, touched) + " rays that meet it (an open, doubled, nested or"
            " self-intersecting surface)";
        estimate.malformed = true;
        return estimate;
    }
    if (inside_mismatches > std::max(2.0, kMalformedFraction * inside_tests)) {
        estimate.reason = "its inside test disagrees with its surface crossings at "
            + count(inside_mismatches, inside_tests) + " checked points";
        estimate.malformed = true;
        return estimate;
    }
    if (thin_chords > kSubResolutionFraction * chords) {
        estimate.reason = count(thin_chords, chords) + " chords through it are shorter than "
            + Format(kChordResolution) + " m, where crossings are merged and material can be lost";
        return estimate;
    }
    if (sparse_directions > kSparseDirections) {
        estimate.reason = "fewer than " + std::to_string(kResolvedHits) + " rays meet it along "
            + std::to_string(sparse_directions) + " of " + std::to_string(kDirections)
            + " directions (a needle, a speck or an empty solid)";
        return estimate;
    }
    double mean = 0.0;
    for (double value : volumes) mean += value;
    mean /= volumes.size();
    double scatter = 0.0;
    for (double value : volumes) scatter += (value - mean) * (value - mean);
    double error = std::sqrt(scatter / (volumes.size() - 1.0) / volumes.size());
    estimate.volume = mean;
    estimate.standard_error = error;
    if (!std::isfinite(mean) || !std::isfinite(error) || !(mean > 0.0)) {
        estimate.reason = "its estimate is not a positive finite number";
        return estimate;
    }
    double unseen = 0.0, largest_missed = 0.0;
    std::string largest;
    for (auto const & part : parts) {
        double missed = 0.0;
        std::string what;
        if (std::isfinite(part.exact)) {
            double own = 0.0, own_scatter = 0.0;
            for (double value : part.volumes) own += value;
            own /= part.volumes.size();
            for (double value : part.volumes) own_scatter += (value - own) * (value - own);
            double own_error = std::sqrt(own_scatter / (part.volumes.size() - 1.0) / part.volumes.size());
            missed = std::max(0.0, std::abs(part.exact - own) - 5.0 * own_error);
            what = part.name + " of exact volume " + Format(part.exact) + " is estimated by the rays at "
                + Format(own) + " +- " + Format(own_error);
        } else if (part.hits < kPartResolvedHits) {
            missed = part.box;
            what = part.name + " is met by only " + std::to_string(part.hits)
                + " rays and its bounding box holds " + Format(part.box);
        }
        unseen += missed;
        if (missed > largest_missed) { largest_missed = missed; largest = what; }
    }
    if (unseen > kUnseenFraction * mean) {
        estimate.reason = "part of it may be missed: " + largest + "; the parts the rays do not"
            " resolve could hold up to " + Format(unseen) + " (the estimate is " + Format(mean) + ")";
        return estimate;
    }
    if (mean > box * (1.0 + AABBRounding(geometry)) + 5.0 * error) {
        estimate.reason = "its estimate " + Format(mean) + " exceeds its bounding-box volume "
            + Format(box);
        return estimate;
    }
    if (error > kMaximumRelativeError * mean) {
        estimate.reason = "its estimate " + Format(mean) + " is too uncertain to check (standard"
            " error " + Format(error) + " from the scatter between ray directions)";
        return estimate;
    }
    double allowed = std::max(kEstimateRelativeTolerance * mean, kEstimateStandardErrors * error);
    estimate.lower = mean - allowed;
    estimate.upper = mean + allowed;
    estimate.resolved = true;
    return estimate;
}

double ResolveDetectorDirectedVolume(
    siren::geometry::Geometry const & geometry,
    bool volume_mode,
    double supplied_volume)
{
    double aabb_volume = AABBVolume(geometry);
    bool has_supplied_volume = supplied_volume > 0.0;
    double exact_volume = ExactGeometryVolume(geometry);
    double target_volume = has_supplied_volume ? supplied_volume : exact_volume;

    if (!volume_mode) return target_volume;
    if (!(aabb_volume > 0.0) || !std::isfinite(aabb_volume)) {
        throw std::runtime_error("Target bounding box has zero volume");
    }
    // An exact volume cannot exceed the bounding box.
    if (std::isfinite(exact_volume) && exact_volume > aabb_volume * (1.0 + AABBRounding(geometry))) {
        throw std::runtime_error(
            "Target's exact volume " + Format(exact_volume) + " exceeds its bounding-box volume "
            + Format(aabb_volume) + "; its shape parameters are inconsistent");
    }
    // A mesh's exact volume (divergence theorem) is the volume of the solid
    // the channel samples only if its surface bounds one solid, which a closed
    // surface need not do (a nested or self-intersecting piece), and rays can
    // only check where they go: it is used only when the rays resolve the mesh
    // and their estimate agrees with it.
    if (std::isfinite(exact_volume) && dynamic_cast<siren::geometry::TriangularMesh const *>(&geometry)) {
        GeometryVolumeEstimate check = EstimateGeometryVolume(geometry);
        std::string problem = check.reason;
        if (check.resolved && (exact_volume < check.lower || exact_volume > check.upper)) {
            problem = "its exact volume " + Format(exact_volume) + " disagrees with the rays' estimate "
                + Format(check.volume) + " (accepted range " + Format(check.lower) + " to "
                + Format(check.upper) + "), so its surface does not bound the solid it samples";
        }
        if (!check.resolved || !problem.empty()) {
            throw std::runtime_error(
                (has_supplied_volume ? "Supplied target volume " + Format(supplied_volume) + " cannot be checked: "
                                     : std::string("Target mesh volume cannot be checked: ")) + problem);
        }
    }
    if (!(target_volume > 0.0) || !std::isfinite(target_volume)) {
        throw std::runtime_error(
            "Volume mode requires an exact caller-supplied volume for this geometry");
    }
    // The volume sets the proposal density while points come from the solid
    // itself, so a wrong value biases every weight. Where the solid has an
    // analytic volume, a supplied one must agree and the analytic one is used.
    // Otherwise it must agree with an independent chord-integration estimate
    // of the solid, and a solid the estimate cannot resolve refuses it.
    if (has_supplied_volume && exact_volume > 0.0 && std::isfinite(exact_volume)) {
        if (std::abs(supplied_volume - exact_volume) > 1e-4 * exact_volume) {
            throw std::runtime_error(
                "Supplied target volume " + Format(supplied_volume)
                + " disagrees with the analytic volume " + Format(exact_volume));
        }
        target_volume = exact_volume;
    } else if (has_supplied_volume) {
        // A solid cannot exceed its bounding box (up to the box's own rounding).
        if (supplied_volume > aabb_volume * (1.0 + AABBRounding(geometry))) {
            throw std::runtime_error(
                "Supplied target volume " + Format(supplied_volume)
                + " exceeds the target's bounding-box volume " + Format(aabb_volume));
        }
        GeometryVolumeEstimate estimate = EstimateGeometryVolume(geometry);
        if (!estimate.resolved || !std::isfinite(estimate.lower) || !std::isfinite(estimate.upper)) {
            // A normalization that cannot be checked is not accepted on trust.
            throw std::runtime_error(
                "Supplied target volume " + Format(supplied_volume) + " cannot be checked: "
                + (estimate.reason.empty() ? std::string("the estimate is not finite") : estimate.reason));
        }
        // Written so that a NaN bound refuses rather than accepts.
        if (!(supplied_volume >= estimate.lower && supplied_volume <= estimate.upper)) {
            throw std::runtime_error(
                "Supplied target volume " + Format(supplied_volume)
                + " disagrees with the chord-integration estimate " + Format(estimate.volume)
                + " of the solid (accepted range " + Format(estimate.lower) + " to "
                + Format(estimate.upper) + ")");
        }
    }
    if (!has_supplied_volume &&
        target_volume / aabb_volume <= kMinimumViableFillFraction) {
        throw std::runtime_error(
            "Target volume is too small relative to its bounding box for sampling to be viable");
    }
    double layer = InsideTestLayer(geometry);
    if (layer > kInsideLayerFraction * target_volume) {
        throw std::runtime_error(
            "Target volume " + Format(target_volume) + " cannot be sampled consistently: the inside test"
            " the channel samples points with treats the last " + Format(siren::geometry::Geometry::GEOMETRY_PRECISION)
            + " m before each exit along +z as outside, about " + Format(layer) + " of this solid (more than "
            + Format(100.0 * kInsideLayerFraction) + "%); it is too thin");
    }
    return target_volume;
}

void ValidateArchivedDetectorDirectedVolume(
    siren::geometry::Geometry const * geometry,
    bool volume_mode,
    double archived_volume)
{
    if (!volume_mode) return;
    try {
        if (!geometry) throw std::runtime_error("no target geometry");
        if (!(archived_volume > 0.0) || !std::isfinite(archived_volume))
            throw std::runtime_error("the volume is not positive and finite");
        ResolveDetectorDirectedVolume(*geometry, true, archived_volume);
    } catch (std::runtime_error const & error) {
        throw std::runtime_error(
            std::string("Archived directed-channel target volume rejected: ")
            + error.what() + "; rebuild the channel with a validated volume");
    }
}

} // namespace injection
} // namespace siren
