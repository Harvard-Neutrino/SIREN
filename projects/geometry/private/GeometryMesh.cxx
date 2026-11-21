#include "SIREN/geometry/GeometryMesh.h"

#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <sstream>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

// =========================================================================
// Woop's watertight ray-triangle intersection
//
// Ensures that a ray hitting a shared edge is assigned to exactly one
// adjacent triangle via consistent tie-breaking. No epsilon tolerance needed.
//
// Sven Woop, Carsten Benthin, Ingo Wald, "Watertight Ray/Triangle
// Intersection," Journal of Computer Graphics Techniques 2(1), 2013.
// https://jcgt.org/published/0002/01/05/
// =========================================================================

namespace {

struct WoopRay {
    double org[3];
    double dir[3];
    int kx, ky, kz;
    double Sx, Sy, Sz;
};

WoopRay PrepareRay(double ox, double oy, double oz, double dx, double dy, double dz) {
    WoopRay r;
    r.org[0] = ox; r.org[1] = oy; r.org[2] = oz;
    r.dir[0] = dx; r.dir[1] = dy; r.dir[2] = dz;

    double abs_dx = std::fabs(dx);
    double abs_dy = std::fabs(dy);
    double abs_dz = std::fabs(dz);

    r.kz = 0;
    if(abs_dy > abs_dx) { r.kz = 1; }
    if(abs_dz > (r.kz == 0 ? abs_dx : abs_dy)) { r.kz = 2; }
    r.kx = (r.kz + 1) % 3;
    r.ky = (r.kz + 2) % 3;

    if(r.dir[r.kz] < 0.0) {
        int tmp = r.kx; r.kx = r.ky; r.ky = tmp;
    }

    r.Sx = r.dir[r.kx] / r.dir[r.kz];
    r.Sy = r.dir[r.ky] / r.dir[r.kz];
    r.Sz = 1.0 / r.dir[r.kz];

    return r;
}

struct TriHit {
    double t;
    bool hit;
    bool front_face;
};

// Prevent inlining so the compiler cannot apply FMA contraction or
// vectorization across the call boundary. The edge function computation
// requires independent rounding of each multiply for watertight results.
#ifdef _MSC_VER
__declspec(noinline)
#else
__attribute__((noinline))
#endif
TriHit WoopIntersect(WoopRay const & ray, double const * v0, double const * v1, double const * v2) {
    TriHit result;
    result.hit = false;
    result.t = 0;
    result.front_face = false;

    double A[3] = {v0[0] - ray.org[0], v0[1] - ray.org[1], v0[2] - ray.org[2]};
    double B[3] = {v1[0] - ray.org[0], v1[1] - ray.org[1], v1[2] - ray.org[2]};
    double C[3] = {v2[0] - ray.org[0], v2[1] - ray.org[1], v2[2] - ray.org[2]};

    double Ax = A[ray.kx] - ray.Sx * A[ray.kz];
    double Ay = A[ray.ky] - ray.Sy * A[ray.kz];
    double Bx = B[ray.kx] - ray.Sx * B[ray.kz];
    double By = B[ray.ky] - ray.Sy * B[ray.kz];
    double Cx = C[ray.kx] - ray.Sx * C[ray.kz];
    double Cy = C[ray.ky] - ray.Sy * C[ray.kz];

    double U = Cx * By - Cy * Bx;
    double V = Ax * Cy - Ay * Cx;
    double W = Bx * Ay - By * Ax;

    if((U < 0.0 || V < 0.0 || W < 0.0) && (U > 0.0 || V > 0.0 || W > 0.0)) {
        return result;
    }

    double det = U + V + W;
    if(det == 0.0) return result;

    double Az = ray.Sz * A[ray.kz];
    double Bz = ray.Sz * B[ray.kz];
    double Cz = ray.Sz * C[ray.kz];
    double T = U * Az + V * Bz + W * Cz;

    result.hit = true;
    result.t = T / det;
    result.front_face = (det > 0.0);
    return result;
}

// Ray-AABB slab test (full-line, handles direction=0 correctly)
bool RayAABBIntersect(float const * bounds, double const * org, double const * inv_dir) {
    double tmin, tmax;

    // For each axis: if direction is 0, check containment directly
    // Otherwise use standard slab method
    if(inv_dir[0] == 0.0) {
        // Parallel to x slabs: miss if origin outside
        if(org[0] < bounds[0] || org[0] > bounds[3]) return false;
        tmin = -1e308;
        tmax = 1e308;
    } else {
        double t1 = (bounds[0] - org[0]) * inv_dir[0];
        double t2 = (bounds[3] - org[0]) * inv_dir[0];
        tmin = std::fmin(t1, t2);
        tmax = std::fmax(t1, t2);
    }

    if(inv_dir[1] == 0.0) {
        if(org[1] < bounds[1] || org[1] > bounds[4]) return false;
    } else {
        double t1 = (bounds[1] - org[1]) * inv_dir[1];
        double t2 = (bounds[4] - org[1]) * inv_dir[1];
        tmin = std::fmax(tmin, std::fmin(t1, t2));
        tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    if(inv_dir[2] == 0.0) {
        if(org[2] < bounds[2] || org[2] > bounds[5]) return false;
    } else {
        double t1 = (bounds[2] - org[2]) * inv_dir[2];
        double t2 = (bounds[5] - org[2]) * inv_dir[2];
        tmin = std::fmax(tmin, std::fmin(t1, t2));
        tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    return tmax >= tmin;
}

} // anonymous namespace

// =========================================================================
// Constructors
// =========================================================================

TriangularMesh::TriangularMesh()
    : Geometry("TriangularMesh")
{ RecomputeWorldAABB(); }

TriangularMesh::TriangularMesh(std::vector<std::array<math::Vector3D, 3>> const & triangles)
    : Geometry("TriangularMesh")
{
    triangles_.reserve(triangles.size());
    for(auto const & tri : triangles) {
        Triangle t = {{
            {{tri[0].GetX(), tri[0].GetY(), tri[0].GetZ()}},
            {{tri[1].GetX(), tri[1].GetY(), tri[1].GetZ()}},
            {{tri[2].GetX(), tri[2].GetY(), tri[2].GetZ()}}
        }};
        triangles_.push_back(t);
    }
    BuildBVH();
    RecomputeWorldAABB();
}

TriangularMesh::TriangularMesh(Placement const & placement, std::vector<std::array<math::Vector3D, 3>> const & triangles)
    : Geometry("TriangularMesh", placement)
{
    triangles_.reserve(triangles.size());
    for(auto const & tri : triangles) {
        Triangle t = {{
            {{tri[0].GetX(), tri[0].GetY(), tri[0].GetZ()}},
            {{tri[1].GetX(), tri[1].GetY(), tri[1].GetZ()}},
            {{tri[2].GetX(), tri[2].GetY(), tri[2].GetZ()}}
        }};
        triangles_.push_back(t);
    }
    BuildBVH();
    RecomputeWorldAABB();
}

TriangularMesh::TriangularMesh(Placement const & placement)
    : Geometry("TriangularMesh", placement)
{ RecomputeWorldAABB(); }

TriangularMesh::TriangularMesh(TriangularMesh const & other)
    : Geometry(other)
    , triangles_(other.triangles_)
    , bvh_nodes_(other.bvh_nodes_)
    , tri_indices_(other.tri_indices_)
{ RecomputeWorldAABB(); }

// =========================================================================
// Geometry interface
// =========================================================================

void TriangularMesh::swap(Geometry& geometry) {
    TriangularMesh* tmesh = dynamic_cast<TriangularMesh*>(&geometry);
    if(!tmesh) throw std::runtime_error("Cannot swap TriangularMesh with non-TriangularMesh");
    Geometry::swap(*tmesh);
    std::swap(triangles_, tmesh->triangles_);
    std::swap(bvh_nodes_, tmesh->bvh_nodes_);
    std::swap(tri_indices_, tmesh->tri_indices_);
}

TriangularMesh& TriangularMesh::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const TriangularMesh* tmesh = dynamic_cast<const TriangularMesh*>(&geometry);
        if(!tmesh) throw std::runtime_error("Cannot assign non-TriangularMesh to TriangularMesh");
        Geometry::operator=(geometry);
        triangles_ = tmesh->triangles_;
        bvh_nodes_ = tmesh->bvh_nodes_;
        tri_indices_ = tmesh->tri_indices_;
    }
    return *this;
}

bool TriangularMesh::equal(const Geometry& geometry) const {
    const TriangularMesh* other = dynamic_cast<const TriangularMesh*>(&geometry);
    if(!other) return false;
    return triangles_ == other->triangles_;
}

bool TriangularMesh::less(const Geometry& geometry) const {
    const TriangularMesh* other = dynamic_cast<const TriangularMesh*>(&geometry);
    if(!other) return false;
    return triangles_ < other->triangles_;
}

void TriangularMesh::print(std::ostream& os) const {
    os << "TriangularMesh(" << triangles_.size() << " triangles)";
}

// =========================================================================
// BVH construction (top-down, median split, flat array)
// =========================================================================

void TriangularMesh::BuildBVH() {
    bvh_nodes_.clear();
    tri_indices_.clear();

    size_t n = triangles_.size();
    if(n == 0) return;

    tri_indices_.resize(n);
    std::iota(tri_indices_.begin(), tri_indices_.end(), 0);

    // Precompute centroids and bounds per triangle
    struct TriInfo {
        float centroid[3];
        float bounds[6]; // min_x, min_y, min_z, max_x, max_y, max_z
    };
    std::vector<TriInfo> info(n);
    for(size_t i = 0; i < n; ++i) {
        auto const & tri = triangles_[i];
        for(int a = 0; a < 3; ++a) {
            float mn = (float)std::fmin(tri[0][a], std::fmin(tri[1][a], tri[2][a]));
            float mx = (float)std::fmax(tri[0][a], std::fmax(tri[1][a], tri[2][a]));
            info[i].bounds[a] = mn;
            info[i].bounds[a + 3] = mx;
            info[i].centroid[a] = (float)((tri[0][a] + tri[1][a] + tri[2][a]) / 3.0);
        }
    }

    // Reserve space for BVH nodes (at most 2n-1 for a full binary tree)
    bvh_nodes_.reserve(2 * n);

    // Recursive build using a stack to avoid deep recursion
    struct BuildTask {
        uint32_t node_idx;
        uint32_t begin;
        uint32_t end;
    };
    std::vector<BuildTask> stack;

    // Create root node
    bvh_nodes_.push_back(BVHNode{});
    stack.push_back({0, 0, (uint32_t)n});

    static const uint32_t MAX_LEAF_SIZE = 8;

    while(!stack.empty()) {
        BuildTask task = stack.back();
        stack.pop_back();

        uint32_t begin = task.begin;
        uint32_t end = task.end;
        uint32_t count = end - begin;
        BVHNode& node = bvh_nodes_[task.node_idx];

        // Compute bounds for this node
        float b[6] = {1e30f, 1e30f, 1e30f, -1e30f, -1e30f, -1e30f};
        for(uint32_t i = begin; i < end; ++i) {
            uint32_t ti = tri_indices_[i];
            for(int a = 0; a < 3; ++a) {
                b[a] = std::fmin(b[a], info[ti].bounds[a]);
                b[a + 3] = std::fmax(b[a + 3], info[ti].bounds[a + 3]);
            }
        }
        for(int i = 0; i < 6; ++i) node.bounds[i] = b[i];

        if(count <= MAX_LEAF_SIZE) {
            node.offset = begin;
            node.count = (uint16_t)count;
            node.axis = 0;
            continue;
        }

        // Find longest axis of centroid bounds
        float cmin[3] = {1e30f, 1e30f, 1e30f};
        float cmax[3] = {-1e30f, -1e30f, -1e30f};
        for(uint32_t i = begin; i < end; ++i) {
            uint32_t ti = tri_indices_[i];
            for(int a = 0; a < 3; ++a) {
                cmin[a] = std::fmin(cmin[a], info[ti].centroid[a]);
                cmax[a] = std::fmax(cmax[a], info[ti].centroid[a]);
            }
        }

        int axis = 0;
        float extent = cmax[0] - cmin[0];
        if(cmax[1] - cmin[1] > extent) { axis = 1; extent = cmax[1] - cmin[1]; }
        if(cmax[2] - cmin[2] > extent) { axis = 2; }

        // Degenerate case: all centroids coincident
        if(cmax[axis] - cmin[axis] < 1e-10f) {
            node.offset = begin;
            node.count = (uint16_t)std::min(count, (uint32_t)65535);
            node.axis = 0;
            continue;
        }

        // Partition at median centroid
        uint32_t mid = begin + count / 2;
        std::nth_element(
            tri_indices_.begin() + begin,
            tri_indices_.begin() + mid,
            tri_indices_.begin() + end,
            [&info, axis](uint32_t a, uint32_t b) {
                return info[a].centroid[axis] < info[b].centroid[axis];
            });

        node.axis = (uint16_t)axis;
        node.count = 0; // internal node

        // Allocate child nodes
        uint32_t left_idx = (uint32_t)bvh_nodes_.size();
        bvh_nodes_.push_back(BVHNode{});
        bvh_nodes_.push_back(BVHNode{});
        node.offset = left_idx + 1; // right child index

        // Push right first (processed second = depth-first left)
        stack.push_back({left_idx + 1, mid, end});
        stack.push_back({left_idx, begin, mid});
    }
}

// =========================================================================
// Ray intersection (BVH traversal + Woop's watertight test)
// =========================================================================

std::vector<Geometry::Intersection> TriangularMesh::ComputeIntersections(
    math::Vector3D const & position, math::Vector3D const & direction) const {

    std::vector<Intersection> result;
    if(triangles_.empty() || bvh_nodes_.empty()) return result;

    double org[3] = {position.GetX(), position.GetY(), position.GetZ()};
    double dir[3] = {direction.GetX(), direction.GetY(), direction.GetZ()};
    double inv_dir[3] = {
        dir[0] == 0.0 ? 0.0 : 1.0 / dir[0],
        dir[1] == 0.0 ? 0.0 : 1.0 / dir[1],
        dir[2] == 0.0 ? 0.0 : 1.0 / dir[2]
    };

    WoopRay wray = PrepareRay(org[0], org[1], org[2], dir[0], dir[1], dir[2]);

    // Iterative BVH traversal
    static constexpr int BVH_STACK_CAP = 64;
    uint32_t node_stack[BVH_STACK_CAP];
    int stack_ptr = 0;
    node_stack[0] = 0;
    std::vector<uint32_t> heap_stack;

    while(stack_ptr >= 0 || !heap_stack.empty()) {
        uint32_t idx;
        if(!heap_stack.empty()) {
            idx = heap_stack.back();
            heap_stack.pop_back();
        } else {
            idx = node_stack[stack_ptr--];
        }
        BVHNode const & node = bvh_nodes_[idx];

        if(!RayAABBIntersect(node.bounds, org, inv_dir)) {
            continue;
        }

        if(node.count > 0) {
            // Leaf node: test triangles
            for(uint32_t i = 0; i < node.count; ++i) {
                uint32_t ti = tri_indices_[node.offset + i];
                auto const & tri = triangles_[ti];
                TriHit hit = WoopIntersect(wray, tri[0].data(), tri[1].data(), tri[2].data());
                if(hit.hit) {
                    Intersection isect;
                    isect.distance = hit.t;
                    isect.hierarchy = 0;
                    isect.entering = hit.front_face;
                    isect.matID = 0;
                    isect.position = position + direction * hit.t;
                    result.push_back(isect);
                }
            }
        } else {
            // Internal node: push children, spilling to heap if stack is full
            if(stack_ptr + 2 < BVH_STACK_CAP - 1) {
                node_stack[++stack_ptr] = node.offset;     // right child
                node_stack[++stack_ptr] = node.offset - 1; // left child
            } else {
                heap_stack.push_back(node.offset);
                heap_stack.push_back(node.offset - 1);
            }
        }
    }

    // Sort by distance
    std::sort(result.begin(), result.end(),
        [](Intersection const & a, Intersection const & b) {
            return a.distance < b.distance;
        });

    // Merge coincident hits from adjacent triangles sharing an edge/vertex.
    // With double-precision Woop, shared edges produce edge functions of exactly
    // zero in both adjacent triangles (same vertex data, same arithmetic), so
    // both report a hit at the same distance. Collapse consecutive duplicates.
    if(result.size() > 1) {
        std::vector<Intersection> merged;
        merged.reserve(result.size());
        merged.push_back(result[0]);
        for(size_t i = 1; i < result.size(); ++i) {
            double prev_t = merged.back().distance;
            double curr_t = result[i].distance;
            double scale = std::fmax(1.0, std::fmax(std::fabs(prev_t), std::fabs(curr_t)));
            if(std::fabs(curr_t - prev_t) < scale * 1e-12) {
                continue;
            }
            merged.push_back(result[i]);
        }
        result = std::move(merged);
    }

    return result;
}

// =========================================================================
// Bounding box
// =========================================================================

AABB TriangularMesh::GetBoundingBox() const {
    AABB box;
    for(auto const & tri : triangles_) {
        for(int v = 0; v < 3; ++v) {
            box.ExpandToInclude(math::Vector3D(tri[v][0], tri[v][1], tri[v][2]));
        }
    }
    return box;
}

// =========================================================================
// Closed mesh validation
//
// A closed 2-manifold mesh has the property that every edge is shared by
// exactly 2 triangles. Boundary edges (shared by 1) indicate holes.
// Non-manifold edges (shared by 3+) indicate self-intersection.
// =========================================================================

std::string TriangularMesh::ValidateClosed() const {
    if(triangles_.empty()) return "mesh has no triangles";

    // Quantize vertex positions to detect shared vertices.
    // Two vertices are "same" if they are within 1e-10 in all coordinates.
    // Use integer grid for exact comparison.
    static const double GRID = 1e10;
    struct VertexKey {
        int64_t x, y, z;
        bool operator<(VertexKey const & o) const {
            if(x != o.x) return x < o.x;
            if(y != o.y) return y < o.y;
            return z < o.z;
        }
    };

    std::map<VertexKey, uint32_t> vertex_map;
    auto getVertexId = [&](double const * v) -> uint32_t {
        VertexKey key = {
            (int64_t)std::round(v[0] * GRID),
            (int64_t)std::round(v[1] * GRID),
            (int64_t)std::round(v[2] * GRID)
        };
        auto it = vertex_map.find(key);
        if(it != vertex_map.end()) return it->second;
        uint32_t id = (uint32_t)vertex_map.size();
        vertex_map[key] = id;
        return id;
    };

    // Build edge-face adjacency
    struct EdgeKey {
        uint32_t v0, v1;
        bool operator<(EdgeKey const & o) const {
            if(v0 != o.v0) return v0 < o.v0;
            return v1 < o.v1;
        }
    };

    std::map<EdgeKey, int> edge_count;
    for(size_t i = 0; i < triangles_.size(); ++i) {
        uint32_t ids[3];
        for(int v = 0; v < 3; ++v) {
            ids[v] = getVertexId(triangles_[i][v].data());
        }

        // Check for degenerate triangle (two or more coincident vertices)
        if(ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0]) {
            continue; // skip degenerate triangles
        }

        for(int e = 0; e < 3; ++e) {
            uint32_t a = ids[e];
            uint32_t b = ids[(e + 1) % 3];
            EdgeKey ek = {std::min(a, b), std::max(a, b)};
            edge_count[ek]++;
        }
    }

    int boundary_edges = 0;
    int nonmanifold_edges = 0;
    for(auto const & pair : edge_count) {
        if(pair.second == 1) ++boundary_edges;
        else if(pair.second > 2) ++nonmanifold_edges;
    }

    if(boundary_edges > 0 || nonmanifold_edges > 0) {
        std::ostringstream oss;
        if(boundary_edges > 0)
            oss << boundary_edges << " boundary edge(s) (mesh has holes)";
        if(nonmanifold_edges > 0) {
            if(boundary_edges > 0) oss << "; ";
            oss << nonmanifold_edges << " non-manifold edge(s)";
        }
        return oss.str();
    }

    return "";
}

} // namespace geometry
} // namespace siren
