#include "SIREN/geometry/BooleanGeometry.h"

#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <functional>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

// Helper: evaluate CSG inside predicate from child states
static bool IsInsideResult(BooleanOperation op, bool in_left, bool in_right) {
    switch(op) {
        case BooleanOperation::UNION:
            return in_left || in_right;
        case BooleanOperation::INTERSECTION:
            return in_left && in_right;
        case BooleanOperation::SUBTRACTION:
            return in_left && !in_right;
    }
    return false;
}

BooleanGeometry::BooleanGeometry()
    : Geometry((std::string)("BooleanGeometry"))
    , op_(BooleanOperation::UNION)
    , left_(nullptr)
      , right_(nullptr) {
}

BooleanGeometry::BooleanGeometry(BooleanOperation op, std::shared_ptr<const Geometry> left, std::shared_ptr<const Geometry> right)
    : Geometry((std::string)("BooleanGeometry"))
    , op_(op)
    , left_(left)
      , right_(right) {
}

BooleanGeometry::BooleanGeometry(Placement const & placement, BooleanOperation op, std::shared_ptr<const Geometry> left, std::shared_ptr<const Geometry> right)
    : Geometry((std::string)("BooleanGeometry"), placement)
    , op_(op)
    , left_(left)
      , right_(right) {
}

BooleanGeometry::BooleanGeometry(const BooleanGeometry& boolean_geometry)
    : Geometry(boolean_geometry)
    , op_(boolean_geometry.op_)
    , left_(boolean_geometry.left_)
      , right_(boolean_geometry.right_) {
}

// ------------------------------------------------------------------------- //
void BooleanGeometry::swap(Geometry& geometry) {
    BooleanGeometry* bg = dynamic_cast<BooleanGeometry*>(&geometry);
    if(!bg) {
        return;
    }

    Geometry::swap(*bg);

    std::swap(op_, bg->op_);
    std::swap(left_, bg->left_);
    std::swap(right_, bg->right_);
}

// ------------------------------------------------------------------------- //
BooleanGeometry& BooleanGeometry::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const BooleanGeometry* bg = dynamic_cast<const BooleanGeometry*>(&geometry);
        if(!bg) {
            return *this;
        }

        BooleanGeometry tmp(*bg);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool BooleanGeometry::equal(const Geometry& geometry) const
{
    const BooleanGeometry* bg = dynamic_cast<const BooleanGeometry*>(&geometry);

    if(!bg)
        return false;
    if(op_ != bg->op_)
        return false;
    // Compare children by value if both non-null
    if((left_ == nullptr) != (bg->left_ == nullptr))
        return false;
    if((right_ == nullptr) != (bg->right_ == nullptr))
        return false;
    if(left_ && bg->left_ && *left_ != *(bg->left_))
        return false;
    if(right_ && bg->right_ && *right_ != *(bg->right_))
        return false;
    return true;
}

// ------------------------------------------------------------------------- //
bool BooleanGeometry::less(const Geometry& geometry) const
{
    const BooleanGeometry* bg = dynamic_cast<const BooleanGeometry*>(&geometry);

    if(op_ != bg->op_)
        return op_ < bg->op_;
    // Compare children: null < non-null
    if(left_ && bg->left_) {
        if(*left_ != *(bg->left_))
            return *left_ < *(bg->left_);
    } else if(!left_ && bg->left_) {
        return true;
    } else if(left_ && !bg->left_) {
        return false;
    }
    if(right_ && bg->right_) {
        return *right_ < *(bg->right_);
    } else if(!right_ && bg->right_) {
        return true;
    }
    return false;
}

// ------------------------------------------------------------------------- //
void BooleanGeometry::print(std::ostream& os) const
{
    os << "BooleanGeometry: op=";
    switch(op_) {
        case BooleanOperation::UNION:        os << "UNION"; break;
        case BooleanOperation::SUBTRACTION:  os << "SUBTRACTION"; break;
        case BooleanOperation::INTERSECTION: os << "INTERSECTION"; break;
    }
    os << '\n';
    if(left_) {
        os << "  Left: " << *left_;
    }
    if(right_) {
        os << "  Right: " << *right_;
    }
}

// ------------------------------------------------------------------------- //
// CSG ray-casting algorithm.
//
// ComputeIntersections receives position/direction in BooleanGeometry's
// LOCAL coordinates (the base class Intersections() already transforms
// from global to local before calling this).
//
// The children's Intersections() method expects GLOBAL coordinates and
// handles their own Placement transforms internally.  From BooleanGeometry's
// perspective, its local frame IS the children's global frame, so we pass
// the local-frame position/direction directly to children's Intersections().
// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> BooleanGeometry::ComputeIntersections(
        siren::math::Vector3D const & position,
        siren::math::Vector3D const & direction) const
{
    // Tag to track which child produced each intersection
    struct TaggedIntersection {
        Intersection isect;
        int child; // 0 = left, 1 = right
    };

    std::vector<TaggedIntersection> all;

    // Collect intersections from left child
    if(left_) {
        std::vector<Intersection> left_hits = left_->Intersections(position, direction);
        for(auto & h : left_hits) {
            // Transform positions back to BooleanGeometry local coords
            // (children's Intersections() returns positions in the children's
            // global frame, which is our local frame, so no transform needed)
            TaggedIntersection ti;
            ti.isect = h;
            ti.child = 0;
            all.push_back(ti);
        }
    }

    // Collect intersections from right child
    if(right_) {
        std::vector<Intersection> right_hits = right_->Intersections(position, direction);
        for(auto & h : right_hits) {
            TaggedIntersection ti;
            ti.isect = h;
            ti.child = 1;
            all.push_back(ti);
        }
    }

    if(all.empty()) {
        return std::vector<Intersection>();
    }

    // Sort by distance
    std::sort(all.begin(), all.end(),
        [](TaggedIntersection const & a, TaggedIntersection const & b) {
            return a.isect.distance < b.isect.distance;
        });

    // Start with both children "outside" -- the state at t = -infinity.
    // The walk through all intersections (including negative-distance ones)
    // will build up the correct state naturally, producing proper enter/exit
    // pairs for the full line.
    bool in_left = false;
    bool in_right = false;

    // Now walk through all intersections in order
    bool was_inside_result = IsInsideResult(op_, in_left, in_right);

    std::vector<Intersection> result;

    for(auto const & ti : all) {
        // Update the appropriate child state
        if(ti.child == 0) {
            in_left = ti.isect.entering;
        } else {
            in_right = ti.isect.entering;
        }

        bool now_inside_result = IsInsideResult(op_, in_left, in_right);

        if(was_inside_result != now_inside_result) {
            Intersection out;
            out.distance = ti.isect.distance;
            out.position = ti.isect.position;
            out.hierarchy = ti.isect.hierarchy;
            out.entering = now_inside_result;
            result.push_back(out);
        }

        was_inside_result = now_inside_result;
    }

    return result;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> BooleanGeometry::ComputeDistanceToBorder(
        const siren::math::Vector3D& position,
        const siren::math::Vector3D& direction) const
{
    // CSG shapes are non-convex: there can be multiple enter/exit pairs.
    // Find the first positive exit (if inside) or the first positive
    // enter-then-exit pair (if outside).
    std::vector<Intersection> intersections = Intersections(position, direction);

    std::pair<double, double> distance;
    distance.first = -1;
    distance.second = -1;

    bool found_first = false;
    for(unsigned int i = 0; i < intersections.size(); ++i) {
        if(intersections[i].distance > GEOMETRY_PRECISION) {
            if(!found_first) {
                found_first = true;
                distance.first = intersections[i].distance;
                if(!intersections[i].entering) {
                    // First positive hit is an exit: we're inside
                    break;
                }
            } else {
                if(!intersections[i].entering) {
                    distance.second = intersections[i].distance;
                    break;
                }
            }
        }
    }

    return distance;
}

// ------------------------------------------------------------------------- //
AABB BooleanGeometry::GetBoundingBox() const
{
    // Children's bounding boxes are in their own local frames.
    // GetWorldBoundingBox() transforms them through their Placements
    // into the shared coordinate frame (which is BooleanGeometry's local frame).
    AABB left_box;
    AABB right_box;

    if(left_) {
        left_box = left_->GetWorldBoundingBox();
    }
    if(right_) {
        right_box = right_->GetWorldBoundingBox();
    }

    switch(op_) {
        case BooleanOperation::UNION:
        {
            // Envelope of both AABBs
            AABB result;
            if(left_) result.ExpandToInclude(left_box);
            if(right_) result.ExpandToInclude(right_box);
            return result;
        }
        case BooleanOperation::INTERSECTION:
        {
            // Intersection of both AABBs (clamped overlap region)
            if(!left_ || !right_) return AABB();
            double min_x = std::max(left_box.min_corner.GetX(), right_box.min_corner.GetX());
            double min_y = std::max(left_box.min_corner.GetY(), right_box.min_corner.GetY());
            double min_z = std::max(left_box.min_corner.GetZ(), right_box.min_corner.GetZ());
            double max_x = std::min(left_box.max_corner.GetX(), right_box.max_corner.GetX());
            double max_y = std::min(left_box.max_corner.GetY(), right_box.max_corner.GetY());
            double max_z = std::min(left_box.max_corner.GetZ(), right_box.max_corner.GetZ());
            // If no overlap, return an empty (invalid) AABB
            if(min_x > max_x || min_y > max_y || min_z > max_z) return AABB();
            return AABB(
                math::Vector3D(min_x, min_y, min_z),
                math::Vector3D(max_x, max_y, max_z)
            );
        }
        case BooleanOperation::SUBTRACTION:
        {
            // Conservative: subtraction cannot make the left shape larger
            if(!left_) return AABB();
            return left_box;
        }
    }
    return AABB();
}

} // namespace geometry
} // namespace siren

CEREAL_REGISTER_DYNAMIC_INIT(siren_BooleanGeometry);
