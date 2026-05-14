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
        siren::math::Vector3D const & direction) const {
    // Collect child intersections into a fixed-size tagged array.
    // Each child produces at most ~8 intersections for typical shapes;
    // 32 total handles deeply nested CSG.
    struct TaggedHit {
        double distance;
        Intersection isect;
        int child; // 0 = left, 1 = right
    };

    TaggedHit all[32];
    int n_all = 0;

    if(left_) {
        auto left_hits = left_->Intersections(position, direction);
        for(auto & h : left_hits) {
            if(n_all < 32) {
                all[n_all].distance = h.distance;
                all[n_all].isect = h;
                all[n_all].child = 0;
                n_all++;
            }
        }
    }

    if(right_) {
        auto right_hits = right_->Intersections(position, direction);
        for(auto & h : right_hits) {
            if(n_all < 32) {
                all[n_all].distance = h.distance;
                all[n_all].isect = h;
                all[n_all].child = 1;
                n_all++;
            }
        }
    }

    if(n_all == 0) return {};

    std::sort(all, all + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // Walk from t=-infinity: both children start "outside"
    bool in_left = false;
    bool in_right = false;
    bool was_inside = IsInsideResult(op_, in_left, in_right);

    Intersection result[32];
    int n_result = 0;

    for(int i = 0; i < n_all; ++i) {
        if(all[i].child == 0) {
            in_left = all[i].isect.entering;
        } else {
            in_right = all[i].isect.entering;
        }

        bool now_inside = IsInsideResult(op_, in_left, in_right);

        if(was_inside != now_inside && n_result < 32) {
            result[n_result].distance = all[i].isect.distance;
            result[n_result].position = all[i].isect.position;
            result[n_result].hierarchy = all[i].isect.hierarchy;
            result[n_result].entering = now_inside;
            n_result++;
        }

        was_inside = now_inside;
    }

    return {result, result + n_result};
}

// ------------------------------------------------------------------------- //
AABB BooleanGeometry::GetBoundingBox() const {
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
