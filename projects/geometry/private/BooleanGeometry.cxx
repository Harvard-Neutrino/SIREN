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
    if(!bg) return false;

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
    // Collect child intersections into a tagged array.
    // Use a fixed-size stack buffer for typical cases; fall back to
    // heap allocation when deeply nested CSG produces many hits.
    struct TaggedHit {
        double distance;
        Intersection isect;
        int child; // 0 = left, 1 = right
    };

    static constexpr int STACK_CAPACITY = 32;
    TaggedHit stack_all[STACK_CAPACITY];
    std::vector<TaggedHit> heap_all;
    int n_all = 0;
    bool using_heap_all = false;

    auto add_tagged = [&](Intersection const & h, int child) {
        TaggedHit th;
        th.distance = h.distance;
        th.isect = h;
        th.child = child;
        if(!using_heap_all) {
            if(n_all < STACK_CAPACITY) {
                stack_all[n_all++] = th;
            } else {
                using_heap_all = true;
                heap_all.assign(stack_all, stack_all + STACK_CAPACITY);
                heap_all.push_back(th);
                n_all++;
            }
        } else {
            heap_all.push_back(th);
            n_all++;
        }
    };

    if(left_) {
        auto left_hits = left_->Intersections(position, direction);
        for(auto & h : left_hits) {
            add_tagged(h, 0);
        }
    }

    if(right_) {
        auto right_hits = right_->Intersections(position, direction);
        for(auto & h : right_hits) {
            add_tagged(h, 1);
        }
    }

    if(n_all == 0) return {};

    auto tag_cmp = [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    };
    if(using_heap_all) {
        std::sort(heap_all.begin(), heap_all.end(), tag_cmp);
    } else {
        std::sort(stack_all, stack_all + n_all, tag_cmp);
    }

    TaggedHit* all_ptr = using_heap_all ? heap_all.data() : stack_all;

    // Walk from t=-infinity: both children start "outside"
    bool in_left = false;
    bool in_right = false;
    bool was_inside = IsInsideResult(op_, in_left, in_right);

    Intersection stack_result[STACK_CAPACITY];
    std::vector<Intersection> heap_result;
    int n_result = 0;
    bool using_heap_result = false;

    auto add_result = [&](Intersection const & isect) {
        if(!using_heap_result) {
            if(n_result < STACK_CAPACITY) {
                stack_result[n_result++] = isect;
            } else {
                using_heap_result = true;
                heap_result.assign(stack_result, stack_result + STACK_CAPACITY);
                heap_result.push_back(isect);
                n_result++;
            }
        } else {
            heap_result.push_back(isect);
            n_result++;
        }
    };

    for(int i = 0; i < n_all; ++i) {
        if(all_ptr[i].child == 0) {
            in_left = all_ptr[i].isect.entering;
        } else {
            in_right = all_ptr[i].isect.entering;
        }

        bool now_inside = IsInsideResult(op_, in_left, in_right);

        if(was_inside != now_inside) {
            Intersection r;
            r.distance = all_ptr[i].isect.distance;
            r.position = all_ptr[i].isect.position;
            r.hierarchy = all_ptr[i].isect.hierarchy;
            r.entering = now_inside;
            add_result(r);
        }

        was_inside = now_inside;
    }

    if(using_heap_result) {
        return heap_result;
    }
    return {stack_result, stack_result + n_result};
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
