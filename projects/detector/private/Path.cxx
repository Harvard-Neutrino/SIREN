#include "SIREN/detector/Path.h"

#include <cmath>                                  // for copysign, abs
#include <vector>                                 // for vector
#include <utility>                                // for swap
#include <assert.h>                               // for assert
#include <stdlib.h>                               // for abs
#include <algorithm>                              // for max

#include "SIREN/dataclasses/Particle.h"  // for Particle
#include "SIREN/detector/DetectorModel.h"   // for DetectorModel
#include "SIREN/geometry/Geometry.h"     // for Geometry

namespace siren {
namespace detector {

void Path::UpdatePoints() {
    if((not set_points_) and set_det_points_ and set_detector_model_) {
        first_point_ = detector_model_->ToGeo(first_point_det_);
        last_point_ = detector_model_->ToGeo(last_point_det_);
        direction_ = detector_model_->ToGeo(direction_det_);
        set_points_ = true;
    } else if ((not set_det_points_) and set_points_ and set_detector_model_) {
        first_point_det_ = detector_model_->ToDet(first_point_);
        last_point_det_ = detector_model_->ToDet(last_point_);
        direction_det_ = detector_model_->ToDet(direction_);
        set_det_points_ = true;
    }
}

bool Path::IsInfinite(siren::math::Vector3D const & vec) {
    return std::isinf(vec.GetX()) or std::isinf(vec.GetY()) or std::isinf(vec.GetZ());
}

void Path::RequireFirstFinite() {
    if(first_inf_)
        throw std::runtime_error("First point is required to be finite here");
}

void Path::RequireLastFinite() {
    if(last_inf_)
        throw std::runtime_error("Last point is required to be finite here");
}

void Path::RequireBothFinite() {
    if(first_inf_ or last_inf_)
        throw std::runtime_error("Both points are required to be finite here");
}

void Path::RequireOneFinite() {
    if(first_inf_ and last_inf_)
        throw std::runtime_error("At least one point is required to be finite here");
}

Path::Path() {

}

Path::Path(std::shared_ptr<const DetectorModel> detector_model) {
    SetDetectorModel(detector_model);
}

Path::Path(std::shared_ptr<const DetectorModel> detector_model, GeometryPosition const & first_point, GeometryPosition const & last_point) {
    SetDetectorModel(detector_model);
    SetPoints(first_point, last_point);
}

Path::Path(std::shared_ptr<const DetectorModel> detector_model, GeometryPosition const & first_point, GeometryDirection const & direction, double distance) {
    SetDetectorModel(detector_model);
    SetPointsWithRay(first_point, direction, distance);
}

Path::Path(std::shared_ptr<const DetectorModel> detector_model, DetectorPosition const & first_point, DetectorPosition const & last_point) {
    SetDetectorModel(detector_model);
    SetPoints(first_point, last_point);
}

Path::Path(std::shared_ptr<const DetectorModel> detector_model, DetectorPosition const & first_point, DetectorDirection const & direction, double distance) {
    SetDetectorModel(detector_model);
    SetPointsWithRay(first_point, direction, distance);
}

bool Path::HasDetectorModel() {
    return set_detector_model_;
}

bool Path::HasPoints() {
    return set_points_;
}

bool Path::HasIntersections() {
    return set_intersections_;
}

bool Path::HasColumnDepth() {
    return set_column_depth_;
}

std::shared_ptr<const DetectorModel> Path::GetDetectorModel() {
    return detector_model_;
}

DetectorPosition const & Path::GetFirstPoint() {
    UpdatePoints();
    return first_point_det_;
}

DetectorPosition const & Path::GetLastPoint() {
    UpdatePoints();
    return last_point_det_;
}

DetectorDirection const & Path::GetDirection() {
    UpdatePoints();
    return direction_det_;
}

GeometryPosition const & Path::GetGeoFirstPoint() {
    UpdatePoints();
    return first_point_;
}

GeometryPosition const & Path::GetGeoLastPoint() {
    UpdatePoints();
    return last_point_;
}

GeometryDirection const & Path::GetGeoDirection() {
    UpdatePoints();
    return direction_;
}

double Path::GetDistance() {
    return distance_;
}

geometry::Geometry::IntersectionList const & Path::GetIntersections() {
    return intersections_;
}

void Path::SetDetectorModel(std::shared_ptr<const DetectorModel> detector_model) {
    if(set_detector_model_ and set_det_points_) {
        set_points_ = false;
    }
    detector_model_ = detector_model;
    set_detector_model_ = true;
    UpdatePoints();
}

void Path::EnsureDetectorModel() {
    if(not set_detector_model_) {
        throw(std::runtime_error("Detector model not set!"));
    }
}

void Path::SetPoints(GeometryPosition first_point, GeometryPosition last_point) {
    first_point_ = first_point;
    last_point_ = last_point;
    direction_ = GeometryDirection(last_point_ - first_point_);
    distance_ = direction_->magnitude();
    direction_->normalize();
    set_points_ = true;
    set_det_points_ = false;
    set_intersections_ = false;
    set_column_depth_ = false;
    first_inf_ = IsInfinite(first_point);
    last_inf_ = IsInfinite(last_point);
    RequireBothFinite();
    UpdatePoints();
}

void Path::SetPoints(DetectorPosition first_point, DetectorPosition last_point) {
    first_point_det_ = first_point;
    last_point_det_ = last_point;
    direction_det_ = DetectorDirection(last_point_det_ - first_point_det_);
    distance_ = direction_det_->magnitude();
    direction_det_->normalize();
    set_points_ = false;
    set_det_points_ = true;
    set_intersections_ = false;
    set_column_depth_ = false;
    first_inf_ = IsInfinite(first_point);
    last_inf_ = IsInfinite(last_point);
    RequireBothFinite();
    UpdatePoints();
}

void Path::SetPointsWithRay(GeometryPosition first_point, GeometryDirection direction, double distance) {
    first_point_ = first_point;
    direction_ = direction;
    direction_->normalize();
    //double dif = std::abs(direction_.magnitude() - direction.magnitude()) / std::max(direction_.magnitude(), direction.magnitude());
    //if(not std::isnan(dif)) assert(dif < 1e-12);
    distance_ = distance;
    last_point_ = first_point + GeometryPosition(direction * distance);
    set_points_ = true;
    set_det_points_ = false;
    set_intersections_ = false;
    set_column_depth_ = false;
    first_inf_ = IsInfinite(first_point_); // Set using GeometryPosition
    last_inf_ = IsInfinite(last_point_); // Set using GeometryPosition
    RequireFirstFinite();
    UpdatePoints();
}

void Path::SetPointsWithRay(DetectorPosition first_point, DetectorDirection direction, double distance) {
    first_point_det_ = first_point;
    direction_det_ = direction;
    direction_det_->normalize();
    //double dif = std::abs(direction_.magnitude() - direction.magnitude()) / std::max(direction_.magnitude(), direction.magnitude());
    //if(not std::isnan(dif)) assert(dif < 1e-12);
    distance_ = distance;
    last_point_det_ = first_point + DetectorPosition(direction * distance);
    set_points_ = false;
    set_det_points_ = true;
    set_intersections_ = false;
    set_column_depth_ = false;
    first_inf_ = IsInfinite(first_point_det_); // Set using DetectorPosition
    last_inf_ = IsInfinite(last_point_det_); // Set using DetectorPosition
    RequireFirstFinite();
    UpdatePoints();
}

void Path::EnsurePoints() {
    UpdatePoints();
    if(not set_points_) {
        throw(std::runtime_error("Points not set!"));
    }
}

void Path::SetIntersections(geometry::Geometry::IntersectionList const & intersections) {
    intersections_ = intersections;
    set_intersections_ = true;
}

void Path::ComputeIntersections() {
    EnsureDetectorModel();
    EnsurePoints();
    intersections_ = detector_model_->GetIntersections(first_point_, direction_);
    set_intersections_ = true;
}

void Path::EnsureIntersections() {
    if(not set_intersections_) {
        ComputeIntersections();
    }
}

void Path::ClipToOuterBounds() {
    EnsureIntersections();
    EnsurePoints();
    geometry::Geometry::IntersectionList bounds = DetectorModel::GetOuterBounds(intersections_);
    if(bounds.intersections.size() > 0) {
        assert(bounds.intersections.size() == 2);
        math::Vector3D p0 = bounds.intersections[0].position;
        math::Vector3D p1 = bounds.intersections[1].position;
        math::Vector3D direction = p1 - p0;
        direction.normalize();
        double dot = direction_ * direction;
        assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
        if(dot < 0) {
            p0.swap(p1);
        }
        bool clip_0 = (first_inf_) || ((p0 - first_point_) * direction_ > 0);
        bool clip_1 = (last_inf_) || ((p1 - last_point_) * direction_ < 0);
        bool clip = clip_0 or clip_1;
        if(clip_0) {
            first_point_ = GeometryPosition(p0);
            first_inf_ = IsInfinite(first_point_); // re-check infinite
        }
        if(clip_1) {
            last_point_ = GeometryPosition(p1);
            last_inf_ = IsInfinite(last_point_); // re-check infinite
        }
        if(clip) {
            distance_ = (last_point_ - first_point_)->magnitude();
            set_column_depth_ = false;
        }
        set_det_points_ = false;
    } else {
        return;
    }
}

void Path::Flip() {
    std::swap(first_point_, last_point_);
    std::swap(first_point_det_, last_point_det_);
    std::swap(first_inf_, last_inf_);
    direction_ *= -1;
    direction_det_ *= -1;
}


////
// Extend / Shrink By Distance
////
void Path::ExtendFromEndByDistance(double distance) {
    EnsurePoints();
    RequireLastFinite();
    distance_ += distance;
    last_point_ += direction_ * distance;
    if(distance_ < 0) {
        distance_ = 0;
        last_point_ = first_point_;
    }
    set_column_depth_ = false;
    set_det_points_ = false;
}

void Path::ExtendFromStartByDistance(double distance) {
    EnsurePoints();
    RequireFirstFinite();
    distance_ += distance;
    first_point_ += direction_ * -distance;
    if(distance_ < 0) {
        distance_ = 0;
        first_point_ = last_point_;
    }
    set_column_depth_ = false;
    set_det_points_ = false;
}

void Path::ShrinkFromEndByDistance(double distance) {
    ExtendFromEndByDistance(-distance);
}

void Path::ShrinkFromStartByDistance(double distance) {
    ExtendFromStartByDistance(-distance);
}


////
// Extend / Shrink By ColumnDepth
////
void Path::ExtendFromEndByColumnDepth(double column_depth) {
    double distance = GetDistanceFromEndAlongPath(column_depth);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromStartByColumnDepth(double column_depth) {
    double distance = GetDistanceFromStartInReverse(column_depth);
    ExtendFromStartByDistance(distance);
}

void Path::ShrinkFromEndByColumnDepth(double column_depth) {
    double distance = GetDistanceFromEndInReverse(column_depth);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromStartByColumnDepth(double column_depth) {
    double distance = GetDistanceFromStartAlongPath(column_depth);
    ShrinkFromStartByDistance(distance);
}


////
// Extend / Shrink By InteractionDepth
////
void Path::ExtendFromEndByInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double distance = GetDistanceFromEndAlongPath(interaction_depth, targets, total_cross_sections, total_decay_length);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromStartByInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double distance = GetDistanceFromStartInReverse(interaction_depth, targets, total_cross_sections, total_decay_length);
    ExtendFromStartByDistance(distance);
}

void Path::ShrinkFromEndByInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double distance = GetDistanceFromEndInReverse(interaction_depth, targets, total_cross_sections, total_decay_length);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromStartByInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double distance = GetDistanceFromStartAlongPath(interaction_depth, targets, total_cross_sections, total_decay_length);
    ShrinkFromStartByDistance(distance);
}


////
// Extend / Shrink To Distance
////
void Path::ExtendFromEndToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromEndByDistance(shift);
    }
}

void Path::ExtendFromStartToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromStartByDistance(shift);
    }
}

void Path::ShrinkFromEndToDistance(double distance) {
    double shift = distance_ - distance;
    if(shift > 0) {
        ShrinkFromEndByDistance(shift);
    }
}

void Path::ShrinkFromStartToDistance(double distance) {
    double shift = distance_ - distance;
    if(shift > 0) {
        ShrinkFromStartByDistance(shift);
    }
}


////
// Extend / Shrink To ColumnDepth
////
void Path::ExtendFromEndToColumnDepth(double column_depth) {
    double shift = column_depth - GetColumnDepthInBounds();
    if(shift > 0) {
        ExtendFromEndByColumnDepth(shift);
    }
}

void Path::ExtendFromStartToColumnDepth(double column_depth) {
    double shift = column_depth - GetColumnDepthInBounds();
    if(shift > 0) {
        ExtendFromStartByColumnDepth(shift);
    }
}

void Path::ShrinkFromEndToColumnDepth(double column_depth) {
    double shift = GetColumnDepthInBounds() - column_depth;
    if(shift > 0) {
        ShrinkFromEndByColumnDepth(shift);
    }
}

void Path::ShrinkFromStartToColumnDepth(double column_depth) {
    double shift = GetColumnDepthInBounds() - column_depth;
    if(shift > 0) {
        ShrinkFromStartByColumnDepth(shift);
    }
}


////
// Extend / Shrink To InteractionDepth
////
void Path::ExtendFromEndToInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double shift = interaction_depth - GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);
    if(shift > 0) {
        ExtendFromEndByInteractionDepth(shift, targets, total_cross_sections, total_decay_length);
    }
}

void Path::ExtendFromStartToInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double shift = interaction_depth - GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length);
    if(shift > 0) {
        ExtendFromStartByInteractionDepth(shift, targets, total_cross_sections, total_decay_length);
    }
}

void Path::ShrinkFromEndToInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double shift = GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length) - interaction_depth;
    if(shift > 0) {
        ShrinkFromEndByInteractionDepth(shift, targets, total_cross_sections, total_decay_length);
    }
}

void Path::ShrinkFromStartToInteractionDepth(double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) {
    double shift = GetInteractionDepthInBounds(targets, total_cross_sections, total_decay_length) - interaction_depth;
    if(shift > 0) {
        ShrinkFromStartByInteractionDepth(shift, targets, total_cross_sections, total_decay_length);
    }
}


////
// Get ColumnDepth
////
double Path::GetColumnDepthInBounds() {
    EnsureIntersections();
    EnsurePoints();
    RequireBothFinite();
    if(HasColumnDepth()) {
        return column_depth_cached_;
    } else {
        double column_depth = detector_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_);
        column_depth_cached_ = column_depth;
        return column_depth;
    }
}


////
// Get InteractionDepth
////
double Path::GetInteractionDepthInBounds(
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireBothFinite();
    double interaction_depth = detector_model_->GetInteractionDepthInCGS(intersections_, first_point_, last_point_, targets, total_cross_sections, total_decay_length);
    return interaction_depth;
}


////
// Get ColumnDepth From
////
double Path::GetColumnDepthFromStartInBounds(double distance) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return detector_model_->GetColumnDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * distance));
}

double Path::GetColumnDepthFromEndInBounds(double distance) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return detector_model_->GetColumnDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * -distance));
}

double Path::GetColumnDepthFromStartAlongPath(double distance) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return std::copysign(detector_model_->GetColumnDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * distance)), distance);
}

double Path::GetColumnDepthFromEndAlongPath(double distance) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return std::copysign(detector_model_->GetColumnDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * distance)), distance);
}

double Path::GetColumnDepthFromStartInReverse(double distance) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return std::copysign(detector_model_->GetColumnDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * -distance)), distance);
}

double Path::GetColumnDepthFromEndInReverse(double distance) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return std::copysign(detector_model_->GetColumnDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * -distance)), distance);
}


////
// Get InteractionDepth From
////
double Path::GetInteractionDepthFromStartInBounds(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return detector_model_->GetInteractionDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * distance), targets, total_cross_sections, total_decay_length);
}

double Path::GetInteractionDepthFromEndInBounds(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return detector_model_->GetInteractionDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * -distance), targets, total_cross_sections, total_decay_length);
}

double Path::GetInteractionDepthFromStartAlongPath(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return std::copysign(detector_model_->GetInteractionDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * distance), targets, total_cross_sections, total_decay_length), distance);
}

double Path::GetInteractionDepthFromEndAlongPath(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return std::copysign(detector_model_->GetInteractionDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * distance), targets, total_cross_sections, total_decay_length), distance);
}

double Path::GetInteractionDepthFromStartInReverse(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    return std::copysign(detector_model_->GetInteractionDepthInCGS(intersections_, first_point_, GeometryPosition(first_point_ + direction_ * -distance), targets, total_cross_sections, total_decay_length), distance);
}

double Path::GetInteractionDepthFromEndInReverse(double distance,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    return std::copysign(detector_model_->GetInteractionDepthInCGS(intersections_, last_point_, GeometryPosition(last_point_ + direction_ * -distance), targets, total_cross_sections, total_decay_length), distance);
}


////
// Get Distance From (ColumnDepth)
///
double Path::GetDistanceFromStartInBounds(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, -direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double column_depth) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    return distance;
}


////
// Get Distance From (InteractionDepth)
///
double Path::GetDistanceFromStartInBounds(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    if(distance > distance_) {
        distance = distance_;
    } else if(interaction_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, -direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    if(distance > distance_) {
        distance = distance_;
    } else if(interaction_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireFirstFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, -direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) {
    EnsureIntersections();
    EnsurePoints();
    RequireLastFinite();
    double distance = detector_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, -direction_, interaction_depth, targets, total_cross_sections, total_decay_length);
    return distance;
}


bool Path::IsWithinBounds(GeometryPosition point) {
    UpdatePoints();
    RequireBothFinite();
    if(set_points_) {
        double d0 = siren::math::scalar_product(direction_, first_point_ - point);
        double d1 = siren::math::scalar_product(direction_, last_point_ - point);
        return d0 <= 0 and d1 >= 0;
    } else {
        EnsurePoints();
        return false;
    }
}

double Path::GetDistanceFromStartInBounds(GeometryPosition point) {
    UpdatePoints();
    RequireFirstFinite();
    if(set_points_) {
        double d0 = siren::math::scalar_product(direction_, point - first_point_);
        return std::max(0.0, d0);
    } else {
        EnsurePoints();
        return false;
    }
}

bool Path::IsWithinBounds(DetectorPosition point) {
    UpdatePoints();
    RequireBothFinite();
    if(set_det_points_) {
        double d0 = siren::math::scalar_product(direction_det_, first_point_det_ - point);
        double d1 = siren::math::scalar_product(direction_det_, last_point_det_ - point);
        return d0 <= 0 and d1 >= 0;
    } else if(set_points_ and set_detector_model_) {
        return IsWithinBounds(detector_model_->ToGeo(point));
    } else {
        throw(std::runtime_error("Detector points not set!"));
    }
}

double Path::GetDistanceFromStartInBounds(DetectorPosition point) {
    UpdatePoints();
    RequireFirstFinite();
    if(set_det_points_) {
        double d0 = siren::math::scalar_product(direction_det_, point - first_point_det_);
        return std::max(0.0, d0);
    } else if (set_points_ and set_detector_model_) {
        return GetDistanceFromStartInBounds(detector_model_->ToGeo(point));
    } else {
        throw(std::runtime_error("Detector points not set!"));
    }
}

} // namespace detector
} // namespace siren
