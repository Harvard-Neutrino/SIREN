#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "earthmodel-service/Path.h"

using namespace earthmodel;

Path::Path() {

}

Path::Path(std::shared_ptr<const EarthModel> earth_model) {
    SetEarthModel(earth_model);
}

void Path::SetEarthModel(std::shared_ptr<const EarthModel> earth_model) {
    earth_model_ = earth_model;
    set_earth_model_ = true;
}

void Path::EnsureEarthModel() {
    assert(set_earth_model_);
}

void Path::SetPoints(Vector3D first_point, Vector3D last_point) {
    first_point_ = first_point;
    last_point_ = last_point;
    direction_ = last_point_ - first_point_;
    distance_ = direction_.magnitude();
    direction_.normalize();
    set_points_ = true;
}

void Path::SetPointsWithRay(Vector3D first_point, Vector3D direction, double distance) {
    first_point_ = first_point;
    direction_ = direction;
    direction_.normalize();
    assert(direction_.magnitude() == direction.magnitude());
    distance_ = distance;
    last_point_ = first_point + direction * distance;
    distance_ = direction.magnitude();
    set_points_ = true;
}

void Path::EnsurePoints() {
    assert(set_points_);
}

void Path::SetIntersections(Geometry::IntersectionList const & intersections) {
    intersections_ = intersections;
    set_intersections_ = true;
}

void Path::ComputeIntersections() {
    EnsureEarthModel();
    EnsurePoints();
    intersections_ = earth_model_->GetIntersections(first_point_, direction_);
    set_intersections_ = true;
}

void Path::EnsureIntersections() {
    if(not set_intersections_) {
        ComputeIntersections();
    }
}

void Path::Flip() {
    EnsurePoints();
    std::swap(first_point_, last_point_);
    direction_ *= -1;
}

void Path::ExtendFromEndByDistance(double distance) {
    distance_ += distance;
    last_point_ += direction_ * distance;
    if(distance < 0) {
        distance = 0;
        last_point_ = first_point_;
    }
    set_column_depth_nucleon_ = false;
    set_column_depth_electron_ = false;
}

void Path::ExtendFromEndByColumnDepth(double column_depth, bool use_electron_density) {
    double distance = GetDistanceFromEndAlongPath(column_depth, use_electron_density);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromStartByDistance(double distance) {
    distance_ -= distance;
    first_point_ += direction_ * -distance;
    if(distance < 0) {
        distance = 0;
        first_point_ = last_point_;
    }
    set_column_depth_nucleon_ = false;
    set_column_depth_electron_ = false;
}

void Path::ExtendFromStartByColumnDepth(double column_depth, bool use_electron_density) {
    double distance = GetDistanceFromStartInReverse(column_depth, use_electron_density);
    ExtendFromStartByDistance(distance);
}

void Path::ShrinkFromEndByDistance(double distance) {
    ExtendFromEndByDistance(-distance);
}

void Path::ShrinkFromEndByColumnDepth(double column_depth, bool use_electron_density) {
    double distance = GetDistanceFromEndInReverse(column_depth, use_electron_density);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromStartByDistance(double distance) {
    ExtendFromStartByDistance(-distance);
}

void Path::ShrinkFromStartByColumnDepth(double column_depth, bool use_electron_density) {
    double distance = GetDistanceFromStartAlongPath(column_depth, use_electron_density);
    ShrinkFromStartByDistance(distance);
}

void Path::ExtendFromEndToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromEndByDistance(shift);
    }
}

void Path::ExtendFromEndToColumnDepth(double column_depth, bool use_electron_density) {
    double shift = column_depth - GetColumnDepthInBounds(use_electron_density);
    if(shift > 0) {
        ExtendFromEndByColumnDepth(shift);
    }
}

void Path::ExtendFromStartToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromStartByDistance(shift);
    }
}

void Path::ExtendFromStartToColumnDepth(double column_depth, bool use_electron_density) {
    double shift = column_depth - GetColumnDepthInBounds(use_electron_density);
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

void Path::ShrinkFromEndToColumnDepth(double column_depth, bool use_electron_density) {
    double shift = GetColumnDepthInBounds(use_electron_density) - column_depth;
    if(shift > 0) {
        ShrinkFromEndByColumnDepth(shift);
    }
}

void Path::ShrinkFromStartToDistance(double distance) {
    double shift = distance_ - distance;
    if(shift > 0) {
        ShrinkFromStartByDistance(shift);
    }
}

void Path::ShrinkFromStartToColumnDepth(double column_depth, bool use_electron_density) {
    double shift = GetColumnDepthInBounds(use_electron_density) - column_depth;
    if(shift > 0) {
        ShrinkFromStartByDistance(shift);
    }
}

double Path::GetColumnDepthInBounds(bool use_electron_density) {
    if((not use_electron_density) and (not set_column_depth_nucleon_)) {
        column_depth_nucleon_ = earth_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_, use_electron_density);
        set_column_depth_nucleon_ = true;
    } else if(use_electron_density and (not set_column_depth_electron_)) {
        column_depth_electron_ = earth_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_, use_electron_density);
        set_column_depth_electron_ = true;
    }

    if(not use_electron_density) {
        return column_depth_nucleon_;
    } else {
        return column_depth_electron_;
    }
}

double Path::GetColumnDepthFromStartInBounds(double distance, bool use_electron_density) {
    if(distance > distance_) {
        distance = distance_;
    }
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, use_electron_density);
}

double Path::GetColumnDepthFromEndInBounds(double distance, bool use_electron_density) {
    if(distance > distance_) {
        distance = distance_;
    }
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, use_electron_density);
}

double Path::GetColumnDepthFromStartAlongPath(double distance, bool use_electron_density) {
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, use_electron_density);
}

double Path::GetColumnDepthFromEndAlongPath(double distance, bool use_electron_density) {
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * distance, use_electron_density);
}

double Path::GetColumnDepthFromStartInReverse(double distance, bool use_electron_density) {
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * -distance, use_electron_density);
}

double Path::GetColumnDepthFromEndInReverse(double distance, bool use_electron_density) {
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, use_electron_density);
}

double Path::GetDistanceFromStartInBounds(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth, use_electron_density);
    if(distance > distance_) {
        distance = distance_;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth, use_electron_density);
    if(distance > distance_) {
        distance = distance_;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth, use_electron_density);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, direction_, column_depth, use_electron_density);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, -direction_, column_depth, use_electron_density);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double column_depth, bool use_electron_density) {
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth, use_electron_density);
    return distance;
}

