#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Geometry.h"

using namespace earthmodel;

Path::Path() {

}

Path::Path(std::shared_ptr<const EarthModel> earth_model) {
    SetEarthModel(earth_model);
}

Path::Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & last_point) {
    SetEarthModel(earth_model);
    SetPoints(first_point, last_point);
}

Path::Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & direction, double distance) {
    SetEarthModel(earth_model);
    SetPointsWithRay(first_point, direction, distance);
}

bool Path::HasEarthModel() {
    return set_earth_model_;
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

std::shared_ptr<const EarthModel> Path::GetEarthModel() {
    return earth_model_;
}

Vector3D Path::GetFirstPoint() {
    return first_point_;
}

Vector3D Path::GetLastPoint() {
    return last_point_;
}

Vector3D Path::GetDirection() {
    return direction_;
}

double Path::GetDistance() {
    return distance_;
}

Geometry::IntersectionList Path::GetIntersections() {
    return intersections_;
}

void Path::SetEarthModel(std::shared_ptr<const EarthModel> earth_model) {
    earth_model_ = earth_model;
    set_earth_model_ = true;
}

void Path::EnsureEarthModel() {
    if(not set_earth_model_) {
        throw(std::runtime_error("Earth model not set!"));
    }
}

void Path::SetPoints(Vector3D first_point, Vector3D last_point) {
    first_point_ = first_point;
    last_point_ = last_point;
    direction_ = last_point_ - first_point_;
    distance_ = direction_.magnitude();
    direction_.normalize();
    set_points_ = true;
    set_intersections_ = false;
    set_column_depth_ = false;
}

void Path::SetPointsWithRay(Vector3D first_point, Vector3D direction, double distance) {
    first_point_ = first_point;
    direction_ = direction;
    direction_.normalize();
    //double dif = std::abs(direction_.magnitude() - direction.magnitude()) / std::max(direction_.magnitude(), direction.magnitude());
    //if(not std::isnan(dif)) assert(dif < 1e-12);
    distance_ = distance;
    last_point_ = first_point + direction * distance;
    set_points_ = true;
    set_intersections_ = false;
    set_column_depth_ = false;
}

void Path::EnsurePoints() {
    if(not set_points_) {
        throw(std::runtime_error("Points not set!"));
    }
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

void Path::ClipToOuterBounds() {
    EnsureIntersections();
    EnsurePoints();
    Geometry::IntersectionList bounds = EarthModel::GetOuterBounds(intersections_);
    if(bounds.intersections.size() > 0) {
        assert(bounds.intersections.size() == 2);
        Vector3D p0 = bounds.intersections[0].position;
        Vector3D p1 = bounds.intersections[1].position;
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        double dot = direction_ * direction;
        assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
        if(dot < 0) {
            p0.swap(p1);
        }
        bool clip_0 = (p0 - first_point_) * direction_ > 0;
        bool clip_1 = (p1 - last_point_) * direction_ < 0;
        bool clip = clip_0 or clip_1;
        if(clip_0) {
            first_point_ = p0;
        }
        if(clip_1) {
            last_point_ = p1;
        }
        if(clip) {
            distance_ = (last_point_ - first_point_).magnitude();
            set_column_depth_ = false;
        }
    } else {
        return;
    }
}

void Path::Flip() {
    EnsurePoints();
    std::swap(first_point_, last_point_);
    direction_ *= -1;
}


////
// Extend / Shrink By Distance
////
void Path::ExtendFromEndByDistance(double distance) {
    distance_ += distance;
    last_point_ += direction_ * distance;
    if(distance_ < 0) {
        distance_ = 0;
        last_point_ = first_point_;
    }
    set_column_depth_ = false;
}

void Path::ExtendFromStartByDistance(double distance) {
    distance_ += distance;
    first_point_ += direction_ * -distance;
    if(distance_ < 0) {
        distance_ = 0;
        first_point_ = last_point_;
    }
    set_column_depth_ = false;
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
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double distance = GetDistanceFromEndAlongPath(interaction_depth, targets, total_cross_sections);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromStartByInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double distance = GetDistanceFromStartInReverse(interaction_depth, targets, total_cross_sections);
    ExtendFromStartByDistance(distance);
}

void Path::ShrinkFromEndByInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double distance = GetDistanceFromEndInReverse(interaction_depth, targets, total_cross_sections);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromStartByInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double distance = GetDistanceFromStartAlongPath(interaction_depth, targets, total_cross_sections);
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
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double shift = interaction_depth - GetInteractionDepthInBounds(targets, total_cross_sections);
    if(shift > 0) {
        ExtendFromEndByInteractionDepth(shift, targets, total_cross_sections);
    }
}

void Path::ExtendFromStartToInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double shift = interaction_depth - GetInteractionDepthInBounds(targets, total_cross_sections);
    if(shift > 0) {
        ExtendFromStartByInteractionDepth(shift, targets, total_cross_sections);
    }
}

void Path::ShrinkFromEndToInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double shift = GetInteractionDepthInBounds(targets, total_cross_sections) - interaction_depth;
    if(shift > 0) {
        ShrinkFromEndByInteractionDepth(shift, targets, total_cross_sections);
    }
}

void Path::ShrinkFromStartToInteractionDepth(double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) {
    double shift = GetInteractionDepthInBounds(targets, total_cross_sections) - interaction_depth;
    if(shift > 0) {
        ShrinkFromStartByInteractionDepth(shift, targets, total_cross_sections);
    }
}


////
// Get ColumnDepth
////
double Path::GetColumnDepthInBounds() {
    EnsureIntersections();
    if(HasColumnDepth()) {
        return column_depth_cached_;
    } else {
        double column_depth = earth_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_);
        column_depth_cached_ = column_depth;
        return column_depth;
    }
}


////
// Get InteractionDepth
////
double Path::GetInteractionDepthInBounds(
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double interaction_depth = earth_model_->GetInteractionDepthInCGS(intersections_, first_point_, last_point_, targets, total_cross_sections);
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
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance);
}

double Path::GetColumnDepthFromEndInBounds(double distance) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance);
}

double Path::GetColumnDepthFromStartAlongPath(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance), distance);
}

double Path::GetColumnDepthFromEndAlongPath(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * distance), distance);
}

double Path::GetColumnDepthFromStartInReverse(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * -distance), distance);
}

double Path::GetColumnDepthFromEndInReverse(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance), distance);
}


////
// Get InteractionDepth From
////
double Path::GetInteractionDepthFromStartInBounds(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetInteractionDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, targets, total_cross_sections);
}

double Path::GetInteractionDepthFromEndInBounds(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetInteractionDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, targets, total_cross_sections);
}

double Path::GetInteractionDepthFromStartAlongPath(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetInteractionDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, targets, total_cross_sections), distance);
}

double Path::GetInteractionDepthFromEndAlongPath(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetInteractionDepthInCGS(intersections_, last_point_, last_point_ + direction_ * distance, targets, total_cross_sections), distance);
}

double Path::GetInteractionDepthFromStartInReverse(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetInteractionDepthInCGS(intersections_, first_point_, first_point_ + direction_ * -distance, targets, total_cross_sections), distance);
}

double Path::GetInteractionDepthFromEndInReverse(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetInteractionDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, targets, total_cross_sections), distance);
}


////
// Get Distance From (ColumnDepth)
///
double Path::GetDistanceFromStartInBounds(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, -direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    return distance;
}


////
// Get Distance From (InteractionDepth)
///
double Path::GetDistanceFromStartInBounds(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, direction_, interaction_depth, targets, total_cross_sections);
    if(distance > distance_) {
        distance = distance_;
    } else if(interaction_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, -direction_, interaction_depth, targets, total_cross_sections);
    if(distance > distance_) {
        distance = distance_;
    } else if(interaction_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, direction_, interaction_depth, targets, total_cross_sections);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, direction_, interaction_depth, targets, total_cross_sections);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, first_point_, -direction_, interaction_depth, targets, total_cross_sections);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForInteractionDepthFromPoint(intersections_, last_point_, -direction_, interaction_depth, targets, total_cross_sections);
    return distance;
}


bool Path::IsWithinBounds(Vector3D point) {
    EnsurePoints();
    double d0 = earthmodel::scalar_product(direction_, first_point_ - point);
    double d1 = earthmodel::scalar_product(direction_, last_point_ - point);
    return d0 <= 0 and d1 >= 0;
}

