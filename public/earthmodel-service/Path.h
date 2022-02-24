#pragma once
#ifndef LI_Path_H
#define LI_Path_H

#include <memory>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>
#include "serialization/array.h"

#include "earthmodel-service/EarthModel.h"

namespace earthmodel {

class Path {
private:
    std::shared_ptr<const EarthModel> earth_model_;
    bool set_earth_model_ = false;

    Vector3D first_point_;
    Vector3D last_point_;
    Vector3D direction_;
    double distance_ = 0;
    bool set_points_ = false;

    double column_depth_cached_;
    bool set_column_depth_ = false;
    std::map<std::set<LeptonInjector::Particle::ParticleType>, double> column_depth_cache_;

    Geometry::IntersectionList intersections_;
    bool set_intersections_ = false;
public:
    Path();
    Path(std::shared_ptr<const EarthModel> earth_model);
    Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & last_point);
    Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & direction, double distance);

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("EarthModel", earth_model_));
            archive(cereal::make_nvp("SetEarthModel", set_earth_model_));
            archive(cereal::make_nvp("FirstPoint", first_point_));
            archive(cereal::make_nvp("LastPoint", last_point_));
            archive(cereal::make_nvp("Distance", distance_));
            archive(cereal::make_nvp("SetPoints", set_points_));
            archive(cereal::make_nvp("ColumnDepthCache", column_depth_cache_));
            archive(cereal::make_nvp("Intersections", intersections_));
            archive(cereal::make_nvp("SetIntersections", set_intersections_));
        } else {
            throw std::runtime_error("Path only supports version <= 0!");
        }
    }

    bool HasEarthModel();
    bool HasPoints();
    bool HasIntersections();
    bool HasColumnDepth();
    bool HasTargetColumnDepth(std::set<LeptonInjector::Particle::ParticleType> const & targets);

    std::shared_ptr<const EarthModel> GetEarthModel();
    Vector3D GetFirstPoint();
    Vector3D GetLastPoint();
    Vector3D GetDirection();
    double GetDistance();
    Geometry::IntersectionList GetIntersections();

    void SetEarthModel(std::shared_ptr<const EarthModel> earth_model);
    void EnsureEarthModel();

    void SetPoints(Vector3D first_point, Vector3D last_point);
    void SetPointsWithRay(Vector3D first_point, Vector3D direction, double distance);
    void EnsurePoints();

    void SetIntersections(Geometry::IntersectionList const & intersections);
    void ComputeIntersections();
    void EnsureIntersections();

    void ClipToOuterBounds();

    void Flip();
    void ExtendFromEndByDistance(double distance);
    void ExtendFromEndByColumnDepth(double column_depth);
    void ExtendFromEndByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ExtendFromStartByDistance(double distance);
    void ExtendFromStartByColumnDepth(double column_depth);
    void ExtendFromStartByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ShrinkFromEndByDistance(double distance);
    void ShrinkFromEndByColumnDepth(double column_depth);
    void ShrinkFromEndByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ShrinkFromStartByDistance(double distance);
    void ShrinkFromStartByColumnDepth(double column_depth);
    void ShrinkFromStartByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);

    void ExtendFromEndToDistance(double distance);
    void ExtendFromEndToColumnDepth(double column_depth);
    void ExtendFromEndToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ExtendFromStartToDistance(double distance);
    void ExtendFromStartToColumnDepth(double column_depth);
    void ExtendFromStartToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ShrinkFromEndToDistance(double distance);
    void ShrinkFromEndToColumnDepth(double column_depth);
    void ShrinkFromEndToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    void ShrinkFromStartToDistance(double distance);
    void ShrinkFromStartToColumnDepth(double column_depth);
    void ShrinkFromStartToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);

    double GetColumnDepthInBounds();
    double GetColumnDepthInBounds(std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromStartInBounds(double distance);
    double GetColumnDepthFromStartInBounds(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromEndInBounds(double distance);
    double GetColumnDepthFromEndInBounds(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromStartAlongPath(double distance);
    double GetColumnDepthFromStartAlongPath(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromEndAlongPath(double distance);
    double GetColumnDepthFromEndAlongPath(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromStartInReverse(double distance);
    double GetColumnDepthFromStartInReverse(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetColumnDepthFromEndInReverse(double distance);
    double GetColumnDepthFromEndInReverse(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets);

    double GetDistanceFromStartInBounds(double column_depth);
    double GetDistanceFromStartInBounds(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetDistanceFromEndInBounds(double column_depth);
    double GetDistanceFromEndInBounds(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetDistanceFromStartAlongPath(double column_depth);
    double GetDistanceFromStartAlongPath(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetDistanceFromEndAlongPath(double column_depth);
    double GetDistanceFromEndAlongPath(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetDistanceFromStartInReverse(double column_depth);
    double GetDistanceFromStartInReverse(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);
    double GetDistanceFromEndInReverse(double column_depth);
    double GetDistanceFromEndInReverse(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets);

    bool IsWithinBounds(Vector3D point);
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::Path, 0);

# endif // LI_Path_H

