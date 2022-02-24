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
            archive(cereal::make_nvp("Intersections", intersections_));
        } else {
            throw std::runtime_error("Path only supports version <= 0!");
        }
    }

    bool HasEarthModel();
    bool HasPoints();
    bool HasIntersections();
    bool HasColumnDepth();

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

    //TODO Clean up cxx and then EarthModel

    // Extend/Shrink By
    void ExtendFromEndByDistance(double distance);
    void ExtendFromStartByDistance(double distance);
    void ShrinkFromEndByDistance(double distance);
    void ShrinkFromStartByDistance(double distance);

    void ExtendFromEndByColumnDepth(double column_depth);
    void ExtendFromStartByColumnDepth(double column_depth);
    void ShrinkFromEndByColumnDepth(double column_depth);
    void ShrinkFromStartByColumnDepth(double column_depth);

    void ExtendFromEndByInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ExtendFromStartByInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ShrinkFromEndByInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ShrinkFromStartByInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);

    // Extend/Shrink To
    void ExtendFromEndToDistance(double distance);
    void ExtendFromStartToDistance(double distance);
    void ShrinkFromEndToDistance(double distance);
    void ShrinkFromStartToDistance(double distance);

    void ExtendFromEndToColumnDepth(double column_depth);
    void ExtendFromStartToColumnDepth(double column_depth);
    void ShrinkFromEndToColumnDepth(double column_depth);
    void ShrinkFromStartToColumnDepth(double column_depth);

    void ExtendFromEndToInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ExtendFromStartToInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ShrinkFromEndToInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    void ShrinkFromStartToInteractionDepth(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    //

    // Get
    double GetColumnDepthInBounds();
    double GetInteractionDepthInBounds(
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    //

    // Get * From
    double GetColumnDepthFromStartInBounds(double distance);
    double GetColumnDepthFromEndInBounds(double distance);
    double GetColumnDepthFromStartAlongPath(double distance);
    double GetColumnDepthFromEndAlongPath(double distance);
    double GetColumnDepthFromStartInReverse(double distance);
    double GetColumnDepthFromEndInReverse(double distance);

    double GetInteractionDepthFromStartInBounds(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetInteractionDepthFromEndInBounds(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetInteractionDepthFromStartAlongPath(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetInteractionDepthFromEndAlongPath(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetInteractionDepthFromStartInReverse(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetInteractionDepthFromEndInReverse(double distance,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    //

    // Get Distance From
    double GetDistanceFromStartInBounds(double column_depth);
    double GetDistanceFromEndInBounds(double column_depth);
    double GetDistanceFromStartAlongPath(double column_depth);
    double GetDistanceFromEndAlongPath(double column_depth);
    double GetDistanceFromStartInReverse(double column_depth);
    double GetDistanceFromEndInReverse(double column_depth);

    double GetDistanceFromStartInBounds(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetDistanceFromEndInBounds(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetDistanceFromStartAlongPath(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetDistanceFromEndAlongPath(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetDistanceFromStartInReverse(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    double GetDistanceFromEndInReverse(double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections);
    //

    bool IsWithinBounds(Vector3D point);
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::Path, 0);

# endif // LI_Path_H

