#pragma once
#ifndef LI_Path_H
#define LI_Path_H

#include <vector>                                 // for vector
#include <memory>                                 // for shared_ptr
#include <cstdint>                                // for uint32_t
#include <stdexcept>                              // for runtime_error

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

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/dataclasses/Particle.h"  // for Particle
#include "LeptonInjector/geometry/Geometry.h"     // for Geometry
#include "LeptonInjector/detector/DetectorModel.h"   // for DetectorModel
#include "LeptonInjector/math/Vector3D.h"         // for Vector3D

namespace LI {
namespace detector {

class Path {
private:
    std::shared_ptr<const DetectorModel> earth_model_;
    bool set_earth_model_ = false;

    math::Vector3D first_point_;
    math::Vector3D last_point_;
    math::Vector3D direction_;
    double distance_ = 0;
    bool set_points_ = false;

    double column_depth_cached_;
    bool set_column_depth_ = false;

    geometry::Geometry::IntersectionList intersections_;
    bool set_intersections_ = false;
public:
    Path();
    Path(std::shared_ptr<const DetectorModel> earth_model);
    Path(std::shared_ptr<const DetectorModel> earth_model, math::Vector3D const & first_point, math::Vector3D const & last_point);
    Path(std::shared_ptr<const DetectorModel> earth_model, math::Vector3D const & first_point, math::Vector3D const & direction, double distance);

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("DetectorModel", earth_model_));
            archive(cereal::make_nvp("SetDetectorModel", set_earth_model_));
            archive(cereal::make_nvp("FirstPoint", first_point_));
            archive(cereal::make_nvp("LastPoint", last_point_));
            archive(cereal::make_nvp("Distance", distance_));
            archive(cereal::make_nvp("Intersections", intersections_));
        } else {
            throw std::runtime_error("Path only supports version <= 0!");
        }
    }

    bool HasDetectorModel();
    bool HasPoints();
    bool HasIntersections();
    bool HasColumnDepth();

    std::shared_ptr<const DetectorModel> GetDetectorModel();
    math::Vector3D const & GetFirstPoint();
    math::Vector3D const & GetLastPoint();
    math::Vector3D const & GetDirection();
    double GetDistance();
    geometry::Geometry::IntersectionList const & GetIntersections();

    void SetDetectorModel(std::shared_ptr<const DetectorModel> earth_model);
    void EnsureDetectorModel();

    void SetPoints(math::Vector3D first_point, math::Vector3D last_point);
    void SetPointsWithRay(math::Vector3D first_point, math::Vector3D direction, double distance);
    void EnsurePoints();

    void SetIntersections(geometry::Geometry::IntersectionList const & intersections);
    void ComputeIntersections();
    void EnsureIntersections();

    void ClipToOuterBounds();

    void Flip();

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
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ExtendFromStartByInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ShrinkFromEndByInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ShrinkFromStartByInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);

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
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ExtendFromStartToInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ShrinkFromEndToInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    void ShrinkFromStartToInteractionDepth(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    //

    // Get
    double GetColumnDepthInBounds();
    double GetInteractionDepthInBounds(
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    //

    // Get * From
    double GetColumnDepthFromStartInBounds(double distance);
    double GetColumnDepthFromEndInBounds(double distance);
    double GetColumnDepthFromStartAlongPath(double distance);
    double GetColumnDepthFromEndAlongPath(double distance);
    double GetColumnDepthFromStartInReverse(double distance);
    double GetColumnDepthFromEndInReverse(double distance);

    double GetInteractionDepthFromStartInBounds(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetInteractionDepthFromEndInBounds(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetInteractionDepthFromStartAlongPath(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetInteractionDepthFromEndAlongPath(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetInteractionDepthFromStartInReverse(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetInteractionDepthFromEndInReverse(double distance,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    //

    // Get Distance From
    double GetDistanceFromStartInBounds(double column_depth);
    double GetDistanceFromEndInBounds(double column_depth);
    double GetDistanceFromStartAlongPath(double column_depth);
    double GetDistanceFromEndAlongPath(double column_depth);
    double GetDistanceFromStartInReverse(double column_depth);
    double GetDistanceFromEndInReverse(double column_depth);

    double GetDistanceFromStartInBounds(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetDistanceFromEndInBounds(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetDistanceFromStartAlongPath(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetDistanceFromEndAlongPath(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetDistanceFromStartInReverse(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    double GetDistanceFromEndInReverse(double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length);
    //

    bool IsWithinBounds(math::Vector3D point);
    double GetDistanceFromStartInBounds(math::Vector3D point);
};

} // namespace detector
} // namespace LI

CEREAL_CLASS_VERSION(LI::detector::Path, 0);

# endif // LI_Path_H

