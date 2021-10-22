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

    std::map<std::vector<LeptonInjector::Particle::ParticleType>, double> column_depth_cache_;

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
    bool HasTargetColumnDepth(std::vector<LeptonInjector::Particle::ParticleType> const & targets);

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
    void ExtendFromEndByColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ExtendFromStartByDistance(double distance);
    void ExtendFromStartByColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ShrinkFromEndByDistance(double distance);
    void ShrinkFromEndByColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ShrinkFromStartByDistance(double distance);
    void ShrinkFromStartByColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});

    void ExtendFromEndToDistance(double distance);
    void ExtendFromEndToColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ExtendFromStartToDistance(double distance);
    void ExtendFromStartToColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ShrinkFromEndToDistance(double distance);
    void ShrinkFromEndToColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    void ShrinkFromStartToDistance(double distance);
    void ShrinkFromStartToColumnDepth(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});

    double GetColumnDepthInBounds(std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromStartInBounds(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromEndInBounds(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromStartAlongPath(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromEndAlongPath(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromStartInReverse(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetColumnDepthFromEndInReverse(double distance, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});

    double GetDistanceFromStartInBounds(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetDistanceFromEndInBounds(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetDistanceFromStartAlongPath(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetDistanceFromEndAlongPath(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetDistanceFromStartInReverse(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
    double GetDistanceFromEndInReverse(double column_depth, std::vector<LeptonInjector::Particle::ParticleType> const & targets=std::vector<LeptonInjector::Particle::ParticleType>{LeptonInjector::Particle::ParticleType::Nucleon});
};

} // namespace earthmodel

CEREAL_CLASS_VERSION(earthmodel::Path, 0);

# endif // LI_Path_H

