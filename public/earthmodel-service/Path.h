#ifndef LI_Path_H
#define LI_Path_H

#include <memory>
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

    double column_depth_nucleon_ = 0;
    double column_depth_electron_ = 0;
    bool set_column_depth_nucleon_ = false;
    bool set_column_depth_electron_ = false;

    Geometry::IntersectionList intersections_;
    bool set_intersections_ = false;
public:
    Path();
    Path(std::shared_ptr<const EarthModel> earth_model);
    Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & last_point);
    Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & direction, double distance);

    bool HasEarthModel();
    bool HasPoints();
    bool HasIntersections();
    bool HasNucleonColumnDepth();
    bool HasElectronColumnDepth();

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
    void ExtendFromEndByColumnDepth(double column_depth, bool use_electron_density=false);
    void ExtendFromStartByDistance(double distance);
    void ExtendFromStartByColumnDepth(double column_depth, bool use_electron_density=false);
    void ShrinkFromEndByDistance(double distance);
    void ShrinkFromEndByColumnDepth(double column_depth, bool use_electron_density=false);
    void ShrinkFromStartByDistance(double distance);
    void ShrinkFromStartByColumnDepth(double column_depth, bool use_electron_density=false);

    void ExtendFromEndToDistance(double distance);
    void ExtendFromEndToColumnDepth(double column_depth, bool use_electron_density=false);
    void ExtendFromStartToDistance(double distance);
    void ExtendFromStartToColumnDepth(double column_depth, bool use_electron_density=false);
    void ShrinkFromEndToDistance(double distance);
    void ShrinkFromEndToColumnDepth(double column_depth, bool use_electron_density=false);
    void ShrinkFromStartToDistance(double distance);
    void ShrinkFromStartToColumnDepth(double column_depth, bool use_electron_density=false);

    double GetColumnDepthInBounds(bool use_electron_density=false);
    double GetColumnDepthFromStartInBounds(double distance, bool use_electron_density=false);
    double GetColumnDepthFromEndInBounds(double distance, bool use_electron_density=false);
    double GetColumnDepthFromStartAlongPath(double distance, bool use_electron_density=false);
    double GetColumnDepthFromEndAlongPath(double distance, bool use_electron_density=false);
    double GetColumnDepthFromStartInReverse(double distance, bool use_electron_density=false);
    double GetColumnDepthFromEndInReverse(double distance, bool use_electron_density=false);

    double GetDistanceFromStartInBounds(double column_depth, bool use_electron_density=false);
    double GetDistanceFromEndInBounds(double column_depth, bool use_electron_density=false);
    double GetDistanceFromStartAlongPath(double column_depth, bool use_electron_density=false);
    double GetDistanceFromEndAlongPath(double column_depth, bool use_electron_density=false);
    double GetDistanceFromStartInReverse(double column_depth, bool use_electron_density=false);
    double GetDistanceFromEndInReverse(double column_depth, bool use_electron_density=false);
};

} // namespace earthmodel

# endif // LI_Path_H

