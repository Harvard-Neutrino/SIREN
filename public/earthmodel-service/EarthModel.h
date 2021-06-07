#ifndef LI_EarthModel_H
#define LI_EarthModel_H

#include <memory>
#include <string>
#include <vector>

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/DensityDist.h"
#include "earthmodel-service/MaterialModel.h"

#include "LeptonInjector/Constants.h"

namespace earthmodel {

struct EarthSector {
    std::string name;
    int material_id;
    int level;
    std::shared_ptr<const Geometry> geo;
    std::shared_ptr<const DensityDistribution> density;
    bool operator==(EarthSector const & o) const {
        return name == o.name and material_id == o.material_id and level == o.level and geo == o.geo and density == o.density;
    }
};

class EarthModel {
private:
    std::string path_;
    MaterialModel materials_;
    std::vector<EarthSector> sectors_;
    std::map<int, unsigned int> sector_map_;
    Vector3D detector_origin_;
public:
    EarthModel();
    EarthModel(std::string const & earth_model, std::string const & material_model);
    EarthModel(std::string const & path, std::string const & earth_model, std::string const & material_model);

    void LoadEarthModel(std::string const & earth_model);
    void LoadMaterialModel(std::string const & material_model);

    double GetDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, bool use_electron_density=false) const;
    double GetDensity(Vector3D const & p0, bool use_electron_density=false) const;
    double GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1, bool use_electron_density=false) const;
    double GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1, bool use_electron_density=false) const;
    double DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double column_depth, bool use_electron_density=false) const;
    double DistanceForColumnDepthToPoint(Vector3D const & end_point, Vector3D const & direction, double column_depth, bool use_electron_density=false) const;
    Vector3D GetEarthCoordPosFromDetCoordPos(Vector3D const & point) const;
    Vector3D GetEarthCoordDirFromDetCoordDir(Vector3D const & direction) const;
    Vector3D GetDetCoordPosFromEarthCoordPos(Vector3D const & point) const;
    Vector3D GetDetCoordDirFromEarthCoordDir(Vector3D const & direction) const;

    std::string GetPath() const;
    void SetPath(std::string const & path);

    MaterialModel const & GetMaterials() const;
    void SetMaterials(MaterialModel const & materials);

    std::vector<EarthSector> const & GetSectors() const;
    void SetSectors(std::vector<EarthSector> const & sectors);

    Vector3D GetDetectorOrigin() const;
    void SetDetectorOrigin(Vector3D const & detector_origin);

    void AddSector(EarthSector sector);
    EarthSector GetSector(int level) const;

    void ClearSectors();

    Geometry::IntersectionList GetIntersections(Vector3D const & p0, Vector3D const & direction) const;
    static void SortIntersections(Geometry::IntersectionList & intersections);
    static void SortIntersections(std::vector<Geometry::Intersection> & intersections);
    void SectorLoop(std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback, Geometry::IntersectionList const & intersections) const;

private:
    void LoadDefaultMaterials();
    void LoadDefaultSectors();
public:
    void LoadConcentricShellsFromLegacyFile(std::string fname, double detector_depth, double ice_angle=-1);
};

}; // namespace earthmodel

#endif // LI_EarthModel_H
