#pragma once
#ifndef LI_EarthModel_H
#define LI_EarthModel_H

#include <memory>
#include <string>
#include <vector>

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

#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/Placement.h"
#include "earthmodel-service/DensityDist.h"
#include "earthmodel-service/MaterialModel.h"

#include "LeptonInjector/Constants.h"
#include "LeptonInjector/Particle.h"

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
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Name", name));
            archive(cereal::make_nvp("MaterialID", material_id));
            archive(cereal::make_nvp("Level", level));
            archive(cereal::make_nvp("Geometry", geo));
            archive(cereal::make_nvp("density", density));
        } else {
            throw std::runtime_error("EarthSector only supports version <= 0!");
        }
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

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Path", path_));
            archive(cereal::make_nvp("MaterialModel", materials_));
            archive(cereal::make_nvp("Sectors", sectors_));
            archive(cereal::make_nvp("SectorMap", sector_map_));
            archive(cereal::make_nvp("DetectorOrigin", detector_origin_));
        } else {
            throw std::runtime_error("EarthModel only supports version <= 0!");
        }
    }

    void LoadEarthModel(std::string const & earth_model);
    void LoadMaterialModel(std::string const & material_model);

    double GetMassDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0) const;
    double GetMassDensity(Vector3D const & p0) const;
    double GetParticleDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, LeptonInjector::Particle::ParticleType target) const;
    double GetParticleDensity(Vector3D const & p0, LeptonInjector::Particle::ParticleType target) const;

    double GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1) const;
    double GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1) const;
    double DistanceForColumnDepthFromPoint(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthFromPoint(Vector3D const & end_point, Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(Vector3D const & end_point, Vector3D const & direction, double column_depth) const;

    // Density/CD calculations with general target list, not just nucleon/electron
    double GetMassDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, std::set<LeptonInjector::Particle::ParticleType> targets) const;
    double GetMassDensity(Vector3D const & p0, std::set<LeptonInjector::Particle::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LeptonInjector::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LeptonInjector::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(Vector3D const & p0, Iterator begin, Iterator end) const;
    std::vector<double> GetParticleDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, std::set<LeptonInjector::Particle::ParticleType> targets) const;
    std::vector<double> GetParticleDensity(Vector3D const & p0, std::set<LeptonInjector::Particle::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LeptonInjector::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LeptonInjector::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(Vector3D const & p0, Iterator begin, Iterator end) const;

    double GetInteractionDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double GetInteractionDepthInCGS(Vector3D const & p0, Vector3D const & p1,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthFromPoint(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthFromPoint(Vector3D const & end_point, Vector3D const & direction, double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & end_point, Vector3D const & direction, double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthToPoint(Vector3D const & end_point, Vector3D const & direction, double interaction_depth,
            std::vector<LeptonInjector::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;

    std::vector<double> GetParticleColumnDepth(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1, std::vector<LeptonInjector::Particle::ParticleType> const & targets) const;

    EarthSector GetContainingSector(Geometry::IntersectionList const & intersections, Vector3D const & p0) const;
    EarthSector GetContainingSector(Vector3D const & p0) const;
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
    static void SectorLoop(std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback, Geometry::IntersectionList const & intersections, bool reverse=false);

    static Geometry::IntersectionList GetOuterBounds(Geometry::IntersectionList const & intersections);
    Geometry::IntersectionList GetOuterBounds(Vector3D const & p0, Vector3D const & direction);
    std::set<LeptonInjector::Particle::ParticleType> GetAvailableTargets(std::array<double,3> const & vertex);

private:
    void LoadDefaultMaterials();
    void LoadDefaultSectors();
public:
    void LoadConcentricShellsFromLegacyFile(std::string fname, double detector_depth, double ice_angle=-1);

    double GetTargetMass(LeptonInjector::Particle::ParticleType target) const;
};

}; // namespace earthmodel

#include "earthmodel-service/EarthModel.tcc"

CEREAL_CLASS_VERSION(earthmodel::EarthModel, 0);

#endif // LI_EarthModel_H
