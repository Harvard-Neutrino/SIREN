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

#include "LeptonInjector/serialization/array.h"

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/detector/DensityDistribution.h"
#include "LeptonInjector/detector/MaterialModel.h"

#include "LeptonInjector/utilities/Constants.h"
#include "LeptonInjector/dataclasses/Particle.h"

namespace LI {
namespace detector {

struct EarthSector {
    std::string name;
    int material_id;
    int level;
    std::shared_ptr<const geometry::Geometry> geo;
    std::shared_ptr<const DensityDistribution> density;
    bool operator==(EarthSector const & o) const;
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
    math::Vector3D detector_origin_;
public:
    EarthModel();
    EarthModel(std::string const & earth_model, std::string const & material_model);
    EarthModel(std::string const & path, std::string const & earth_model, std::string const & material_model);

    bool operator==(EarthModel const & o) const;

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

    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0) const;
    double GetMassDensity(math::Vector3D const & p0) const;
    double GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, LI::dataclasses::Particle::ParticleType target) const;
    double GetParticleDensity(math::Vector3D const & p0, LI::dataclasses::Particle::ParticleType target) const;
    double GetInteractionDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double GetInteractionDensity(math::Vector3D const & p0,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;

    double GetColumnDepthInCGS(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, math::Vector3D const & p1) const;
    double GetColumnDepthInCGS(math::Vector3D const & p0, math::Vector3D const & p1) const;
    double DistanceForColumnDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthFromPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double column_depth) const;

    // Density/CD calculations with general target list, not just nucleon/electron
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, std::set<LI::dataclasses::Particle::ParticleType> targets) const;
    double GetMassDensity(math::Vector3D const & p0, std::set<LI::dataclasses::Particle::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LI::dataclasses::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LI::dataclasses::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(math::Vector3D const & p0, Iterator begin, Iterator end) const;
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, std::set<LI::dataclasses::Particle::ParticleType> targets) const;
    std::vector<double> GetParticleDensity(math::Vector3D const & p0, std::set<LI::dataclasses::Particle::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LI::dataclasses::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<LI::dataclasses::Particle::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(math::Vector3D const & p0, Iterator begin, Iterator end) const;

    double GetInteractionDepthInCGS(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, math::Vector3D const & p1,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double GetInteractionDepthInCGS(math::Vector3D const & p0, math::Vector3D const & p1,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthFromPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthToPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;
    double DistanceForInteractionDepthToPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections) const;

    std::vector<double> GetParticleColumnDepth(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, math::Vector3D const & p1, std::vector<LI::dataclasses::Particle::ParticleType> const & targets) const;

    EarthSector GetContainingSector(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0) const;
    EarthSector GetContainingSector(math::Vector3D const & p0) const;
    math::Vector3D GetEarthCoordPosFromDetCoordPos(math::Vector3D const & point) const;
    math::Vector3D GetEarthCoordDirFromDetCoordDir(math::Vector3D const & direction) const;
    math::Vector3D GetDetCoordPosFromEarthCoordPos(math::Vector3D const & point) const;
    math::Vector3D GetDetCoordDirFromEarthCoordDir(math::Vector3D const & direction) const;

    std::string GetPath() const;
    void SetPath(std::string const & path);

    MaterialModel const & GetMaterials() const;
    void SetMaterials(MaterialModel const & materials);

    std::vector<EarthSector> const & GetSectors() const;
    void SetSectors(std::vector<EarthSector> const & sectors);

    math::Vector3D GetDetectorOrigin() const;
    void SetDetectorOrigin(math::Vector3D const & detector_origin);

    void AddSector(EarthSector sector);
    EarthSector GetSector(int level) const;

    void ClearSectors();

    geometry::Geometry::IntersectionList GetIntersections(math::Vector3D const & p0, math::Vector3D const & direction) const;
    static void SortIntersections(geometry::Geometry::IntersectionList & intersections);
    static void SortIntersections(std::vector<geometry::Geometry::Intersection> & intersections);
    static void SectorLoop(std::function<bool(std::vector<geometry::Geometry::Intersection>::const_iterator, std::vector<geometry::Geometry::Intersection>::const_iterator, double)> callback, geometry::Geometry::IntersectionList const & intersections, bool reverse=false);

    static geometry::Geometry::IntersectionList GetOuterBounds(geometry::Geometry::IntersectionList const & intersections);
    geometry::Geometry::IntersectionList GetOuterBounds(math::Vector3D const & p0, math::Vector3D const & direction) const;
    std::set<LI::dataclasses::Particle::ParticleType> GetAvailableTargets(std::array<double,3> const & vertex) const;
    std::set<LI::dataclasses::Particle::ParticleType> GetAvailableTargets(geometry::Geometry::IntersectionList const & intersections, std::array<double,3> const & vertex) const;

private:
    void LoadDefaultMaterials();
    void LoadDefaultSectors();
public:
    void LoadConcentricShellsFromLegacyFile(std::string fname, double detector_depth, double ice_angle=-1);

    double GetTargetMass(LI::dataclasses::Particle::ParticleType target) const;
};

} // namespace detector
} // namespace LI

#include "LeptonInjector/detector/EarthModel.tcc"

CEREAL_CLASS_VERSION(LI::detector::EarthModel, 0);

#endif // LI_EarthModel_H
