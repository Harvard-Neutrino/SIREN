#pragma once
#ifndef LI_DetectorModel_H
#define LI_DetectorModel_H

#include <map>                                      // for map
#include <set>                                      // for set
#include <array>                                    // for array
#include <memory>                                   // for shared_ptr
#include <cstdint>                                  // for uint32_t
#include <stdexcept>                                // for runtime_error
#include <functional>                               // for function
#include <type_traits>                              // for enable_if, is_same

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

#include "LeptonInjector/dataclasses/Particle.h"    // for Particle
#include "LeptonInjector/detector/MaterialModel.h"  // for MaterialModel
#include "LeptonInjector/geometry/Geometry.h"       // for Geometry
#include "LeptonInjector/math/Vector3D.h"           // for Vector3D

#include "LeptonInjector/dataclasses/Particle.h"

namespace LI { namespace detector { class DensityDistribution; } }

namespace LI {
namespace detector {

struct DetectorSector {
    std::string name;
    int material_id;
    int level;
    std::shared_ptr<const geometry::Geometry> geo;
    std::shared_ptr<const DensityDistribution> density;
    bool operator==(DetectorSector const & o) const;
    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Name", name));
            archive(cereal::make_nvp("MaterialID", material_id));
            archive(cereal::make_nvp("Level", level));
            archive(cereal::make_nvp("Geometry", geo));
            archive(cereal::make_nvp("density", density));
        } else {
            throw std::runtime_error("DetectorSector only supports version <= 0!");
        }
    }
	std::ostream & Print(std::ostream& oss) const;
};

class DetectorModel {
private:
    std::string path_;
    MaterialModel materials_;
    std::vector<DetectorSector> sectors_;
    std::map<int, unsigned int> sector_map_;
    math::Vector3D detector_origin_;
public:
    DetectorModel();
    DetectorModel(std::string const & earth_model, std::string const & material_model);
    DetectorModel(std::string const & path, std::string const & earth_model, std::string const & material_model);

    bool operator==(DetectorModel const & o) const;

    template<class Archive>
    void serialize(Archive & archive, std::uint32_t const version) {
        if(version == 0) {
            archive(cereal::make_nvp("Path", path_));
            archive(cereal::make_nvp("MaterialModel", materials_));
            archive(cereal::make_nvp("Sectors", sectors_));
            archive(cereal::make_nvp("SectorMap", sector_map_));
            archive(cereal::make_nvp("DetectorOrigin", detector_origin_));
        } else {
            throw std::runtime_error("DetectorModel only supports version <= 0!");
        }
    }

    void LoadDetectorModel(std::string const & earth_model);
    void LoadMaterialModel(std::string const & material_model);

    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0) const;
    double GetMassDensity(math::Vector3D const & p0) const;
    double GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, LI::dataclasses::Particle::ParticleType target) const;
    double GetParticleDensity(math::Vector3D const & p0, LI::dataclasses::Particle::ParticleType target) const;
    double GetInteractionDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDensity(math::Vector3D const & p0,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

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
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDepthInCGS(math::Vector3D const & p0, math::Vector3D const & p1,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(math::Vector3D const & end_point, math::Vector3D const & direction, double interaction_depth,
            std::vector<LI::dataclasses::Particle::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

    std::vector<double> GetParticleColumnDepth(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, math::Vector3D const & p1, std::vector<LI::dataclasses::Particle::ParticleType> const & targets) const;

    DetectorSector GetContainingSector(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0) const;
    DetectorSector GetContainingSector(math::Vector3D const & p0) const;
    math::Vector3D GetEarthCoordPosFromDetCoordPos(math::Vector3D const & point) const;
    math::Vector3D GetEarthCoordDirFromDetCoordDir(math::Vector3D const & direction) const;
    math::Vector3D GetDetCoordPosFromEarthCoordPos(math::Vector3D const & point) const;
    math::Vector3D GetDetCoordDirFromEarthCoordDir(math::Vector3D const & direction) const;

    std::string GetPath() const;
    void SetPath(std::string const & path);

    MaterialModel const & GetMaterials() const;
    void SetMaterials(MaterialModel const & materials);

    std::vector<DetectorSector> const & GetSectors() const;
    void SetSectors(std::vector<DetectorSector> const & sectors);

    math::Vector3D GetDetectorOrigin() const;
    void SetDetectorOrigin(math::Vector3D const & detector_origin);

    void AddSector(DetectorSector sector);
    DetectorSector GetSector(int level) const;

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


std::ostream& operator<<(std::ostream& oss, LI::detector::DetectorSector const & bcm);
std::ostream& operator<<(std::ostream& oss, LI::detector::DetectorSector & bcm);

#include "LeptonInjector/detector/DetectorModel.tcc"

CEREAL_CLASS_VERSION(LI::detector::DetectorModel, 0);

#endif // LI_DetectorModel_H