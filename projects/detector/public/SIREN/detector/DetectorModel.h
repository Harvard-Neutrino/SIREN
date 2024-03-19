#pragma once
#ifndef SIREN_DetectorModel_H
#define SIREN_DetectorModel_H

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

#include "SIREN/serialization/array.h"

#include "SIREN/dataclasses/Particle.h"    // for Particle
#include "SIREN/detector/MaterialModel.h"  // for MaterialModel
#include "SIREN/geometry/Geometry.h"       // for Geometry
#include "SIREN/math/Vector3D.h"           // for Vector3D
#include "SIREN/math/Quaternion.h"         // for Quaternion
#include "SIREN/detector/Coordinates.h"

#include "SIREN/dataclasses/Particle.h"

namespace siren { namespace detector { class Path; } }
namespace siren { namespace detector { class DensityDistribution; } }

namespace siren {
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
friend siren::detector::Path;
    std::string path_;
    MaterialModel materials_;
    std::vector<DetectorSector> sectors_;
    std::map<int, unsigned int> sector_map_;
    math::Vector3D detector_origin_;
    math::Quaternion detector_rotation_;
public:
    DetectorModel();
    DetectorModel(std::string const & detector_model, std::string const & material_model);
    DetectorModel(std::string const & path, std::string const & detector_model, std::string const & material_model);

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

    void LoadDetectorModel(std::string const & detector_model);
    void LoadMaterialModel(std::string const & material_model);
private:
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0) const;
    double GetMassDensity(GeometryPosition const & p0) const;
    double GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, siren::dataclasses::ParticleType target) const;
    double GetParticleDensity(GeometryPosition const & p0, siren::dataclasses::ParticleType target) const;
    double GetInteractionDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDensity(GeometryPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

    double GetColumnDepthInCGS(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1) const;
    double GetColumnDepthInCGS(GeometryPosition const & p0, GeometryPosition const & p1) const;
    double DistanceForColumnDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & end_point, GeometryDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthFromPoint(GeometryPosition const & end_point, GeometryDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & end_point, GeometryDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(GeometryPosition const & end_point, GeometryDirection const & direction, double column_depth) const;

    // Density/CD calculations with general target list, not just nucleon/electron
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    double GetMassDensity(GeometryPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(GeometryPosition const & p0, Iterator begin, Iterator end) const;
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    std::vector<double> GetParticleDensity(GeometryPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(GeometryPosition const & p0, Iterator begin, Iterator end) const;

    double GetInteractionDepthInCGS(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDepthInCGS(GeometryPosition const & p0, GeometryPosition const & p1,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & end_point, GeometryDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(GeometryPosition const & end_point, GeometryDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & end_point, GeometryDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(GeometryPosition const & end_point, GeometryDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

    std::vector<double> GetParticleColumnDepth(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1, std::vector<siren::dataclasses::ParticleType> const & targets) const;

    DetectorSector GetContainingSector(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & p0) const;
    DetectorSector GetContainingSector(GeometryPosition const & p0) const;

    geometry::Geometry::IntersectionList GetIntersections(GeometryPosition const & p0, GeometryDirection const & direction) const;
    geometry::Geometry::IntersectionList GetOuterBounds(GeometryPosition const & p0, GeometryDirection const & direction) const;

    std::set<siren::dataclasses::ParticleType> GetAvailableTargets(GeometryPosition const & vertex) const;
    std::set<siren::dataclasses::ParticleType> GetAvailableTargets(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & vertex) const;

public:
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0) const;
    double GetMassDensity(DetectorPosition const & p0) const;
    double GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, siren::dataclasses::ParticleType target) const;
    double GetParticleDensity(DetectorPosition const & p0, siren::dataclasses::ParticleType target) const;
    double GetInteractionDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDensity(DetectorPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

    double GetColumnDepthInCGS(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1) const;
    double GetColumnDepthInCGS(DetectorPosition const & p0, DetectorPosition const & p1) const;
    double DistanceForColumnDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & end_point, DetectorDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthFromPoint(DetectorPosition const & end_point, DetectorDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & end_point, DetectorDirection const & direction, double column_depth) const;
    double DistanceForColumnDepthToPoint(DetectorPosition const & end_point, DetectorDirection const & direction, double column_depth) const;

    // Density/CD calculations with general target list, not just nucleon/electron
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    double GetMassDensity(DetectorPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    double GetMassDensity(DetectorPosition const & p0, Iterator begin, Iterator end) const;
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    std::vector<double> GetParticleDensity(DetectorPosition const & p0, std::set<siren::dataclasses::ParticleType> targets) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, Iterator begin, Iterator end) const;
    template<typename Iterator, typename = typename std::enable_if<std::is_same<siren::dataclasses::ParticleType, typename Iterator::value_type>::value, Iterator>::type>
    std::vector<double> GetParticleDensity(DetectorPosition const & p0, Iterator begin, Iterator end) const;

    double GetInteractionDepthInCGS(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double GetInteractionDepthInCGS(DetectorPosition const & p0, DetectorPosition const & p1,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & end_point, DetectorDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthFromPoint(DetectorPosition const & end_point, DetectorDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & end_point, DetectorDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;
    double DistanceForInteractionDepthToPoint(DetectorPosition const & end_point, DetectorDirection const & direction, double interaction_depth,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const;

    std::vector<double> GetParticleColumnDepth(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1, std::vector<siren::dataclasses::ParticleType> const & targets) const;

    DetectorSector GetContainingSector(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & p0) const;
    DetectorSector GetContainingSector(DetectorPosition const & p0) const;

    geometry::Geometry::IntersectionList GetIntersections(DetectorPosition const & p0, DetectorDirection const & direction) const;
    geometry::Geometry::IntersectionList GetOuterBounds(DetectorPosition const & p0, DetectorDirection const & direction) const;

    std::set<siren::dataclasses::ParticleType> GetAvailableTargets(DetectorPosition const & vertex) const;
    std::set<siren::dataclasses::ParticleType> GetAvailableTargets(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & vertex) const;

    DetectorPosition ToDet(GeometryPosition const &) const;
    DetectorDirection ToDet(GeometryDirection const &) const;
    DetectorPosition ToDet(GeometryPosition &&) const;
    DetectorDirection ToDet(GeometryDirection &&) const;
    GeometryPosition ToGeo(DetectorPosition const &) const;
    GeometryDirection ToGeo(DetectorDirection const &) const;
    GeometryPosition ToGeo(DetectorPosition &&) const;
    GeometryDirection ToGeo(DetectorDirection &&) const;


    std::string GetPath() const;
    void SetPath(std::string const & path);

    MaterialModel const & GetMaterials() const;
    void SetMaterials(MaterialModel const & materials);

    std::vector<DetectorSector> const & GetSectors() const;
    void SetSectors(std::vector<DetectorSector> const & sectors);

    GeometryPosition GetDetectorOrigin() const;
    void SetDetectorOrigin(GeometryPosition const & detector_origin);
    siren::math::Quaternion GetDetectorRotation() const;
    void SetDetectorRotation(siren::math::Quaternion const & detector_rotation);

    void AddSector(DetectorSector sector);
    DetectorSector GetSector(int level) const;

    void ClearSectors();

    static void SortIntersections(geometry::Geometry::IntersectionList & intersections);
    static void SortIntersections(std::vector<geometry::Geometry::Intersection> & intersections);
    static void SectorLoop(std::function<bool(std::vector<geometry::Geometry::Intersection>::const_iterator, std::vector<geometry::Geometry::Intersection>::const_iterator, double)> callback, geometry::Geometry::IntersectionList const & intersections, bool reverse=false);

    static geometry::Geometry::IntersectionList GetOuterBounds(geometry::Geometry::IntersectionList const & intersections);

    static std::tuple<siren::math::Vector3D, siren::math::Quaternion> ParseDetector(std::stringstream & ss);
    static std::shared_ptr<geometry::Geometry> ParseGeometryObject(std::stringstream & line);
    static int ParseMaterialID(std::stringstream & line, MaterialModel const & materials);
    static std::shared_ptr<detector::DensityDistribution> ParseDensityDistribution(std::stringstream & line);
    static std::shared_ptr<geometry::Geometry> ParseFiducialVolume(std::string fiducial_line, std::string origin_line);
    static std::shared_ptr<geometry::Geometry> ParseFiducialVolume(std::string line, siren::math::Vector3D detector_origin, siren::math::Quaternion detector_quaternion);

private:
    void LoadDefaultMaterials();
    void LoadDefaultSectors();
public:
    void LoadConcentricShellsFromLegacyFile(std::string fname, double detector_depth, double ice_angle=-1);

    double GetTargetMass(siren::dataclasses::ParticleType target) const;
};

} // namespace detector
} // namespace siren


std::ostream& operator<<(std::ostream& oss, siren::detector::DetectorSector const & bcm);
std::ostream& operator<<(std::ostream& oss, siren::detector::DetectorSector & bcm);

#include "SIREN/detector/DetectorModel.tcc"

CEREAL_CLASS_VERSION(siren::detector::DetectorModel, 0);

#endif // SIREN_DetectorModel_H
