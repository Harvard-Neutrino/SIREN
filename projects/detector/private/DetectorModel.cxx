#include "SIREN/detector/DetectorModel.h"

#include <tuple>
#include <cmath>
#include <cctype>
#include <memory>
#include <vector>
#include <string>
#include <limits>
#include <numeric>
#include <utility>
#include <assert.h>
#include <iostream>
#include <iterator>
#include <stddef.h>
#include <stdlib.h>
#include <algorithm>
#include <initializer_list>

#include "SIREN/math/Vector3D.h"
#include "SIREN/math/EulerQuaternionConversions.h"

#include "SIREN/detector/MaterialModel.h"

#include "SIREN/detector/RadialAxis1D.h"
#include "SIREN/detector/DensityDistribution.h"
#include "SIREN/detector/ConstantDensityDistribution.h"
#include "SIREN/detector/RadialAxisPolynomialDensityDistribution.h"

#include "SIREN/geometry/Box.h"
#include "SIREN/geometry/Sphere.h"
#include "SIREN/geometry/Cylinder.h"
#include "SIREN/geometry/ExtrPoly.h"

#include "SIREN/geometry/Placement.h"

#include "SIREN/utilities/Constants.h"

using namespace siren::math;
using namespace siren::geometry;
using namespace siren::detector;

namespace {
void string_to_lower(std::string & data) {
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c){ return std::tolower(c); });
}

template <class InIt>
typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
    typedef typename std::iterator_traits<InIt>::value_type real;
    real sum = real(0);
    real running_error = real(0);
    real temp;
    real difference;

    for (; begin != end; ++begin) {
        difference = *begin;
        difference -= running_error;
        temp = sum;
        temp += difference;
        running_error = temp;
        running_error -= sum;
        running_error -= difference;
        sum = std::move(temp);
    }
    return sum;
};

template<typename T>
T accumulate(std::initializer_list<T> list) {
    return accumulate(list.begin(), list.end());
}

template <typename T>
class accumulator {
private:
    T sum;
    T running_error;
    T temp;
    T difference;
public:
    accumulator() : sum(0), running_error(0) {};
    accumulator(std::initializer_list<T> const & list) : sum(0), running_error(0) {
        for(auto const & val : list) {
            accumulate(val);
        }
    }
    void accumulate(T const & val) {
        difference = val;
        difference -= running_error;
        temp = sum;
        temp += difference;
        running_error = temp;
        running_error -= sum;
        running_error -= difference;
        sum = std::move(temp);
    }
    T result() {
        return sum;
    }
    operator T() {
        return sum;
    }
};
}

bool DetectorSector::operator==(DetectorSector const & o) const {
    return name == o.name and material_id == o.material_id and level == o.level and geo == o.geo and density == o.density;
}

std::ostream & DetectorSector::Print(std::ostream& oss) const {
    oss << "[DetectorSector:\n"
        << "         Name : " << name << '\n'
        << "   MaterialID : " << material_id << '\n'
        << "        Level : " << level << '\n'
        << "          Geo : " << geo << '\n'
        << "      Density : " << density << "\n]";
    return oss;
}

DetectorPosition DetectorModel::ToDet(GeometryPosition const & pos) const {
    return DetectorPosition(detector_rotation_.rotate(pos - detector_origin_, true));
}

DetectorDirection DetectorModel::ToDet(GeometryDirection const & dir) const {
    return DetectorDirection(detector_rotation_.rotate(dir, true));
}

DetectorPosition DetectorModel::ToDet(GeometryPosition && pos) const {
    return DetectorPosition(detector_rotation_.rotate(pos - detector_origin_, true));
}

DetectorDirection DetectorModel::ToDet(GeometryDirection && dir) const {
    return DetectorDirection(detector_rotation_.rotate(dir, true));
}

GeometryPosition DetectorModel::ToGeo(DetectorPosition const & pos) const {
    return GeometryPosition(detector_rotation_.rotate(pos, false) + detector_origin_);
}

GeometryDirection DetectorModel::ToGeo(DetectorDirection const & dir) const {
    return GeometryDirection(detector_rotation_.rotate(dir, false));
}

GeometryPosition DetectorModel::ToGeo(DetectorPosition && pos) const {
    return GeometryPosition(detector_rotation_.rotate(pos, false) + detector_origin_);
}

GeometryDirection DetectorModel::ToGeo(DetectorDirection && dir) const {
    return GeometryDirection(detector_rotation_.rotate(dir, false));
}

std::ostream& operator<<(std::ostream& oss, DetectorSector const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, DetectorSector & bcm) {
    return(bcm.Print(oss));
}

DetectorModel::DetectorModel() {
    LoadDefaultMaterials();
    LoadDefaultSectors();
}

DetectorModel::DetectorModel(std::string const & detector_model, std::string const & material_model) {
    LoadDefaultMaterials();
    LoadDefaultSectors();
    LoadMaterialModel(material_model);
    LoadDetectorModel(detector_model);
}

DetectorModel::DetectorModel(std::string const & path, std::string const & detector_model, std::string const & material_model) : path_(path) {
    LoadDefaultMaterials();
    LoadDefaultSectors();
    LoadMaterialModel(material_model);
    LoadDetectorModel(detector_model);
}

bool DetectorModel::operator==(DetectorModel const & o) const {
    return
        std::tie(materials_, sectors_, sector_map_, detector_origin_)
        ==
        std::tie(o.materials_, o.sectors_, o.sector_map_, o.detector_origin_);
}

std::string DetectorModel::GetPath() const {
    return path_;
}

void DetectorModel::SetPath(std::string const & path) {
    path_ = path;
}

MaterialModel const & DetectorModel::GetMaterials() const {
    return materials_;
}

void DetectorModel::SetMaterials(MaterialModel const & materials) {
    materials_ = materials;
}

std::vector<DetectorSector> const & DetectorModel::GetSectors() const {
    return sectors_;
}

void DetectorModel::SetSectors(std::vector<DetectorSector> const & sectors) {
    sectors_ = sectors;
}

GeometryPosition DetectorModel::GetDetectorOrigin() const {
    return GeometryPosition(detector_origin_);
}

void DetectorModel::SetDetectorOrigin(GeometryPosition const & detector_origin) {
    detector_origin_ = detector_origin;
}

void DetectorModel::AddSector(DetectorSector sector) {
    if(sector_map_.count(sector.level) > 0) {
        throw(std::runtime_error("Already have a sector of that heirarchy!"));
    }
    else {
        sector_map_[sector.level] = sectors_.size();
        sectors_.push_back(sector);
    }
}

DetectorSector DetectorModel::GetSector(int heirarchy) const {
    auto const iter = sector_map_.find(heirarchy);
    assert(iter != sector_map_.end());
    unsigned int index = sector_map_.at(heirarchy);
    assert(index < sectors_.size());
    unsigned int alt_index = sector_map_.find(heirarchy)->second;
    assert(index == alt_index);
    return sectors_[index];
}

void DetectorModel::ClearSectors() {
    sectors_.clear();
    sector_map_.clear();
}

namespace {
bool fexists(const std::string filename)
{
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}
}

std::shared_ptr<siren::geometry::Geometry> DetectorModel::ParseGeometryObject(std::stringstream & ss) {
    std::string shape;
    double xc, yc, zc; // Coordinates of the center of the shape
    double alpha, beta, gamma; // Euler Angles of shape rotation

    ss >> shape;
    ss >> xc >> yc >> zc;
    ss >> alpha >> beta >> gamma;
    Placement placement(Vector3D(xc,yc,zc), QFromZXZr(alpha,beta,gamma));

    std::shared_ptr<geometry::Geometry> geo;

    if(shape.find("sphere")!=std::string::npos) {
        double radius; // For Sphere shapes
        ss >> radius;
        geo = Sphere(placement, radius, 0).create();
    }
    else if(shape.find("box")!=std::string::npos) {
        double dx, dy, dz; // For Box shapes
        ss >> dx >> dy >> dz;
        geo = Box(placement, dx, dy, dz).create();
    }
    else if(shape.find("cylinder")!=std::string::npos) {
        double _or, ir, z; // For Cylinder shapes
        ss >> _or >> ir >> z;
        geo = Cylinder(placement, _or, ir, z).create();
    }
    else if(shape.find("extr")!=std::string::npos) {
        int nverts;
        double v1,v2; // For Extr Poly vertices
        int nzsec;
        double zpos, off1, off2, scale; // For Extr Poly zsections
        double offset[2];
        std::vector<std::vector<double>> poly;
        std::vector<double> polyVert;
        std::vector<ExtrPoly::ZSection> zsecs;
        ss >> nverts;
        for (int i = 0; i < nverts; ++i){
            ss >> v1 >> v2;
            polyVert.push_back(v1);
            polyVert.push_back(v2);
            poly.push_back(polyVert);
            polyVert.clear();
        }
        ss >> nzsec;
        for (int i = 0; i < nzsec; ++i){
            ss >> zpos >> off1 >> off2 >> scale;
            offset[0] = off1;
            offset[1] = off2;
            zsecs.push_back(ExtrPoly::ZSection(zpos,offset,scale));
        }
        geo = ExtrPoly(placement, poly, zsecs).create();
    }
    else {
        std::stringstream ss_err;
        ss_err
            << "Shape \""
            << shape
            << "\" not recognized on line:\n"
            << ss.str();
        throw(std::runtime_error(ss_err.str()));
    }

    return geo;
}

int DetectorModel::ParseMaterialID(std::stringstream & ss, MaterialModel const & materials) {
    std::string medtype;
    ss >> medtype;

    if(not materials.HasMaterial(medtype)) {
        std::stringstream ss_err;
        ss_err
            << "Detector model uses undefined material \""
            << medtype
            << "\" on line:\n"
            << ss.str();
        throw(std::runtime_error(ss_err.str()));
    }

    return materials.GetMaterialId(medtype);
}

std::shared_ptr<siren::detector::DensityDistribution> DetectorModel::ParseDensityDistribution(std::stringstream & ss) {
    std::string distribution_type;
    ss >> distribution_type;

    std::shared_ptr<detector::DensityDistribution> density;

    if(distribution_type.find("constant") != std::string::npos) {
        double param;
        ss >> param;
        density = ConstantDensityDistribution(param).create();
    } else if (distribution_type.find("radial_polynomial") != std::string::npos) {
        double xc, yc, zc;
        ss >> xc >> yc >> zc;
        Vector3D center(xc, yc, zc);
        RadialAxis1D radial_ax(center);

        int nparams;
        ss >> nparams;

        double param;
        std::vector<double> params;
        for(int i=0; i<nparams; ++i) {
            ss >> param;
            params.push_back(param);
        }
        density = RadialAxisPolynomialDensityDistribution(radial_ax, params).create();
    } else {
        std::stringstream ss_err;
        ss_err
            << "Density distribution \""
            << distribution_type
            << "\" not recognized on line:\n"
            << ss.str();
        throw(std::runtime_error(ss_err.str()));
    }
    return density;
}

std::tuple<siren::math::Vector3D, siren::math::Quaternion> DetectorModel::ParseDetector(std::stringstream & ss) {
    std::string type, line;
    std::getline(ss, line);
    ss.clear(); ss.str(line);

    ss >> type;

    if(type.find("detector") != std::string::npos) {
        std::getline(ss, line);
        ss.clear(); ss.str(line);
    } else {
        ss.clear(); ss.str(line);
    }
    double x0, y0, z0; // Coordinates of the center of the detector
    ss >> x0 >> y0 >> z0;
    // Set the detector origin
    siren::math::Vector3D detector_origin(x0, y0, z0);

    siren::math::Quaternion detector_rotation;
    try {
        double alpha, beta, gamma; // Euler Angles of shape rotation
        ss >> alpha >> beta >> gamma;
        detector_rotation = QFromZXZr(alpha, beta, gamma);
    } catch(std::ios_base::failure e) {
        detector_rotation = siren::math::Quaternion();
    }
    return {detector_origin, detector_rotation};
}

std::shared_ptr<siren::geometry::Geometry> DetectorModel::ParseFiducialVolume(std::string fiducial_line, std::string origin_line) {
    std::string line = origin_line;
    std::stringstream ss(line);

    std::tuple<siren::math::Vector3D, siren::math::Quaternion> detector = ParseDetector(ss);

    return ParseFiducialVolume(fiducial_line, std::get<0>(detector), std::get<1>(detector));
}

std::shared_ptr<siren::geometry::Geometry> DetectorModel::ParseFiducialVolume(std::string fiducial_line, siren::math::Vector3D detector_origin, siren::math::Quaternion detector_quaternion) {
    std::string line = fiducial_line;
    std::stringstream ss(line);

    std::string type;
    ss >> type;

    bool detector_coords = true;

    if(type.find("fiducial") != std::string::npos) {
        std::getline(ss, line);
        ss.clear(); ss.str(line);
    } else {
        ss.clear(); ss.str(line);
    }

    std::string coords_type;
    ss >> coords_type;
    if(coords_type.find("detector_coords") != std::string::npos) {
        detector_coords = true;
        std::getline(ss, line);
        ss.clear(); ss.str(line);
    } else if(coords_type.find("geometry_coords") != std::string::npos) {
        detector_coords = false;
        std::getline(ss, line);
        ss.clear(); ss.str(line);
    } else {
        detector_coords = true;
        ss.clear(); ss.str(line);
    }

    std::shared_ptr<siren::geometry::Geometry> geo = ParseGeometryObject(ss);
    if(not detector_coords) {
        Placement p = geo->GetPlacement();
        p.SetPosition(detector_quaternion.rotate(p.GetPosition() - detector_origin, true));
        p.SetQuaternion(detector_quaternion.rotate(p.GetQuaternion(), true));
        geo->SetPlacement(p);
    }
    return geo;
}

void DetectorModel::LoadDetectorModel(std::string const & detector_model) {
    if(detector_model.empty())
        throw(std::runtime_error("Received empty detector model filename!"));

    std::string fname;

    if(fexists(detector_model)) {
        fname = detector_model;
    }
    else if(fexists(detector_model + ".dat")) {
        fname = detector_model + ".dat";
    }
    else if(fexists(path_ + "/densities/" + detector_model)) {
        fname = path_ + "/densities/" + detector_model;
    }
    else if(fexists(path_ + "/densities/" + detector_model + ".dat")) {
        fname = path_ + "/densities/" + detector_model + ".dat";
    }
    else if(fexists(path_ + "/Detectors/" + detector_model)) {
        fname = path_ + "/Detectors/" + detector_model;
    }
    else if(fexists(path_ + "/Detectors/" + detector_model + ".dat")) {
        fname = path_ + "/Detectors/" + detector_model + ".dat";
    }
    else if(fexists(path_ + "/" + detector_model)) {
        fname = path_ + "/" + detector_model;
    }
    else if(fexists(path_ + "/" + detector_model + ".dat")) {
        fname = path_ + "/" + detector_model + ".dat";
    }
    else {
        throw(std::runtime_error("Cannot open detector model file!"));
    }

    std::ifstream in(fname.c_str());

    // if the detectormodel file doesn't exist, stop simulation
    if(in.fail()){
        throw(std::runtime_error("Failed to open " + fname + " Set correct DetectorsPath."));
    }

    ClearSectors();
    LoadDefaultSectors();

    // read the file
    std::string buf, type;
    int level = 0;

    while(getline(in,buf)) {
        {
            size_t pos;
            // eliminate data after first #
            if((pos=buf.find('#'))!=std::string::npos)
                buf.erase(pos);
            // trim whitespace
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                // if there are no non-whitespace characters pos==std::string::npos, so the entire line is erased
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;
        }

        // density data
        std::stringstream ss(buf);
        ss >> type;

        if(type.find("object") != std::string::npos) {
            DetectorSector sector;
            sector.level = level;
            level += 1;

            sector.geo = ParseGeometryObject(ss);

            std::string label;
            ss >> label;
            sector.name = label;

            sector.material_id = ParseMaterialID(ss, materials_);

            sector.density = ParseDensityDistribution(ss);

            AddSector(sector);
        }
        else if(type.find("detector") != std::string::npos) {
            std::tuple<siren::math::Vector3D, siren::math::Quaternion> detector = ParseDetector(ss);
            detector_origin_ = std::get<0>(detector);
            detector_rotation_ = std::get<1>(detector);
        }
    } // end of the while loop
    in.close();
}

void DetectorModel::LoadDefaultMaterials() {
    // Interstellar medium mass composition from
    // https://arxiv.org/abs/astro-ph/0106359
    // Limited information exists for Z > 2
    // Approximate Z > 2 as carbon
    //
    materials_.AddMaterial(
            "VACUUM",
            std::map<int, double>(
                {
                    {1000010010, 0.704},
                    {1000020040, 0.281},
                    {1000060120, 0.015},
                }
            )
        ); // Assume there are 1 neutrons for every 7 protons in the universe
}

void DetectorModel::LoadDefaultSectors() {
    DetectorSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = std::numeric_limits<int>::min();
    sector.geo = Sphere(std::numeric_limits<double>::infinity(), 0).create();
    sector.density = ConstantDensityDistribution().create(); // Use the universe_mean_density from GEANT4
    AddSector(sector);
}

void DetectorModel::LoadMaterialModel(std::string const & material_model) {
    materials_.SetPath(path_);
    materials_.AddModelFile(material_model);
}


double DetectorModel::GetMassDensity(Geometry::IntersectionList const & intersections, GeometryPosition const & p0) const {
    Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(density >= 0);

    return density;
}

double DetectorModel::GetMassDensity(GeometryPosition const & p0) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetMassDensity(intersections, p0);
}

double DetectorModel::GetParticleDensity(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, siren::dataclasses::ParticleType target) const {
    Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            density *= materials_.GetTargetParticleFraction(sector.material_id, target);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(density >= 0);

    return density;
}

double DetectorModel::GetParticleDensity(GeometryPosition const & p0, siren::dataclasses::ParticleType target) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetParticleDensity(intersections, p0, target);
}

double DetectorModel::GetInteractionDensity(Geometry::IntersectionList const & intersections, GeometryPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const {
    Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    // If we have only decays, avoid the sector loop
    // total_decay_length is in m
    if(targets.empty()) {
        return 1.0 / total_decay_length;
    }

    double interaction_density = std::numeric_limits<double>::quiet_NaN();

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double density = sector.density->Evaluate(p0);
            std::vector<double> particle_fractions = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            interaction_density = 0.0;
            for(unsigned int i=0; i<targets.size(); ++i) {
                interaction_density += density * particle_fractions[i] * total_cross_sections[i];
            }
            interaction_density *= 100; // cm^-1 --> m^-1
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(interaction_density >= 0);

    interaction_density += 1.0 / total_decay_length;

    return interaction_density;
}

double DetectorModel::GetInteractionDensity(GeometryPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetInteractionDensity(intersections, p0, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    double dot = intersections.direction * direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    double column_depth = 0.0;

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = std::min(offset + dot * intersection->distance, distance);
        // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
        double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
        if(end_point > 0) {
            double segment_length = end_point - start_point;
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            column_depth += integral;
        }
        // last_point = end_point;
        bool done = offset + dot * intersection->distance >= distance;
        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    return column_depth * 100;
}

double DetectorModel::GetColumnDepthInCGS(GeometryPosition const & p0, GeometryPosition const & p1) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetColumnDepthInCGS(intersections, p0, p1);
}

double DetectorModel::DistanceForColumnDepthFromPoint(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryDirection const & dir, double column_depth) const {
    Vector3D direction = dir;
    column_depth /= 100;
    bool flip = column_depth < 0;
    if(column_depth < 0) {
        column_depth *= -1;
        direction = -direction;
    }

    double dot = intersections.direction * direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    double total_column_depth = 0.0;
    double total_distance = 0.0;
    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = offset + dot * intersection->distance;
        bool done = false;
        if(end_point > 0) {
            // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
            double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
            double segment_length = end_point - start_point;
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double target = column_depth - total_column_depth;
            double distance = sector.density->InverseIntegral(p0+start_point*direction, direction, target, segment_length);

            done = distance >= 0;
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            total_column_depth += integral;
            if(done) {
                total_distance = start_point + distance;
            } else {
                total_distance = start_point + segment_length;
            }
        }

        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    if(flip) {
        total_distance *= -1;
    }

    return total_distance;
}

double DetectorModel::DistanceForColumnDepthFromPoint(GeometryPosition const & p0, GeometryDirection const & direction, double column_depth) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return DistanceForColumnDepthFromPoint(intersections, p0, direction, column_depth);
}

double DetectorModel::DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryDirection const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(intersections, p0, -direction, column_depth);
}

double DetectorModel::DistanceForColumnDepthToPoint(GeometryPosition const & p0, GeometryDirection const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(p0, -direction, column_depth);
}

double DetectorModel::GetMassDensity(Geometry::IntersectionList const & intersections, GeometryPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            std::vector<double> mass_fractions = materials_.GetTargetMassFraction(sector.material_id, targets.begin(), targets.end());
            density *= std::accumulate(mass_fractions.begin(), mass_fractions.end(), 0.0);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(density >= 0);

    return density;
}

double DetectorModel::GetMassDensity(GeometryPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetMassDensity(intersections, p0, targets);
}

std::vector<double> DetectorModel::GetParticleDensity(Geometry::IntersectionList const & intersections, GeometryPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> particle_fractions;

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            particle_fractions = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    for(unsigned int i=0; i<particle_fractions.size(); ++i) {
        particle_fractions[i] *= density;
    }

    assert(density >= 0);

    return particle_fractions;
}

std::vector<double> DetectorModel::GetParticleDensity(GeometryPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetParticleDensity(intersections, p0, targets);
}

double DetectorModel::GetInteractionDepthInCGS(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    double dot = intersections.direction * direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    // If we have only decays, avoid the sector loop
    if(targets.empty()) {
      return distance/total_decay_length; // m / m --> dimensionless
    }

    std::vector<double> interaction_depths(targets.size(), 0.0);

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = std::min(offset + dot * intersection->distance, distance);
        // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
        double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
        if(end_point > 0) {
            double segment_length = end_point - start_point;
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);

            std::vector<double> particle_fractions = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            for(unsigned int i=0; i<targets.size(); ++i) {
                interaction_depths[i] += (integral * 100) * particle_fractions[i]; // cm^-3 * m --> cm^-2
            }
        }
        // last_point = end_point;
        bool done = offset + dot * intersection->distance >= distance;
        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    for(unsigned int i=0; i<targets.size(); ++i) {
        interaction_depths[i] *= total_cross_sections[i]; // cm^-2 * cm^2 == dimensionless
    }

    double interaction_depth = accumulate(interaction_depths.begin(), interaction_depths.end());

    interaction_depth += distance/total_decay_length;

    return interaction_depth;
}

std::vector<double> DetectorModel::GetParticleColumnDepth(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryPosition const & p1,  std::vector<siren::dataclasses::ParticleType> const & targets) const {
    if(p0 == p1) {
        return std::vector<double>(targets.size(), 0.0);
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return std::vector<double>(targets.size(), 0.0);
    }
    direction.normalize();

    double dot = intersections.direction * direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }


    std::vector<double> target_counts(targets.size(), 0);

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = std::min(offset + dot * intersection->distance, distance);
        // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
        double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
        if(end_point > 0) {
            double segment_length = end_point - start_point;
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            std::vector<double> particle_fractions = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            for(unsigned int i=0; i<target_counts.size(); ++i) {
                target_counts[i] += (integral * 100) * particle_fractions[i];
            }
        }
        // last_point = end_point;
        bool done = offset + dot * intersection->distance >= distance;
        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    return target_counts;
}

double DetectorModel::GetInteractionDepthInCGS(GeometryPosition const & p0, GeometryPosition const & p1,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetInteractionDepthInCGS(intersections, p0, p1, targets, total_cross_sections, total_decay_length);
}

DetectorSector DetectorModel::GetContainingSector(Geometry::IntersectionList const & intersections, GeometryPosition const & p0) const {
    Vector3D direction = intersections.direction;

    double offset = (intersections.position - p0) * direction;
    double dot = (intersections.position - p0) * (intersections.position - p0);

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    DetectorSector sector;

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        double end_point = offset + dot * intersection->distance;
        double start_point = offset + dot * current_intersection->distance;
        bool done = false;
        if((start_point < 0 and end_point > 0) or start_point == 0) {
            sector = GetSector(current_intersection->hierarchy);
            done = true;
        }
        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    return sector;
}

DetectorSector DetectorModel::GetContainingSector(GeometryPosition const & p0) const {
    Vector3D direction(0, 0, 1);
    Geometry::IntersectionList intersections = GetIntersections(p0, GeometryDirection(direction));
    return GetContainingSector(intersections, p0);
}

Geometry::IntersectionList DetectorModel::GetIntersections(GeometryPosition const & p0, GeometryDirection const & direction) const {
    Geometry::IntersectionList intersections;
    intersections.position = p0;
    intersections.direction = direction;

    // Obtain the intersections with each sector geometry
    for(auto const & sector : sectors_) {
        std::vector<Geometry::Intersection> i = sector.geo->Intersections(p0, direction);
        intersections.intersections.reserve(intersections.intersections.size() + std::distance(i.begin(), i.end()));
        intersections.intersections.insert(intersections.intersections.end(), i.begin(), i.end());
        for(unsigned int j=intersections.intersections.size(); j>intersections.intersections.size()-i.size(); --j) {
            intersections.intersections[j-1].hierarchy = sector.level;
            intersections.intersections[j-1].matID = sector.material_id;
        }
    }

    SortIntersections(intersections);

    return intersections;
}

void DetectorModel::SortIntersections(Geometry::IntersectionList & intersections) {
    SortIntersections(intersections.intersections);
}

void DetectorModel::SortIntersections(std::vector<Geometry::Intersection> & intersections) {
    // Intersections should be sorted according to distance and then hierarchy
    std::function<bool(Geometry::Intersection const &, Geometry::Intersection const &)> comp = [](Geometry::Intersection const & a, Geometry::Intersection const & b){
        bool a_enter = a.entering;
        bool b_enter = b.entering;
        if(a.distance < b.distance)
            return true;
        else if(a.distance == b.distance) {
            bool low_high = a.hierarchy < b.hierarchy;
            bool high_low = b.hierarchy < a.hierarchy;
            if(a_enter) {
                if(b_enter)
                    return high_low;
                else
                    return false;
            }
            else {
                if(b_enter)
                    return true;
                else
                    return low_high;
            }
        }
        else
            return false;
    };
    std::sort(intersections.begin(), intersections.end(), comp);
}

Geometry::IntersectionList DetectorModel::GetOuterBounds(Geometry::IntersectionList const & intersections) {
    Geometry::IntersectionList result;
    result.position = intersections.position;
    result.direction = intersections.direction;
    int min_hierarchy = std::numeric_limits<int>::min();
    ptrdiff_t min_index = 0;
    for(ptrdiff_t i=0; i<intersections.intersections.size(); ++i) {
        if(intersections.intersections[i].hierarchy > min_hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
            min_index = i;
            break;
        }
    }
    for(ptrdiff_t i=ptrdiff_t(intersections.intersections.size())-1; (i >= 0 and i > min_index); --i) {
        if(intersections.intersections[i].hierarchy > min_hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
            break;
        }
    }
    return result;
}

Geometry::IntersectionList DetectorModel::GetOuterBounds(GeometryPosition const & p0, GeometryDirection const & direction) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetOuterBounds(intersections);
}

std::set<siren::dataclasses::ParticleType> DetectorModel::GetAvailableTargets(GeometryPosition const & vertex) const {
    Geometry::IntersectionList intersections = GetIntersections(vertex, GeometryDirection(math::Vector3D(0,0,1)));
    return GetAvailableTargets(intersections, vertex);
}

std::set<siren::dataclasses::ParticleType> DetectorModel::GetAvailableTargets(geometry::Geometry::IntersectionList const & intersections, GeometryPosition const & vertex) const {
    int matID = GetContainingSector(intersections, vertex).material_id;
    std::vector<siren::dataclasses::ParticleType> particles = materials_.GetMaterialConstituents(matID);
    return std::set<siren::dataclasses::ParticleType>(particles.begin(), particles.end());
}

void DetectorModel::SectorLoop(std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback, Geometry::IntersectionList const & intersections, bool reverse) {
    // Keep track of the integral progress
    double last_point;

    // Keep track of the sectors our integrand is inside
    std::map<int, std::vector<Geometry::Intersection>::const_iterator> stack;

    // Keep track of the entry point to the relevant sector
    std::vector<Geometry::Intersection>::const_iterator current_intersection = intersections.intersections.begin();

    if(current_intersection == intersections.intersections.end()) {
        return;
    }

    int start_index;
    int end_index;
    int increment;
    if(reverse) {
        end_index = -1;
        start_index = std::max(int(intersections.intersections.size()-2), int(end_index));
        increment = -1;
        current_intersection = intersections.intersections.begin() + start_index + 1;
    } else {
        end_index = intersections.intersections.size();
        start_index = std::min(int(1), int(end_index));
        increment = 1;
        current_intersection = intersections.intersections.begin() + start_index - 1;
    }

    last_point = current_intersection->distance;

    // Integration only begins once we are inside a sector
    stack.insert({current_intersection->hierarchy, current_intersection});

    for(int i=start_index; i!=end_index; i += increment) {
        // The transition point into the next sector
        std::vector<Geometry::Intersection>::const_iterator intersection = intersections.intersections.begin() + i;
        if(intersection == intersections.intersections.end())
            throw(std::runtime_error("Reached end of intersections! Should never reach this point!"));
        if(intersection->entering ^ reverse) {
            // Entering a sector means it is added to the stack
            stack.insert({intersection->hierarchy, intersection});
            // A sector transition only occurs if the intersection is with a sector of larger hierarchy
            if(intersection->hierarchy > current_intersection->hierarchy) {
                bool done = callback(current_intersection, intersection, last_point);
                last_point = intersection->distance;
                if(done) {
                    break;
                }
                current_intersection = intersection;
            }
        }
        else {
            if(intersection->hierarchy <= current_intersection->hierarchy) {
                // Exiting a sector means we remove it from the stack
                stack.erase(intersection->hierarchy);
                // If the intersection boundary is from the same hierarchy, we perform the integral
                if(intersection->hierarchy == current_intersection->hierarchy) {
                    bool done = callback(current_intersection, intersection, last_point);
                    last_point = intersection->distance;
                    if(done) {
                        break;
                    }
                    if(stack.size() > 0) {
                        // Exiting the current sector means we need to move one level down and grab the entry to that sector
                        auto lb = --stack.lower_bound(intersection->hierarchy);
                        assert(lb != stack.end());
                        current_intersection = lb->second;
                    } else {
                        // An empty stack means we have no intersections left to work with
                        break;
                    }
                }
                // If the intersection boundary is from a lower hierarchy,
                // then we can wait to perform the integration,
                // since this intersection does not represent a physical transition to a different sector
            }
            else {
                // If we are exiting a sector with larger hierarchy the current_intersection should have been set to match that sector.
                // Thus, we would not reach this point.
                throw(std::runtime_error("Cannot exit a level that we have not entered!"));
            }
        }
    }
}

double DetectorModel::DistanceForInteractionDepthFromPoint(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryDirection const & dir, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    Vector3D direction = dir;
    bool flip = interaction_depth < 0;
    if(interaction_depth < 0) {
        interaction_depth *= -1;
        direction = -direction;
    }

    double dot = intersections.direction * direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    // If we have only decays, avoid the sector loop
    // decay length in m
    if(targets.empty()) {
      return interaction_depth * total_decay_length;
    }

    // Recast decay length to cm for density integral
    double total_decay_length_cm = total_decay_length / siren::utilities::Constants::cm;

    double total_interaction_depth = 0.0;
    double total_distance = 0.0;
    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = offset + dot * intersection->distance;
        bool done = false;
        if(end_point > 0) {
            // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
            double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
            double segment_length = end_point - start_point;
            DetectorSector sector = GetSector(current_intersection->hierarchy);
            double target = interaction_depth - total_interaction_depth;
            // This next line is because when we evaluate the density integral,
            // we end up calculating an interaction length in units of m/cm.
            // This is a correction
            target /= 100;
            std::vector<double> interaction_depths = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            for(unsigned int i=0; i<targets.size(); ++i) {
                interaction_depths[i] *= total_cross_sections[i];
            }
            double target_composition = accumulate(interaction_depths.begin(), interaction_depths.end(), 0.0); // cm^2 g^-1
            target /= target_composition;
            double distance;
            // total_decay_length now in cm
            if (total_decay_length < std::numeric_limits<double>::infinity()) {
              distance = sector.density->InverseIntegral(p0+start_point*direction, direction, 1./(total_decay_length_cm*target_composition), target, segment_length);
            }
            else {
              distance = sector.density->InverseIntegral(p0+start_point*direction, direction, target, segment_length);
            }
            done = distance >= 0;
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length); // g cm^-3 * m
            integral *= (target_composition*siren::utilities::Constants::m/siren::utilities::Constants::cm); // --> m cm^-1 --> dimensionless
            total_interaction_depth += integral;
            if(done) {
                total_distance = start_point + distance;
            } else {
                total_distance = start_point + segment_length;
            }
        }

        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    if(flip) {
        total_distance *= -1;
    }

    return total_distance;
}

double DetectorModel::DistanceForInteractionDepthFromPoint(GeometryPosition const & p0, GeometryDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return DistanceForInteractionDepthFromPoint(intersections, p0, direction, interaction_depth, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::DistanceForInteractionDepthToPoint(Geometry::IntersectionList const & intersections, GeometryPosition const & p0, GeometryDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(intersections, p0, -direction, interaction_depth, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::DistanceForInteractionDepthToPoint(GeometryPosition const & p0, GeometryDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(p0, -direction, interaction_depth, targets, total_cross_sections, total_decay_length);
}

//////////////////////////////////////////////////////////


double DetectorModel::GetMassDensity(Geometry::IntersectionList const & intersections, DetectorPosition const & p0) const {
    return GetMassDensity(intersections, ToGeo(p0));
}

double DetectorModel::GetMassDensity(DetectorPosition const & p0) const {
    return GetMassDensity(ToGeo(p0));
}

double DetectorModel::GetParticleDensity(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, siren::dataclasses::ParticleType target) const {
    return GetParticleDensity(intersections, ToGeo(p0), target);
}

double DetectorModel::GetParticleDensity(DetectorPosition const & p0, siren::dataclasses::ParticleType target) const {
    return GetParticleDensity(ToGeo(p0), target);
}

double DetectorModel::GetInteractionDensity(Geometry::IntersectionList const & intersections, DetectorPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const {
    return GetInteractionDensity(intersections, ToGeo(p0), targets, total_cross_sections, total_decay_length);
}

double DetectorModel::GetInteractionDensity(DetectorPosition const & p0,
            std::vector<siren::dataclasses::ParticleType> const & targets,
            std::vector<double> const & total_cross_sections,
            double const & total_decay_length) const {
    return GetInteractionDensity(ToGeo(p0), targets, total_cross_sections, total_decay_length);
}

double DetectorModel::GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1) const {
    return GetColumnDepthInCGS(intersections, ToGeo(p0), ToGeo(p1));
}

double DetectorModel::GetColumnDepthInCGS(DetectorPosition const & p0, DetectorPosition const & p1) const {
    return GetColumnDepthInCGS(ToGeo(p0), ToGeo(p1));
}

double DetectorModel::DistanceForColumnDepthFromPoint(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorDirection const & dir, double column_depth) const {
    return DistanceForColumnDepthFromPoint(intersections, ToGeo(p0), ToGeo(dir), column_depth);
}

double DetectorModel::DistanceForColumnDepthFromPoint(DetectorPosition const & p0, DetectorDirection const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(ToGeo(p0), ToGeo(direction), column_depth);
}

double DetectorModel::DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorDirection const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(intersections, ToGeo(p0), ToGeo(direction), column_depth);
}

double DetectorModel::DistanceForColumnDepthToPoint(DetectorPosition const & p0, DetectorDirection const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(ToGeo(p0), ToGeo(direction), column_depth);
}

double DetectorModel::GetMassDensity(Geometry::IntersectionList const & intersections, DetectorPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    return GetMassDensity(intersections, ToGeo(p0), targets);
}

double DetectorModel::GetMassDensity(DetectorPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    return GetMassDensity(ToGeo(p0), targets);
}

std::vector<double> DetectorModel::GetParticleDensity(Geometry::IntersectionList const & intersections, DetectorPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    return GetParticleDensity(intersections, ToGeo(p0), targets);
}

std::vector<double> DetectorModel::GetParticleDensity(DetectorPosition const & p0,  std::set<siren::dataclasses::ParticleType> targets) const {
    return GetParticleDensity(ToGeo(p0), targets);
}

double DetectorModel::GetInteractionDepthInCGS(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return GetInteractionDepthInCGS(intersections, ToGeo(p0), ToGeo(p1), targets, total_cross_sections, total_decay_length);
}

std::vector<double> DetectorModel::GetParticleColumnDepth(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorPosition const & p1,  std::vector<siren::dataclasses::ParticleType> const & targets) const {
    return GetParticleColumnDepth(intersections, ToGeo(p0), ToGeo(p1), targets);
}

double DetectorModel::GetInteractionDepthInCGS(DetectorPosition const & p0, DetectorPosition const & p1,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return GetInteractionDepthInCGS(ToGeo(p0), ToGeo(p1), targets, total_cross_sections, total_decay_length);
}

DetectorSector DetectorModel::GetContainingSector(Geometry::IntersectionList const & intersections, DetectorPosition const & p0) const {
    return GetContainingSector(intersections, ToGeo(p0));
}

DetectorSector DetectorModel::GetContainingSector(DetectorPosition const & p0) const {
    return GetContainingSector(ToGeo(p0));
}

Geometry::IntersectionList DetectorModel::GetIntersections(DetectorPosition const & p0, DetectorDirection const & direction) const {
    return GetIntersections(ToGeo(p0), ToGeo(direction));
}

Geometry::IntersectionList DetectorModel::GetOuterBounds(DetectorPosition const & p0, DetectorDirection const & direction) const {
    return GetOuterBounds(ToGeo(p0), ToGeo(direction));
}

std::set<siren::dataclasses::ParticleType> DetectorModel::GetAvailableTargets(DetectorPosition const & vertex) const {
    return GetAvailableTargets(ToGeo(vertex));
}

std::set<siren::dataclasses::ParticleType> DetectorModel::GetAvailableTargets(geometry::Geometry::IntersectionList const & intersections, DetectorPosition const & vertex) const {
    return GetAvailableTargets(intersections, ToGeo(vertex));
}


double DetectorModel::DistanceForInteractionDepthFromPoint(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorDirection const & dir, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(intersections, ToGeo(p0), ToGeo(dir), interaction_depth, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::DistanceForInteractionDepthFromPoint(DetectorPosition const & p0, DetectorDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(ToGeo(p0), ToGeo(direction), interaction_depth, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::DistanceForInteractionDepthToPoint(Geometry::IntersectionList const & intersections, DetectorPosition const & p0, DetectorDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(intersections, ToGeo(p0), ToGeo(direction), interaction_depth, targets, total_cross_sections, total_decay_length);
}

double DetectorModel::DistanceForInteractionDepthToPoint(DetectorPosition const & p0, DetectorDirection const & direction, double interaction_depth,
        std::vector<siren::dataclasses::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections,
        double const & total_decay_length) const {
    return DistanceForInteractionDepthFromPoint(ToGeo(p0), ToGeo(direction), interaction_depth, targets, total_cross_sections, total_decay_length);
}

//////////////////////////////////////////////////////////

void DetectorModel::LoadConcentricShellsFromLegacyFile(std::string model_fname, double detector_depth, double ice_cap_angle) {
    if(model_fname.empty())
        throw(std::runtime_error("Received empty detector model filename!"));

    std::string fname;

    if(fexists(model_fname)) {
        fname = model_fname;
    }
    else if(fexists(model_fname + ".dat")) {
        fname = model_fname + ".dat";
    }
    else if(fexists(path_ + "/densities/" + model_fname)) {
        fname = path_ + "/densities/" + model_fname;
    }
    else if(fexists(path_ + "/densities/" + model_fname + ".dat")) {
        fname = path_ + "/densities/" + model_fname + ".dat";
    }
    else if(fexists(path_ + "/Detectors/" + model_fname)) {
        fname = path_ + "/Detectors/" + model_fname;
    }
    else if(fexists(path_ + "/Detectors/" + model_fname + ".dat")) {
        fname = path_ + "/Detectors/" + model_fname + ".dat";
    }
    else if(fexists(path_ + "/" + model_fname)) {
        fname = path_ + "/" + model_fname;
    }
    else if(fexists(path_ + "/" + model_fname + ".dat")) {
        fname = path_ + "/" + model_fname + ".dat";
    }
    else {
        throw(std::runtime_error("Cannot open detector model file!"));
    }

    std::ifstream in(fname.c_str());

    // if the detectormodel file doesn't exist, stop simulation
    if(in.fail()){
        throw(std::runtime_error("Failed to open " + fname + " Set correct DetectorsPath."));
    }

    ClearSectors();
    LoadDefaultSectors();

    // read the file
    std::string buf;
    std::string label, medtype;
    double radius, param;
    int nparams;

    int level = -sectors_.size();
    double max_radius = 0;
    while(getline(in,buf)) {
        {
            size_t pos;
            // eliminate data after first #
            if((pos=buf.find('#'))!=std::string::npos)
                buf.erase(pos);
            // trim whitespace
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                // if there are no non-whitespace characters pos==std::string::npos, so the entire line is erased
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;
        }

        // density data
        std::stringstream ss(buf);
        ss >> radius >> label >> medtype >> nparams;

        if(not materials_.HasMaterial(medtype)) {
            std::stringstream ss;
            ss << "Detector model uses undefined material " << medtype;
            throw(std::runtime_error(ss.str()));
        }

        DetectorSector sector;
        sector.material_id = materials_.GetMaterialId(medtype);
        sector.level = level;
        sector.name = label;
        sector.geo = Sphere(radius, 0).create();
        level -= 1;
        if(nparams == 1) {
            ss >> param;
            sector.density = ConstantDensityDistribution(param).create();
        }
        else {
            std::vector<double> params;
            for(int i=0; i<nparams; ++i) {
                ss >> param;
                params.push_back(param);
            }
            RadialAxis1D radial_ax;
            sector.density = RadialAxisPolynomialDensityDistribution(radial_ax, params).create();
        }

        // stop the process if layering assumptions are violated
        if(radius < max_radius) {
            throw(std::runtime_error("Layers must be radially ordered in file!"));
        }
        max_radius = radius;
        AddSector(sector);
    } // end of the while loop
    in.close();

    // Examine the ice
    double earth_radius = 0;
    double ice_radius = 0;
    std::vector<int> ice_layers;
    bool saw_ice = false;
    for(unsigned int i=1; i<sectors_.size(); ++i) {
        DetectorSector const & sector = sectors_[i];
        std::string name = materials_.GetMaterialName(sector.material_id);
        string_to_lower(name);

        bool in_ice = name == "ice";
        bool solid = (name != "air") and (name != "atmosphere") and (name != "vacuum");
        saw_ice |= in_ice;

        if(not saw_ice) {
            // In the Earth, keep increasing the radius
            if(solid)
                earth_radius = ((Sphere *)(sector.geo.get()))->GetRadius();
        }
        else if(in_ice) {
            // In the ice, keep increasing the radius
            ice_radius = ((Sphere *)(sector.geo.get()))->GetRadius();
            ice_layers.push_back(i);
        }
        else {
            // Out of the ice, stop counting layers
            break;
        }
    }

    // Set the detector origin
    // Depth is defined relative to the top solid layer
    detector_origin_ = Vector3D(0,0,std::max(earth_radius, ice_radius)-detector_depth);

    if(ice_cap_angle < 0 || ice_cap_angle > siren::utilities::Constants::pi) {
        // Leave any ice as is
    }
    else { // Setup the ice cap

        if(ice_layers.size() > 0) {
            // calculate radius
            double costheta = cos(ice_cap_angle);
            //double sintheta = sqrt((1-costheta)*(1+costheta));

            // get largest depth of ice
            double h = ice_radius - earth_radius;
            // h : max depth of ice at south pole
            // d + h = radius of icecap
            // r = bed_rock radius
            // r - d : z-position of center of icecap sphere
            //
            // (d + h)^2 = (r - d)^2 + r^2
            //                     - 2*r*(r-d)*costheta
            // d = (r^2*(1-costheta) - 0.5*h^2) / (h + r*(1-costheta))
            //
            double r = earth_radius;
            double d = (r*r*(1 - costheta) - 0.5*h*h) / (h + r*(1 - costheta));

            ice_radius = d + h; // radius of sphere of icecap
            double ice_offset = r - d; // z-pos of center of sphere of icecap

            for(auto const & i : ice_layers) {
                DetectorSector & sector = sectors_[i];
                Sphere const * geo = dynamic_cast<Sphere const *>(sector.geo.get());
                sector.geo = Sphere(Placement(Vector3D(0,0,ice_offset), QFromZXZr(0,0,0)), geo->GetRadius()-ice_offset, 0).create();
                //geo->SetRadius(geo->GetRadius()-ice_offset);
                //geo->SetPosition(Vector3D(0,0,ice_offset));
            }
        }
    }
}

double DetectorModel::GetTargetMass(siren::dataclasses::ParticleType target) const {
    double molar_mass = materials_.GetMolarMass(target); // grams per mole
    return molar_mass * siren::utilities::Constants::GeV_per_amu;
}
