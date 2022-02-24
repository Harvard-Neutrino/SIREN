#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "earthmodel-service/EarthModel.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EulerQuaternionConversions.h"

using namespace earthmodel;

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

EarthModel::EarthModel() {
    LoadDefaultMaterials();
    LoadDefaultSectors();
}

EarthModel::EarthModel(std::string const & earth_model, std::string const & material_model) {
    LoadDefaultMaterials();
    LoadDefaultSectors();
    LoadMaterialModel(material_model);
    LoadEarthModel(earth_model);
}

EarthModel::EarthModel(std::string const & path, std::string const & earth_model, std::string const & material_model) : path_(path) {
    LoadDefaultMaterials();
    LoadDefaultSectors();
    LoadMaterialModel(material_model);
    LoadEarthModel(earth_model);
}

std::string EarthModel::GetPath() const {
    return path_;
}

void EarthModel::SetPath(std::string const & path) {
    path_ = path;
}

MaterialModel const & EarthModel::GetMaterials() const {
    return materials_;
}

void EarthModel::SetMaterials(MaterialModel const & materials) {
    materials_ = materials;
}

std::vector<EarthSector> const & EarthModel::GetSectors() const {
    return sectors_;
}

void EarthModel::SetSectors(std::vector<EarthSector> const & sectors) {
    sectors_ = sectors;
}

Vector3D EarthModel::GetDetectorOrigin() const {
    return detector_origin_;
}

void EarthModel::SetDetectorOrigin(Vector3D const & detector_origin) {
    detector_origin_ = detector_origin;
}

void EarthModel::AddSector(EarthSector sector) {
    if(sector_map_.count(sector.level) > 0) {
        throw(std::runtime_error("Already have a sector of that heirarchy!"));
    }
    else {
        sector_map_[sector.level] = sectors_.size();
        sectors_.push_back(sector);
    }
}

EarthSector EarthModel::GetSector(int heirarchy) const {
    auto const iter = sector_map_.find(heirarchy);
    assert(iter != sector_map_.end());
    unsigned int index = sector_map_.at(heirarchy);
    assert(index < sectors_.size());
    unsigned int alt_index = sector_map_.find(heirarchy)->second;
    assert(index == alt_index);
    return sectors_[index];
}

void EarthModel::ClearSectors() {
    sectors_.clear();
    sector_map_.clear();
}

namespace {
bool fexists(const char *filename)
{
    std::ifstream ifile(filename);
    return (bool)ifile;
}
bool fexists(const std::string filename)
{
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}
}

void EarthModel::LoadEarthModel(std::string const & earth_model) {
    if(earth_model.empty())
        throw(std::runtime_error("Received empty earth model filename!"));

    std::string fname;

    if(fexists(earth_model)) {
        fname = earth_model;
    }
    else if(fexists(earth_model + ".dat")) {
        fname = earth_model + ".dat";
    }
    else if(fexists(path_ + "/densities/" + earth_model)) {
        fname = path_ + "/densities/" + earth_model;
    }
    else if(fexists(path_ + "/densities/" + earth_model + ".dat")) {
        fname = path_ + "/densities/" + earth_model + ".dat";
    }
    else if(fexists(path_ + "/earthparams/" + earth_model)) {
        fname = path_ + "/earthparams/" + earth_model;
    }
    else if(fexists(path_ + "/earthparams/" + earth_model + ".dat")) {
        fname = path_ + "/earthparams/" + earth_model + ".dat";
    }
    else if(fexists(path_ + "/" + earth_model)) {
        fname = path_ + "/" + earth_model;
    }
    else if(fexists(path_ + "/" + earth_model + ".dat")) {
        fname = path_ + "/" + earth_model + ".dat";
    }
    else {
        throw(std::runtime_error("Cannot open earth model file!"));
    }

    std::ifstream in(fname.c_str());

    // if the earthmodel file doesn't exist, stop simulation
    if(in.fail()){
        throw(std::runtime_error("Failed to open " + fname + " Set correct EarthParamsPath."));
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
            EarthSector sector;
            sector.level = level;
            level += 1;
            double xc, yc, zc; // Coordinates of the center of the shape
            double alpha, beta, gamma; // Euler Angles of shape rotation

            std::string shape;
            ss >> shape;
            ss >> xc >> yc >> zc;
            ss >> alpha >> beta >> gamma;
            Placement placement(Vector3D(xc,yc,zc), QFromZXZr(alpha,beta,gamma));

            if(shape.find("sphere")!=std::string::npos) {
                double radius; // For Sphere shapes
                ss >> radius;
                sector.geo = Sphere(placement, radius, 0).create();
            }
            else if(shape.find("box")!=std::string::npos) {
                double dx, dy, dz; // For Box shapes
                ss >> dx >> dy >> dz;
                sector.geo = Box(placement, dx, dy, dz).create();
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
                sector.geo = ExtrPoly(placement, poly, zsecs).create();
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

            std::string label, medtype;
            ss >> label >> medtype;

            if(not materials_.HasMaterial(medtype)) {
                std::stringstream ss_err;
                ss_err
                    << "Earth model uses undefined material \""
                    << medtype
                    << "\" on line:\n"
                    << ss.str();
                throw(std::runtime_error(ss_err.str()));
            }

            sector.material_id = materials_.GetMaterialId(medtype);

            std::string distribution_type;
            ss >> distribution_type;

            if(distribution_type.find("constant") != std::string::npos) {
                double param;
                ss >> param;
                sector.density = DensityDistribution1D<CartesianAxis1D,ConstantDistribution1D>(param).create();
            } else if (distribution_type.find("radial_polynomial") != std::string::npos) {
                double xc, yc, zc;
                ss >> xc, yc, zc;
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
                sector.density = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(radial_ax, params).create();
            } else {
                std::stringstream ss_err;
                ss_err
                    << "Density distribution \""
                    << distribution_type
                    << "\" not recognized on line:\n"
                    << ss.str();
                throw(std::runtime_error(ss_err.str()));
            }

            AddSector(sector);
        }
        else if(type.find("detector") != std::string::npos) {
            double x0, y0, z0; // Coordinates of the center of the detector
            ss >> x0 >> y0 >> z0;
            // Set the detector origin
            detector_origin_ = Vector3D(x0, y0, z0);
        }
    } // end of the while loop
    in.close();
}

void EarthModel::LoadDefaultMaterials() {
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

void EarthModel::LoadDefaultSectors() {
    EarthSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = std::numeric_limits<int>::min();
    sector.geo = Sphere(std::numeric_limits<double>::infinity(), 0).create();
    sector.density = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>().create(); // Use the universe_mean_density from GEANT4
    AddSector(sector);
}

void EarthModel::LoadMaterialModel(std::string const & material_model) {
    materials_.SetPath(path_);
    materials_.AddModelFile(material_model);
}


double EarthModel::GetMassDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

double EarthModel::GetMassDensity(Vector3D const & p0) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetMassDensity(intersections, p0);
}

double EarthModel::GetParticleDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, LeptonInjector::Particle::ParticleType target) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

double EarthModel::GetParticleDensity(Vector3D const & p0, LeptonInjector::Particle::ParticleType target) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetParticleDensity(intersections, p0, target);
}

double EarthModel::GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

double EarthModel::GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetColumnDepthInCGS(intersections, p0, p1);
}

double EarthModel::DistanceForColumnDepthFromPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & dir, double column_depth) const {
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
    double total_distance = std::numeric_limits<double>::quiet_NaN();
    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = offset + dot * intersection->distance;
        bool done = false;
        if(end_point > 0) {
            // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
            double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
            double segment_length = end_point - start_point;
            EarthSector sector = GetSector(current_intersection->hierarchy);
            double target = column_depth - total_column_depth;
            double distance = sector.density->InverseIntegral(p0+start_point*direction, direction, target, segment_length);
            done = distance >= 0;
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            total_column_depth += integral;
            if(done) {
                total_distance = start_point + distance;
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

double EarthModel::DistanceForColumnDepthFromPoint(Vector3D const & p0, Vector3D const & direction, double column_depth) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return DistanceForColumnDepthFromPoint(intersections, p0, direction, column_depth);
}

double EarthModel::DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(intersections, p0, -direction, column_depth);
}

double EarthModel::DistanceForColumnDepthToPoint(Vector3D const & p0, Vector3D const & direction, double column_depth) const {
    return DistanceForColumnDepthFromPoint(p0, -direction, column_depth);
}

double EarthModel::GetMassDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0,  std::set<LeptonInjector::Particle::ParticleType> targets) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

double EarthModel::GetMassDensity(Vector3D const & p0,  std::set<LeptonInjector::Particle::ParticleType> targets) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetMassDensity(intersections, p0, targets);
}

std::vector<double> EarthModel::GetParticleDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0,  std::set<LeptonInjector::Particle::ParticleType> targets) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

std::vector<double> EarthModel::GetParticleDensity(Vector3D const & p0,  std::set<LeptonInjector::Particle::ParticleType> targets) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetParticleDensity(intersections, p0, targets);
}

double EarthModel::GetInteractionDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
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

    std::vector<double> interaction_depths(targets.size(), 0.0);

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = std::min(offset + dot * intersection->distance, distance);
        // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
        double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
        if(end_point > 0) {
            double segment_length = end_point - start_point;
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

    return interaction_depth;
}

std::vector<double> EarthModel::GetParticleColumnDepth(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1,  std::vector<LeptonInjector::Particle::ParticleType> const & targets) const {
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
            EarthSector sector = GetSector(current_intersection->hierarchy);
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

double EarthModel::GetInteractionDepthInCGS(Vector3D const & p0, Vector3D const & p1,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
    if(p0 == p1) {
        return 0.0;
    }
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    if(distance == 0.0) {
        return 0.0;
    }
    direction.normalize();

    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetInteractionDepthInCGS(intersections, p0, p1, targets, total_cross_sections);
}

EarthSector EarthModel::GetContainingSector(Geometry::IntersectionList const & intersections, Vector3D const & p0) const {
    Vector3D direction = intersections.direction;

    double offset = (intersections.position - p0) * direction;
    double dot = (intersections.position - p0) * (intersections.position - p0);

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }

    EarthSector sector;

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

EarthSector EarthModel::GetContainingSector(Vector3D const & p0) const {
    Vector3D direction(0, 0, 1);
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetContainingSector(intersections, p0);
}

Geometry::IntersectionList EarthModel::GetIntersections(Vector3D const & p0, Vector3D const & direction) const {
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

void EarthModel::SortIntersections(Geometry::IntersectionList & intersections) {
    SortIntersections(intersections.intersections);
}

void EarthModel::SortIntersections(std::vector<Geometry::Intersection> & intersections) {
    // Intersections should be sorted according to distance and then hierarchy
    std::function<bool(Geometry::Intersection const &, Geometry::Intersection const &)> comp = [](Geometry::Intersection const & a, Geometry::Intersection const & b){
		bool a_enter = a.entering;
		bool b_enter = b.entering;
        if(a.distance < b.distance)
            return true;
        else if(a.distance == b.distance) {
            bool low_high = a.hierarchy < b.hierarchy;
            if(a_enter) {
                if(b_enter)
                    return not low_high;
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

Geometry::IntersectionList EarthModel::GetOuterBounds(Geometry::IntersectionList const & intersections) {
    Geometry::IntersectionList result;
    result.position = intersections.position;
    result.direction = intersections.direction;
    int min_hierarchy = std::numeric_limits<int>::min();
    int min_index = 0;
    for(unsigned int i=0; i<intersections.intersections.size(); ++i) {
        if(intersections.intersections[i].hierarchy > min_hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
            min_index = i;
            break;
        }
    }
    for(unsigned int i=intersections.intersections.size()-1; (i >= 0 and i > min_index); --i) {
        if(intersections.intersections[i].hierarchy > min_hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
            break;
        }
    }
    return result;
}

Geometry::IntersectionList EarthModel::GetOuterBounds(Vector3D const & p0, Vector3D const & direction) {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetOuterBounds(intersections);
}

std::set<LeptonInjector::Particle::ParticleType> EarthModel::GetAvailableTargets(std::array<double,3> const & vertex) {
		int matID = GetContainingSector(Vector3D(vertex[0],vertex[1],vertex[2])).material_id;
        std::vector<LeptonInjector::Particle::ParticleType> particles = materials_.GetMaterialConstituents(matID);
        return std::set<LeptonInjector::Particle::ParticleType>(particles.begin(), particles.end());
}


void EarthModel::SectorLoop(std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback, Geometry::IntersectionList const & intersections, bool reverse) {
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

    for(unsigned int i=start_index; i!=end_index; i += increment) {
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

double EarthModel::DistanceForInteractionDepthFromPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & dir, double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
    Vector3D direction = dir;
    interaction_depth /= 100;
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

    double total_interaction_depth = 0.0;
    double total_distance = std::numeric_limits<double>::quiet_NaN();
    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection and the global integral boundary
        double end_point = offset + dot * intersection->distance;
        bool done = false;
        if(end_point > 0) {
            // whereas the lower end is bounded by the global start point, the end of the last line segment, and the entry into the sector
            double start_point = std::max(std::max(offset + dot * current_intersection->distance, 0.0), offset + dot * last_point);
            double segment_length = end_point - start_point;
            EarthSector sector = GetSector(current_intersection->hierarchy);
            double target = interaction_depth - total_interaction_depth;
            std::vector<double> interaction_depths = materials_.GetTargetParticleFraction(sector.material_id, targets.begin(), targets.end());
            for(unsigned int i=0; i<targets.size(); ++i) {
                interaction_depths[i] *= total_cross_sections[i];
            }
            double target_composition = accumulate(interaction_depths.begin(), interaction_depths.end(), 0.0); // g * cm^-3
            target /= target_composition;
            double distance = sector.density->InverseIntegral(p0+start_point*direction, direction, target, segment_length);
            done = distance >= 0;
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            integral *= target_composition;
            total_interaction_depth += integral;
            if(done) {
                total_distance = start_point + distance;
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

double EarthModel::DistanceForInteractionDepthFromPoint(Vector3D const & p0, Vector3D const & direction, double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return DistanceForInteractionDepthFromPoint(intersections, p0, direction, interaction_depth, targets, total_cross_sections);
}

double EarthModel::DistanceForInteractionDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & direction, double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
    return DistanceForInteractionDepthFromPoint(intersections, p0, -direction, interaction_depth, targets, total_cross_sections);
}

double EarthModel::DistanceForInteractionDepthToPoint(Vector3D const & p0, Vector3D const & direction, double interaction_depth,
        std::vector<LeptonInjector::Particle::ParticleType> const & targets,
        std::vector<double> const & total_cross_sections) const {
    return DistanceForInteractionDepthFromPoint(p0, -direction, interaction_depth, targets, total_cross_sections);
}

Vector3D EarthModel::GetEarthCoordPosFromDetCoordPos(Vector3D const & point) const {
    return point + detector_origin_;
}

Vector3D EarthModel::GetEarthCoordDirFromDetCoordDir(Vector3D const & direction) const {
    return direction;
}

Vector3D EarthModel::GetDetCoordPosFromEarthCoordPos(Vector3D const & point) const {
    return point - detector_origin_;
}

Vector3D EarthModel::GetDetCoordDirFromEarthCoordDir(Vector3D const & direction) const {
   return direction;
}

void EarthModel::LoadConcentricShellsFromLegacyFile(std::string model_fname, double detector_depth, double ice_cap_angle) {
    if(model_fname.empty())
        throw(std::runtime_error("Received empty earth model filename!"));

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
    else if(fexists(path_ + "/earthparams/" + model_fname)) {
        fname = path_ + "/earthparams/" + model_fname;
    }
    else if(fexists(path_ + "/earthparams/" + model_fname + ".dat")) {
        fname = path_ + "/earthparams/" + model_fname + ".dat";
    }
    else if(fexists(path_ + "/" + model_fname)) {
        fname = path_ + "/" + model_fname;
    }
    else if(fexists(path_ + "/" + model_fname + ".dat")) {
        fname = path_ + "/" + model_fname + ".dat";
    }
    else {
        throw(std::runtime_error("Cannot open earth model file!"));
    }

    std::ifstream in(fname.c_str());

    // if the earthmodel file doesn't exist, stop simulation
    if(in.fail()){
        throw(std::runtime_error("Failed to open " + fname + " Set correct EarthParamsPath."));
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
            ss << "Earth model uses undefined material " << medtype;
            throw(std::runtime_error(ss.str()));
        }

        EarthSector sector;
        sector.material_id = materials_.GetMaterialId(medtype);
        sector.level = level;
        sector.name = label;
        sector.geo = Sphere(radius, 0).create();
        level -= 1;
        if(nparams == 1) {
            ss >> param;
            sector.density = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>(param).create();
        }
        else {
            std::vector<double> params;
            for(int i=0; i<nparams; ++i) {
                ss >> param;
                params.push_back(param);
            }
            RadialAxis1D radial_ax;
            sector.density = DensityDistribution1D<RadialAxis1D,PolynomialDistribution1D>(radial_ax, params).create();
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
        EarthSector const & sector = sectors_[i];
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

    if(ice_cap_angle < 0 || ice_cap_angle > LeptonInjector::Constants::pi) {
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
                EarthSector & sector = sectors_[i];
                Sphere const * geo = dynamic_cast<Sphere const *>(sector.geo.get());
                sector.geo = Sphere(Placement(Vector3D(0,0,ice_offset), QFromZXZr(0,0,0)), geo->GetRadius()-ice_offset, 0).create();
                //geo->SetRadius(geo->GetRadius()-ice_offset);
                //geo->SetPosition(Vector3D(0,0,ice_offset));
            }
        }
    }
}

double EarthModel::GetTargetMass(LeptonInjector::Particle::ParticleType target) const {
    double molar_mass = materials_.GetMolarMass(target); // grams per mole
    return molar_mass * LeptonInjector::Constants::GeV_per_amu;
}
