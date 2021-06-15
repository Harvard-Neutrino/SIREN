#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "earthmodel-service/EarthModel.h"
#include "earthmodel-service/Vector3D.h"

using namespace earthmodel;

namespace {
void string_to_lower(std::string & data) {
    std::transform(data.begin(), data.end(), data.begin(), [](unsigned char c){ return std::tolower(c); });
}
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
        throw("Already have a sector of that heirarchy!");
    }
    else {
        sector_map_[sector.level] = sectors_.size();
        sectors_.push_back(sector);
    }
}

EarthSector EarthModel::GetSector(int heirarchy) const {
    return sectors_[sector_map_.find(heirarchy)->second];
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
        throw("Received empty earth model filename!");

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
        throw("Cannot open earth model file!");
    }

    std::ifstream in(fname.c_str());

    // if the earthmodel file doesn't exist, stop simulation
    if(in.fail()){
        throw("Failed to open " + fname + " Set correct EarthParamsPath.");
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

            std::string shape;
            ss >> shape;
            ss >> xc >> yc >> zc;

            if(shape.find("sphere")!=std::string::npos) {
                double radius; // For Sphere shapes
                ss >> radius;
                sector.geo = Sphere(Vector3D(xc, yc, zc), radius, 0).create();
            }
            else if(shape.find("box")!=std::string::npos) {
                double dx, dy, dz; // For Box shapes
                ss >> dx >> dy >> dz;
                sector.geo = Box(Vector3D(xc, yc, zc), dx, dy, dz).create();
            }
            else {
                std::stringstream ss_err;
                ss_err
                    << "Shape \""
                    << shape
                    << "\" not recognized on line:\n"
                    << ss.str();
                throw(ss_err.str());
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
                throw(ss_err.str());
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
                throw(ss_err.str());
            }

            sectors_.push_back(sector);
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
    materials_.AddMaterial("VACUUM", std::map<int, double>({{1000070080,1.0},})); // Assume there are 1 neutrons for every 7 protons in the universe
}

void EarthModel::LoadDefaultSectors() {
    EarthSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = std::numeric_limits<int>::min();
    sector.geo = Sphere(Vector3D(0,0,0), std::numeric_limits<double>::infinity(), 0).create();
    sector.density = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>().create(); // Use the universe_mean_density from GEANT4
    AddSector(sector);
}

void EarthModel::LoadMaterialModel(std::string const & material_model) {
    materials_.SetPath(path_);
    materials_.AddModelFile(material_model);
}

double EarthModel::GetDensity(Geometry::IntersectionList const & intersections, Vector3D const & p0, bool use_electron_density) const {
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
    double density = -1.0;

    std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<Geometry::Intersection>::const_iterator current_intersection, std::vector<Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            EarthSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            if(use_electron_density)
                density *= materials_.GetPNERatio(sector.material_id);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(density >= 0);

    return density;
}

double EarthModel::GetDensity(Vector3D const & p0, bool use_electron_density) const {
    Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetDensity(intersections, p0, use_electron_density);
}

double EarthModel::GetColumnDepthInCGS(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & p1, bool use_electron_density) const {
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
            if(use_electron_density)
                integral *= materials_.GetPNERatio(sector.material_id);
            column_depth += integral;
        }
        // last_point = end_point;
        bool done = offset + dot * intersection->distance >= distance;
        return done;
    };

    SectorLoop(callback, intersections, dot < 0);

    return column_depth;
}

double EarthModel::GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1, bool use_electron_density) const {
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
    return GetColumnDepthInCGS(intersections, p0, p1, use_electron_density);
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
    int hierarchy = intersections.intersections[0].hierarchy;
    for(unsigned int i=1; i<intersections.intersections.size(); ++i) {
        if(intersections.intersections[i].hierarchy != hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
        }
    }
    for(unsigned int i=intersections.intersections.size()-1; i>=0; --i) {
        if(intersections.intersections[i].hierarchy != hierarchy) {
            result.intersections.push_back(intersections.intersections[i]);
        }
    }
    return result;
}

Geometry::IntersectionList EarthModel::GetOuterBounds(Vector3D const & p0, Vector3D const & direction) {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetOuterBounds(intersections);
}

void EarthModel::SectorLoop(std::function<bool(std::vector<Geometry::Intersection>::const_iterator, std::vector<Geometry::Intersection>::const_iterator, double)> callback, Geometry::IntersectionList const & intersections, bool reverse) const {
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
            throw("Reached end of intersections! Should never reach this point!");
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
                throw("Cannot exit a level that we have not entered!");
            }
        }
    }
}

double EarthModel::DistanceForColumnDepthFromPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & dir, double column_depth, bool use_electron_density) const {
    Vector3D direction = dir;
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
    double total_distance = -1;
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
            if(use_electron_density)
                target /= materials_.GetPNERatio(sector.material_id);
            double distance = sector.density->InverseIntegral(p0+start_point*direction, direction, target, segment_length);
            done = distance >= 0;
            double integral = sector.density->Integral(p0+start_point*direction, direction, segment_length);
            if(use_electron_density)
                integral *= materials_.GetPNERatio(sector.material_id);
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

double EarthModel::DistanceForColumnDepthFromPoint(Vector3D const & p0, Vector3D const & direction, double column_depth, bool use_electron_density) const {
    Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return DistanceForColumnDepthFromPoint(intersections, p0, direction, column_depth, use_electron_density);
}

double EarthModel::DistanceForColumnDepthToPoint(Geometry::IntersectionList const & intersections, Vector3D const & p0, Vector3D const & direction, double column_depth, bool use_electron_density) const {
    return DistanceForColumnDepthFromPoint(intersections, p0, -direction, column_depth, use_electron_density);
}

double EarthModel::DistanceForColumnDepthToPoint(Vector3D const & p0, Vector3D const & direction, double column_depth, bool use_electron_density) const {
    return DistanceForColumnDepthFromPoint(p0, -direction, column_depth, use_electron_density);
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
        throw("Received empty earth model filename!");

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
        throw("Cannot open earth model file!");
    }

    std::ifstream in(fname.c_str());

    // if the earthmodel file doesn't exist, stop simulation
    if(in.fail()){
        throw("Failed to open " + fname + " Set correct EarthParamsPath.");
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
            throw(ss.str());
        }

        EarthSector sector;
        sector.material_id = materials_.GetMaterialId(medtype);
        sector.level = level;
        sector.name = label;
        sector.geo = Sphere(Vector3D(0,0,0), radius, 0).create();
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
            throw("Layers must be radially ordered in file!");
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
                sector.geo = Sphere(Vector3D(0,0,ice_offset), geo->GetRadius()-ice_offset, 0).create();
                //geo->SetRadius(geo->GetRadius()-ice_offset);
                //geo->SetPosition(Vector3D(0,0,ice_offset));
            }
        }
    }
}
