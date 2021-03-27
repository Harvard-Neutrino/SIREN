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
}

void EarthModel::LoadDefaultMaterials() {
    materials_.AddMaterial("VACUUM", std::map<int, double>({{1000070080,1.0},})); // Assume there are 1 neutrons for every 7 protons in the universe
}

void EarthModel::LoadDefaultSectors() {
    EarthSector sector;
    sector.material_id = materials_.GetMaterialId("VACUUM");
    sector.level = std::numeric_limits<int>::min();
    sector.geo = Sphere(Vector3D(0,0,0), std::numeric_limits<double>::infinity(), 0, sector.level).create();
    sector.density = DensityDistribution1D<RadialAxis1D,ConstantDistribution1D>().create(); // Use the universe_mean_density from GEANT4
    AddSector(sector);
}

void EarthModel::LoadMaterialModel(std::string const & material_model) {
    materials_.SetPath(path_);
    materials_.AddModelFile(material_model);
}

double EarthModel::GetColumnDepthInCGS(Vector3D const & p0, Vector3D const & p1) const {
    Vector3D direction = p1 - p0;
    double distance = direction.magnitude();
    direction.normalize();
    std::vector<Geometry::Intersection> intersections;

    for(auto const & sector : sectors_) {
        std::vector<Geometry::Intersection> i = sector.geo->Intersections(p0, direction);
        intersections.reserve(intersections.size() + std::distance(i.begin(), i.end()));
        intersections.insert(intersections.end(), i.begin(), i.end());
    }

    assert(intersections.size() > 0);
    assert(distance > 0);

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

    double column_depth = 0;

    std::map<unsigned int, Geometry::Intersection const *> stack;
    Geometry::Intersection const * current_intersection = &intersections[0];
    stack.insert({intersections[0].hierarchy, current_intersection});
    for(unsigned int i=1; i<intersections.size(); ++i) {
        Geometry::Intersection const & intersection = intersections[i];
        if(intersection.entering) {
            stack.insert({intersection.hierarchy, &intersection});
            if(intersection.hierarchy > current_intersection->hierarchy) {
                if(intersection.distance > 0) {
                    // Store integral between current_intersection and new intersection
                    double end_point = std::min(intersection.distance, distance);
                    double start_point = std::max(current_intersection->distance, 0.0);
                    double segment_length = end_point - start_point;
                    double integral = GetSector(current_intersection->hierarchy).density->Integral(p0+start_point*direction, direction, segment_length);
                    column_depth += integral;
                    if(intersection.distance >= distance) {
                        break;
                    }
                }
                current_intersection = &intersection;
            }
        }
        else {
            if(intersection.hierarchy <= current_intersection->hierarchy) {
                stack.erase(intersection.hierarchy);
                if(intersection.distance > 0) {
                    // Store integral between current_intersection and new intersection
                    double end_point = std::min(intersection.distance, distance);
                    double start_point = std::max(current_intersection->distance, 0.0);
                    double segment_length = end_point - start_point;
                    double integral = GetSector(current_intersection->hierarchy).density->Integral(p0+start_point*direction, direction, segment_length);
                    column_depth += integral;
                    if(intersection.distance >= distance) {
                        break;
                    }
                }
                current_intersection = stack.lower_bound(intersection.hierarchy)->second;
            }
            else {
                throw("Cannot exit a level that we have not entered!");
            }
        }
    }
    return column_depth;
}

double EarthModel::DistanceForColumnDepthToPoint(Vector3D const & end_point, Vector3D const & direction, double column_depth, bool use_electron_density) const {

}

Vector3D EarthModel::GetEarthCoordPosFromDetCoordPos(Vector3D const & point) const {
    point + detector_origin_;
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
        sector.geo = Sphere(Vector3D(0,0,0), radius, 0, level).create();
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
