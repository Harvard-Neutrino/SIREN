#include "SIREN/geometry/ExtrPoly.h"

#include <tuple>
#include <math.h>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <functional>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace SI {
namespace geometry {

ExtrPoly::ExtrPoly()
    : Geometry((std::string)("ExtrPoly"))
      , polygon_({})
, zsections_({})
{
    // Do nothing here
}


ExtrPoly::ExtrPoly(const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry("ExtrPoly")
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(Placement const & placement)
    : Geometry((std::string)("ExtrPoly"), placement)
      , polygon_({})
, zsections_({})
{
    ComputeLateralPlanes();
}


ExtrPoly::ExtrPoly(Placement const & placement, const std::vector<std::vector<double>>& polygon,
        const std::vector<ExtrPoly::ZSection>& zsections)
    : Geometry((std::string)("ExtrPoly"), placement)
    , polygon_(polygon)
      , zsections_(zsections)
{
    if (polygon.size() < 3)
    {
        std::cout << "Need 3 polygon vertices at least!! Give it another shot";
        return;
    }
    ComputeLateralPlanes();
}

ExtrPoly::ExtrPoly(const ExtrPoly& extr)
    : Geometry(extr)
    , polygon_(extr.polygon_)
      , zsections_(extr.zsections_)
{
    ComputeLateralPlanes();
}


// ------------------------------------------------------------------------- //
void ExtrPoly::swap(Geometry& geometry)
{
    ExtrPoly* extr = dynamic_cast<ExtrPoly*>(&geometry);
    if (!extr)
    {
        //log_warn("Cannot swap ExtrPoly!");
        return;
    }

    Geometry::swap(*extr);

    std::swap(polygon_, extr->polygon_);
    std::swap(zsections_, extr->zsections_);
}

//------------------------------------------------------------------------- //
ExtrPoly& ExtrPoly::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);
        if (!extr)
        {
            //log_warn("Cannot assign ExtrPoly!");
            return *this;
        }

        ExtrPoly tmp(*extr);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::equal(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);

    if (!extr)
        return false;
    else if (polygon_ != extr->polygon_)
        return false;
    else if (!(zsections_ == extr->zsections_))
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool ExtrPoly::less(const Geometry& geometry) const
{
    const ExtrPoly* extr = dynamic_cast<const ExtrPoly*>(&geometry);

    return
        std::tie(polygon_, zsections_)
        <
        std::tie(extr->polygon_, extr->zsections_);
}

// ------------------------------------------------------------------------- //
void ExtrPoly::print(std::ostream& os) const
{
    os << "Not implemented yet\n";
    //os << "Polygon: " << polygon_ << "\tZsections: " << zsections_ << '\n';
}

// ------------------------------------------------------------------------- //
void ExtrPoly::ComputeLateralPlanes()
{
    int Nv = polygon_.size();
    planes_.resize(Nv);
    for (int i=0, k=Nv-1; i<Nv; k = i++)
    {
        std::vector<double> dir = {polygon_[i][0] - polygon_[k][0],polygon_[i][1] - polygon_[k][1]};
        double norm = sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
        dir[0]/=norm; dir[1]/=norm;
        planes_[i].a = -dir[1];
        planes_[i].b = dir[0];
        planes_[i].c = 0;
        planes_[i].d = dir[1]*polygon_[i][0] - dir[0]*polygon_[i][1];
    }
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> ExtrPoly::ComputeIntersections(SI::math::Vector3D const & position, SI::math::Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the extr poly
    // Implementation follows that of Geant4, see here:
    //
    // https://gitlab.cern.ch/geant4/geant4/-/blob/master/source/math/solids/specific/src/G4ExtrudedSolid.cc
    //
    // NOTE: Only works for convex right prisms at the moment

    std::vector<Geometry::Intersection> dist;

    SI::math::Vector3D intersection;

    std::function<void(double, bool)> save = [&](double t, bool entering){
        Intersection i;
        i.position = SI::math::Vector3D(position.GetX() + direction.GetX()*t,
                position.GetY() + direction.GetY()*t,
                position.GetZ() + direction.GetZ()*t);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };


    int Nz = zsections_.size();
    double z0 = zsections_[0].zpos;
    double z1 = zsections_[Nz-1].zpos;

    if ((position.GetZ() <= z0 + GEOMETRY_PRECISION) && direction.GetZ() <= 0) return dist;
    if ((position.GetZ() >= z1 - GEOMETRY_PRECISION) && direction.GetZ() >= 0) return dist;

    // Intersection with Z planes
    double dz = (z1 - z0)*0.5;
    double pz = position.GetZ() - dz - z0;

    double invz = (direction.GetZ() == 0) ? DBL_MAX : -1./direction.GetZ();
    double ddz = (invz < 0) ? dz : -dz;
    double tzmin = (pz + ddz)*invz;
    double tzmax = (pz - ddz)*invz;

    // Intersection with lateral planes
    int np = planes_.size();
    double txmin = tzmin, txmax = tzmax;
    for (int i=0; i<np; ++i)
    {
        double cosa = planes_[i].a*direction.GetX()+planes_[i].b*direction.GetY();
        double distnce = planes_[i].a*position.GetX()+planes_[i].b*position.GetY()+planes_[i].d;
        // case 1: particle is outside of outward-facing normal vector of plane in XY projection
        if (distnce >= -GEOMETRY_PRECISION)
        {
            if (cosa >= 0) { return dist; } // If particle is currently moving away from any face, it will never intersect
            double tmp  = -distnce/cosa;
            if (txmin < tmp)  { txmin = tmp; }

        }
        // case 2: particle is inside of outward-facing normal vector of plane
        else if (cosa > 0)
        {
            double tmp  = -distnce/cosa;
            if (txmax > tmp)  { txmax = tmp; }
        }
    }
    double tmin = txmin, tmax = txmax;
		if (tmax <= tmin + GEOMETRY_PRECISION)   // touch or no hit
      {
        return dist;
      }


    save(tmin,true);
    save(tmax,false);

    std::function<bool(Intersection const &, Intersection const &)> comp = [](Intersection const & a, Intersection const & b){
        return a.distance < b.distance;
    };

    std::sort(dist.begin(), dist.end(), comp);
    return dist;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> ExtrPoly::ComputeDistanceToBorder(const SI::math::Vector3D& position, const SI::math::Vector3D& direction) const
{
    // Compute the surface intersections
    std::vector<Intersection> intersections = Intersections(position, direction);
    std::vector<double> dist;
    bool first = true;
    for(unsigned int i=0; i<intersections.size(); ++i) {
        Intersection const & obj = intersections[i];
        if(obj.distance > 0) {
            if(first) {
                first = false;
                dist.push_back(obj.distance);
                if(not obj.entering) {
                    break;
                }
            }
            else {
                if(not obj.entering) {
                    dist.push_back(obj.distance);
                    break;
                }
                else {
                    throw(std::runtime_error("There should never be two \"entering\" intersections in a row!"));
                }
            }
        }
    }

    std::pair<double, double> distance;

    // No intersection with the outer cylinder
    if (dist.size() < 1)
    {
        distance.first  = -1;
        distance.second = -1;
        //    return distance;
    } else if (dist.size() == 1) // particle is inside the cylinder
    {
        distance.first  = dist.at(0);
        distance.second = -1;

    } else if (dist.size() == 2) // cylinder is infront of the particle
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);

        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else
    {
        //log_error("This point should never be reached");
    }
    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside

    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

} // namespace geometry
} // namespace SI
