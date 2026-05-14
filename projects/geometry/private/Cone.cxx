#include "SIREN/geometry/Cone.h"

#include <cmath>
#include <tuple>
#include <math.h>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "SIREN/math/Vector3D.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/Placement.h"

namespace siren {
namespace geometry {

Cone::Cone()
    : Geometry((std::string)("Cone"))
    , rmin1_(0.0)
    , rmax1_(0.0)
    , rmin2_(0.0)
    , rmax2_(0.0)
      , z_(0.0) {
    // Do nothing here
}

Cone::Cone(double rmin1, double rmax1, double rmin2, double rmax2, double z)
    : Geometry((std::string)("Cone"))
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
      , z_(z) {
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
}

Cone::Cone(Placement const & placement)
    : Geometry((std::string)("Cone"), placement)
    , rmin1_(0.0)
    , rmax1_(0.0)
    , rmin2_(0.0)
    , rmax2_(0.0)
      , z_(0.0) {
    // Do nothing here
}

Cone::Cone(Placement const & placement, double rmin1, double rmax1, double rmin2, double rmax2, double z)
    : Geometry((std::string)("Cone"), placement)
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
      , z_(z) {
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
}

Cone::Cone(const Cone& cone)
    : Geometry(cone)
    , rmin1_(cone.rmin1_)
    , rmax1_(cone.rmax1_)
    , rmin2_(cone.rmin2_)
    , rmax2_(cone.rmax2_)
      , z_(cone.z_) {
    // Nothing to do here
}

// ------------------------------------------------------------------------- //
void Cone::swap(Geometry& geometry) {
    Cone* cone = dynamic_cast<Cone*>(&geometry);
    if(!cone) {
        //log_warn("Cannot swap Cone!");
        return;
    }

    Geometry::swap(*cone);

    std::swap(rmin1_, cone->rmin1_);
    std::swap(rmax1_, cone->rmax1_);
    std::swap(rmin2_, cone->rmin2_);
    std::swap(rmax2_, cone->rmax2_);
    std::swap(z_, cone->z_);
}

//------------------------------------------------------------------------- //
Cone& Cone::operator=(const Geometry& geometry) {
    if(this != &geometry) {
        const Cone* cone = dynamic_cast<const Cone*>(&geometry);
        if(!cone) {
            //log_warn("Cannot assign Cone!");
            return *this;
        }
        Cone tmp(*cone);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Cone::equal(const Geometry& geometry) const
{
    const Cone* cone = dynamic_cast<const Cone*>(&geometry);

    if(!cone)
        return false;
    else if(rmin1_ != cone->rmin1_)
        return false;
    else if(rmax1_ != cone->rmax1_)
        return false;
    else if(rmin2_ != cone->rmin2_)
        return false;
    else if(rmax2_ != cone->rmax2_)
        return false;
    else if(z_ != cone->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Cone::less(const Geometry& geometry) const
{
    const Cone* cone = dynamic_cast<const Cone*>(&geometry);

    return
        std::tie(rmin1_, rmax1_, rmin2_, rmax2_, z_)
        <
        std::tie(cone->rmin1_, cone->rmax1_, cone->rmin2_, cone->rmax2_, cone->z_);
}

void Cone::print(std::ostream& os) const
{
    os << "Rmin1: " << rmin1_ << "\tRmax1: " << rmax1_
       << "\tRmin2: " << rmin2_ << "\tRmax2: " << rmax2_
       << "\tHeight: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Cone::ComputeIntersections(siren::math::Vector3D const & position, siren::math::Vector3D const & direction) const {
    // Calculate intersections of a ray with a conical section (frustum).
    // The cone outer surface at height z has radius:
    //   r_outer(z) = rmax1_ + (rmax2_ - rmax1_) * (z / z_ + 0.5)
    // which can be rewritten as r(z) = a + b*z where:
    //   a = (rmax1_ + rmax2_) / 2   (radius at z=0)
    //   b = (rmax2_ - rmax1_) / z_  (slope)
    //
    // Surface equation: x^2 + y^2 = (a + b*z)^2
    // Substituting ray p + t*d and expanding gives quadratic A*t^2 + B*t + C = 0

    double px = position.GetX();
    double py = position.GetY();
    double pz = position.GetZ();
    double dx = direction.GetX();
    double dy = direction.GetY();
    double dz = direction.GetZ();

    double z_calc_pos = 0.5 * z_;
    double z_calc_neg = -0.5 * z_;

    // Max 8 intersections: 2 outer barrel + 2 inner barrel + 2 top cap + 2 bottom cap
    Intersection hits[8];
    int n_hits = 0;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    // --- Outer conical surface ---
    {
        double a_outer = (rmax1_ + rmax2_) * 0.5;
        double b_outer = 0.0;
        if(z_ > 0) {
            b_outer = (rmax2_ - rmax1_) / z_;
        }

        // Quadratic coefficients for outer surface
        // A*t^2 + B*t + C = 0
        double A = dx * dx + dy * dy - b_outer * b_outer * dz * dz;
        double B = 2.0 * (px * dx + py * dy - b_outer * dz * (a_outer + b_outer * pz));
        double C = px * px + py * py - (a_outer + b_outer * pz) * (a_outer + b_outer * pz);

        if(!(dx == 0 && dy == 0 && b_outer == 0)) // not parallel to barrel
        {
            double determinant = B * B - 4.0 * A * C;

            if(std::fabs(A) > GEOMETRY_PRECISION) {
                // Standard quadratic
                if(determinant > 0) {
                    double sqrt_det = std::sqrt(determinant);
                    double t1 = (-B + sqrt_det) / (2.0 * A);
                    double t2 = (-B - sqrt_det) / (2.0 * A);

                    // Computer precision control
                    if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                        t1 = 0;
                    if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                        t2 = 0;

                    intersection_z = pz + t1 * dz;
                    if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                        intersection_x = px + t1 * dx;
                        intersection_y = py + t1 * dy;
                        // Check that the hit is on the correct nappe
                        // r(z) = a + b*z must be >= 0 at the intersection
                        double r_at_z = a_outer + b_outer * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                            hits[n_hits].distance = t1;
                            hits[n_hits].hierarchy = 0;
                            hits[n_hits].entering = entering;
                            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                            n_hits++;
                        }
                    }

                    intersection_z = pz + t2 * dz;
                    if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                        intersection_x = px + t2 * dx;
                        intersection_y = py + t2 * dy;
                        double r_at_z = a_outer + b_outer * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                            hits[n_hits].distance = t2;
                            hits[n_hits].hierarchy = 0;
                            hits[n_hits].entering = entering;
                            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                            n_hits++;
                        }
                    }
                }
            } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                // Linear case (A ~ 0): ray is parallel to cone surface
                double t1 = -C / B;
                if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;

                intersection_z = pz + t1 * dz;
                if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                    intersection_x = px + t1 * dx;
                    intersection_y = py + t1 * dy;
                    double r_at_z = a_outer + b_outer * intersection_z;
                    if(r_at_z >= 0) {
                        bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                        hits[n_hits].distance = t1;
                        hits[n_hits].hierarchy = 0;
                        hits[n_hits].entering = entering;
                        hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                        n_hits++;
                    }
                }
            }
        }
    }

    // --- Inner conical surface (hollow cone) ---
    if(rmin1_ > 0 || rmin2_ > 0) {
        double a_inner = (rmin1_ + rmin2_) * 0.5;
        double b_inner = 0.0;
        if(z_ > 0) {
            b_inner = (rmin2_ - rmin1_) / z_;
        }

        double A = dx * dx + dy * dy - b_inner * b_inner * dz * dz;
        double B = 2.0 * (px * dx + py * dy - b_inner * dz * (a_inner + b_inner * pz));
        double C = px * px + py * py - (a_inner + b_inner * pz) * (a_inner + b_inner * pz);

        if(!(dx == 0 && dy == 0 && b_inner == 0)) {
            double determinant = B * B - 4.0 * A * C;

            if(std::fabs(A) > GEOMETRY_PRECISION) {
                if(determinant > 0) {
                    double sqrt_det = std::sqrt(determinant);
                    double t1 = (-B + sqrt_det) / (2.0 * A);
                    double t2 = (-B - sqrt_det) / (2.0 * A);

                    if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                        t1 = 0;
                    if(t2 > 0 && t2 < GEOMETRY_PRECISION)
                        t2 = 0;

                    intersection_z = pz + t1 * dz;
                    if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                        intersection_x = px + t1 * dx;
                        intersection_y = py + t1 * dy;
                        double r_at_z = a_inner + b_inner * intersection_z;
                        if(r_at_z >= 0) {
                            // Inner surface: entering is inverted (entering inner = exiting volume)
                            bool entering = !((intersection_x * dx + intersection_y * dy - b_inner * r_at_z * dz) < 0);
                            hits[n_hits].distance = t1;
                            hits[n_hits].hierarchy = 0;
                            hits[n_hits].entering = entering;
                            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                            n_hits++;
                        }
                    }

                    intersection_z = pz + t2 * dz;
                    if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                        intersection_x = px + t2 * dx;
                        intersection_y = py + t2 * dy;
                        double r_at_z = a_inner + b_inner * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = !((intersection_x * dx + intersection_y * dy - b_inner * r_at_z * dz) < 0);
                            hits[n_hits].distance = t2;
                            hits[n_hits].hierarchy = 0;
                            hits[n_hits].entering = entering;
                            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                            n_hits++;
                        }
                    }
                }
            } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                double t1 = -C / B;
                if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;

                intersection_z = pz + t1 * dz;
                if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                    intersection_x = px + t1 * dx;
                    intersection_y = py + t1 * dy;
                    double r_at_z = a_inner + b_inner * intersection_z;
                    if(r_at_z >= 0) {
                        bool entering = !((intersection_x * dx + intersection_y * dy - b_inner * r_at_z * dz) < 0);
                        hits[n_hits].distance = t1;
                        hits[n_hits].hierarchy = 0;
                        hits[n_hits].entering = entering;
                        hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
                        n_hits++;
                    }
                }
            }
        }
    }

    // --- Top cap (z = +z_/2) annular disk ---
    if(dz != 0) {
        double t = (z_calc_pos - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = px + t * dx;
        intersection_y = py + t * dy;

        double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
        // At z = +z_/2, outer radius is rmax2_, inner radius is rmin2_
        if(r2_hit <= rmax2_ * rmax2_ && r2_hit >= rmin2_ * rmin2_) {
            intersection_z = pz + t * dz;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = dz < 0;
            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
            n_hits++;
        }
    }

    // --- Bottom cap (z = -z_/2) annular disk ---
    if(dz != 0) {
        double t = (z_calc_neg - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = px + t * dx;
        intersection_y = py + t * dy;

        double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
        // At z = -z_/2, outer radius is rmax1_, inner radius is rmin1_
        if(r2_hit <= rmax1_ * rmax1_ && r2_hit >= rmin1_ * rmin1_) {
            intersection_z = pz + t * dz;
            hits[n_hits].distance = t;
            hits[n_hits].hierarchy = 0;
            hits[n_hits].entering = dz > 0;
            hits[n_hits].position = siren::math::Vector3D(intersection_x, intersection_y, intersection_z);
            n_hits++;
        }
    }

    std::sort(hits, hits + n_hits, [](Intersection const & a, Intersection const & b) {
        return a.distance < b.distance;
    });
    return {hits, hits + n_hits};
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Cone::ComputeDistanceToBorder(const siren::math::Vector3D& position, const siren::math::Vector3D& direction) const {
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

    // No intersection with the cone
    if(dist.size() < 1) {
        distance.first  = -1;
        distance.second = -1;
    } else if(dist.size() == 1) // particle is inside the cone
    {
        distance.first  = dist.at(0);
        distance.second = -1;

    } else if(dist.size() == 2) // cone is infront of the particle
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);

        if(distance.second < distance.first) {
            std::swap(distance.first, distance.second);
        }

    } else {
        //log_error("This point should never be reached");
    }
    // Make a computer precision control!
    // This is necessary cause due to numerical effects it might happen
    // that a particle which is located on a geometry border is treated as
    // inside or outside

    if(distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if(distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if(distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

// ------------------------------------------------------------------------- //
AABB Cone::GetBoundingBox() const {
    double max_r = std::max(rmax1_, rmax2_);
    double hz = z_ * 0.5;
    return AABB(
        math::Vector3D(-max_r, -max_r, -hz),
        math::Vector3D( max_r,  max_r,  hz)
    );
}

} // namespace geometry
} // namespace siren
