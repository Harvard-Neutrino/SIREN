#include "SIREN/geometry/Cone.h"

#include <cmath>
#include <tuple>
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

namespace {

static const double TWO_PI = 2.0 * M_PI;

double NormalizePhi(double phi) {
    phi = std::fmod(phi, TWO_PI);
    if(phi < 0) phi += TWO_PI;
    return phi;
}

bool PhiInRange(double x, double y, double start_phi, double delta_phi) {
    double phi = NormalizePhi(std::atan2(y, x));
    double sp = NormalizePhi(start_phi);
    double ep = sp + delta_phi;
    if(ep <= TWO_PI + 1e-9) {
        return phi >= sp - 1e-9 && phi <= ep + 1e-9;
    } else {
        return phi >= sp - 1e-9 || phi <= NormalizePhi(ep) + 1e-9;
    }
}

} // anonymous namespace

Cone::Cone()
    : Geometry((std::string)("Cone"))
    , rmin1_(0.0)
    , rmax1_(0.0)
    , rmin2_(0.0)
    , rmax2_(0.0)
    , z_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
}

Cone::Cone(double rmin1, double rmax1, double rmin2, double rmax2, double z)
    : Geometry((std::string)("Cone"))
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(z_ <= 0) {
        throw std::invalid_argument("Cone height must be positive!");
    }
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
    if(rmax1_ <= 0 && rmax2_ <= 0) {
        throw std::invalid_argument("Cone must have at least one positive outer radius!");
    }
}

Cone::Cone(double rmin1, double rmax1, double rmin2, double rmax2, double z,
           double start_phi, double delta_phi)
    : Geometry((std::string)("Cone"))
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(z_ <= 0) {
        throw std::invalid_argument("Cone height must be positive!");
    }
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
    if(rmax1_ <= 0 && rmax2_ <= 0) {
        throw std::invalid_argument("Cone must have at least one positive outer radius!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Cone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Cone::Cone(Placement const & placement)
    : Geometry((std::string)("Cone"), placement)
    , rmin1_(0.0)
    , rmax1_(0.0)
    , rmin2_(0.0)
    , rmax2_(0.0)
    , z_(0.0)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
}

Cone::Cone(Placement const & placement, double rmin1, double rmax1, double rmin2, double rmax2, double z)
    : Geometry((std::string)("Cone"), placement)
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
    , z_(z)
    , start_phi_(0.0)
    , delta_phi_(2.0 * M_PI)
    , has_phi_cut_(false) {
    if(z_ <= 0) {
        throw std::invalid_argument("Cone height must be positive!");
    }
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
    if(rmax1_ <= 0 && rmax2_ <= 0) {
        throw std::invalid_argument("Cone must have at least one positive outer radius!");
    }
}

Cone::Cone(Placement const & placement, double rmin1, double rmax1, double rmin2, double rmax2, double z,
           double start_phi, double delta_phi)
    : Geometry((std::string)("Cone"), placement)
    , rmin1_(rmin1)
    , rmax1_(rmax1)
    , rmin2_(rmin2)
    , rmax2_(rmax2)
    , z_(z)
    , start_phi_(start_phi)
    , delta_phi_(delta_phi) {
    if(z_ <= 0) {
        throw std::invalid_argument("Cone height must be positive!");
    }
    if(rmin1_ > rmax1_) {
        std::swap(rmin1_, rmax1_);
    }
    if(rmin2_ > rmax2_) {
        std::swap(rmin2_, rmax2_);
    }
    if(rmax1_ <= 0 && rmax2_ <= 0) {
        throw std::invalid_argument("Cone must have at least one positive outer radius!");
    }
    if(delta_phi_ <= 0 || delta_phi_ > 2.0 * M_PI + 1e-9) {
        throw std::invalid_argument("Cone delta_phi must be in (0, 2*pi]!");
    }
    has_phi_cut_ = (delta_phi_ < 2.0 * M_PI - 1e-9);
}

Cone::Cone(const Cone& cone)
    : Geometry(cone)
    , rmin1_(cone.rmin1_)
    , rmax1_(cone.rmax1_)
    , rmin2_(cone.rmin2_)
    , rmax2_(cone.rmax2_)
    , z_(cone.z_)
    , start_phi_(cone.start_phi_)
    , delta_phi_(cone.delta_phi_)
    , has_phi_cut_(cone.has_phi_cut_) {
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
    std::swap(start_phi_, cone->start_phi_);
    std::swap(delta_phi_, cone->delta_phi_);
    std::swap(has_phi_cut_, cone->has_phi_cut_);
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
    else if(start_phi_ != cone->start_phi_)
        return false;
    else if(delta_phi_ != cone->delta_phi_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Cone::less(const Geometry& geometry) const
{
    const Cone* cone = dynamic_cast<const Cone*>(&geometry);
    if(!cone) return false;

    return
        std::tie(rmin1_, rmax1_, rmin2_, rmax2_, z_, start_phi_, delta_phi_)
        <
        std::tie(cone->rmin1_, cone->rmax1_, cone->rmin2_, cone->rmax2_, cone->z_, cone->start_phi_, cone->delta_phi_);
}

void Cone::print(std::ostream& os) const {
    os << "Rmin1: " << rmin1_ << "\tRmax1: " << rmax1_
       << "\tRmin2: " << rmin2_ << "\tRmax2: " << rmax2_
       << "\tHeight: " << z_;
    if(has_phi_cut_) os << "\tStartPhi: " << start_phi_ << "\tDeltaPhi: " << delta_phi_;
    os << '\n';
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

    struct TaggedHit {
        double distance;
        siren::math::Vector3D position;
        bool entering;
        int source; // 0 = surface, 1 = wedge
    };

    TaggedHit all_hits[12]; // max: 2 outer + 2 inner + 2 top + 2 bottom + 2 wedge + spare
    int n_all = 0;

    double intersection_x;
    double intersection_y;
    double intersection_z;

    // --- Outer conical surface ---
    bool has_outer_apex_top = (rmax2_ <= 0 && rmin2_ <= 0);
    bool has_outer_apex_bot = (rmax1_ <= 0 && rmin1_ <= 0);
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

        int n_before_outer = n_all;
        bool outer_apex_handled = false;

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

                    // Detect near-coincident roots at the apex
                    if(has_outer_apex_top || has_outer_apex_bot) {
                        double apex_z = has_outer_apex_top ? z_calc_pos : z_calc_neg;
                        double z1 = pz + t1 * dz;
                        double z2 = pz + t2 * dz;
                        double apex_tol = GEOMETRY_PRECISION * 1e3;
                        if(std::fabs(z1 - apex_z) < apex_tol && std::fabs(z2 - apex_z) < apex_tol) {
                            outer_apex_handled = true;
                            if(C < -GEOMETRY_PRECISION) {
                                double t_mid = 0.5 * (t1 + t2);
                                if(t_mid > 0 && t_mid < GEOMETRY_PRECISION) t_mid = 0;
                                // The cone body lies below a top apex and
                                // above a bottom apex, so a ray crossing
                                // the tip is entering the solid when it
                                // moves toward that side.
                                bool apex_entering = has_outer_apex_top ? (dz < 0.0) : (dz > 0.0);
                                all_hits[n_all] = {t_mid, siren::math::Vector3D(px + t_mid*dx, py + t_mid*dy, apex_z), apex_entering, 0};
                                n_all++;
                            }
                        }
                    }

                    if(!outer_apex_handled) {
                        intersection_z = pz + t1 * dz;
                        bool z1_ok = (has_outer_apex_bot ? intersection_z >= z_calc_neg : intersection_z > z_calc_neg)
                                  && (has_outer_apex_top ? intersection_z <= z_calc_pos : intersection_z < z_calc_pos);
                        if(z1_ok) {
                            intersection_x = px + t1 * dx;
                            intersection_y = py + t1 * dy;
                            double r_at_z = a_outer + b_outer * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                                all_hits[n_all] = {t1, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                                n_all++;
                            }
                        }

                        intersection_z = pz + t2 * dz;
                        bool z2_ok = (has_outer_apex_bot ? intersection_z >= z_calc_neg : intersection_z > z_calc_neg)
                                  && (has_outer_apex_top ? intersection_z <= z_calc_pos : intersection_z < z_calc_pos);
                        if(z2_ok) {
                            intersection_x = px + t2 * dx;
                            intersection_y = py + t2 * dy;
                            double r_at_z = a_outer + b_outer * intersection_z;
                            if(r_at_z >= 0) {
                                bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                                all_hits[n_all] = {t2, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                                n_all++;
                            }
                        }
                    }
                }
            } else if(std::fabs(B) > GEOMETRY_PRECISION) {
                // Linear case (A ~ 0): ray is parallel to cone surface
                double t1 = -C / B;
                if(t1 > 0 && t1 < GEOMETRY_PRECISION)
                    t1 = 0;

                intersection_z = pz + t1 * dz;
                bool z_ok = (has_outer_apex_bot ? intersection_z >= z_calc_neg : intersection_z > z_calc_neg)
                         && (has_outer_apex_top ? intersection_z <= z_calc_pos : intersection_z < z_calc_pos);
                if(z_ok) {
                    intersection_x = px + t1 * dx;
                    intersection_y = py + t1 * dy;
                    double r_at_z = a_outer + b_outer * intersection_z;
                    if(r_at_z >= 0) {
                        bool entering = (intersection_x * dx + intersection_y * dy - b_outer * r_at_z * dz) < 0;
                        all_hits[n_all] = {t1, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                        n_all++;
                    }
                }
            }
        }

        // Apex fallback: if no outer surface hit was produced at the apex,
        // check if the ray passes through the apex point
        if(!outer_apex_handled && (has_outer_apex_top || has_outer_apex_bot) && std::fabs(dz) > GEOMETRY_PRECISION) {
            bool found_apex_hit = false;
            double apex_z = has_outer_apex_top ? z_calc_pos : z_calc_neg;
            for(int i = n_before_outer; i < n_all; ++i) {
                if(std::fabs(all_hits[i].position.GetZ() - apex_z) < GEOMETRY_PRECISION * 1e3) {
                    found_apex_hit = true;
                    break;
                }
            }
            if(!found_apex_hit) {
                double t_apex = (apex_z - pz) / dz;
                double ix = px + t_apex * dx;
                double iy = py + t_apex * dy;
                double r2_apex = ix*ix + iy*iy;
                if(r2_apex < GEOMETRY_PRECISION * GEOMETRY_PRECISION * 1e6) {
                    if(C <= GEOMETRY_PRECISION) {
                        if(t_apex > 0 && t_apex < GEOMETRY_PRECISION) t_apex = 0;
                        // Entering when the ray moves toward the body
                        // side of the tip (below a top apex, above a
                        // bottom apex).
                        bool apex_entering = has_outer_apex_top ? (dz < 0.0) : (dz > 0.0);
                        all_hits[n_all] = {t_apex, siren::math::Vector3D(ix, iy, apex_z), apex_entering, 0};
                        n_all++;
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
                            all_hits[n_all] = {t1, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                            n_all++;
                        }
                    }

                    intersection_z = pz + t2 * dz;
                    if(intersection_z > z_calc_neg && intersection_z < z_calc_pos) {
                        intersection_x = px + t2 * dx;
                        intersection_y = py + t2 * dy;
                        double r_at_z = a_inner + b_inner * intersection_z;
                        if(r_at_z >= 0) {
                            bool entering = !((intersection_x * dx + intersection_y * dy - b_inner * r_at_z * dz) < 0);
                            all_hits[n_all] = {t2, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                            n_all++;
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
                        all_hits[n_all] = {t1, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), entering, 0};
                        n_all++;
                    }
                }
            }
        }
    }

    // --- Top cap (z = +z_/2) annular disk ---
    // Skip degenerate cap at apex (rmax=0): apex is handled by the surface solver
    if(std::fabs(dz) > GEOMETRY_PRECISION && !has_outer_apex_top) {
        double t = (z_calc_pos - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = px + t * dx;
        intersection_y = py + t * dy;

        double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
        // At z = +z_/2, outer radius is rmax2_, inner radius is rmin2_
        if(r2_hit <= rmax2_ * rmax2_ && r2_hit >= rmin2_ * rmin2_) {
            intersection_z = pz + t * dz;
            all_hits[n_all] = {t, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), dz < 0, 0};
            n_all++;
        }
    }

    // --- Bottom cap (z = -z_/2) annular disk ---
    // Skip degenerate cap at apex (rmax=0): apex is handled by the surface solver
    if(std::fabs(dz) > GEOMETRY_PRECISION && !has_outer_apex_bot) {
        double t = (z_calc_neg - pz) / dz;

        if(t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        intersection_x = px + t * dx;
        intersection_y = py + t * dy;

        double r2_hit = intersection_x * intersection_x + intersection_y * intersection_y;
        // At z = -z_/2, outer radius is rmax1_, inner radius is rmin1_
        if(r2_hit <= rmax1_ * rmax1_ && r2_hit >= rmin1_ * rmin1_) {
            intersection_z = pz + t * dz;
            all_hits[n_all] = {t, siren::math::Vector3D(intersection_x, intersection_y, intersection_z), dz > 0, 0};
            n_all++;
        }
    }

    if(!has_phi_cut_) {
        // No phi cut: surface hits are the final result
        if(n_all == 0) return {};
        std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
            return a.distance < b.distance;
        });
        std::vector<Intersection> result;
        result.reserve(n_all);
        for(int i = 0; i < n_all; ++i) {
            Intersection isect;
            isect.distance = all_hits[i].distance;
            isect.hierarchy = 0;
            isect.entering = all_hits[i].entering;
            isect.position = all_hits[i].position;
            result.push_back(isect);
        }
        return result;
    }

    // Phi cut: compute infinite wedge intersections (two half-planes from z-axis).
    // See Polycone.cxx for method description; same pattern in all phi-cut shapes.
    for(int face = 0; face < 2; ++face) {
        double alpha = start_phi_ + face * delta_phi_;
        double ca = std::cos(alpha), sa = std::sin(alpha);
        double nx, ny;
        if(face == 0) { nx = sa; ny = -ca; }
        else { nx = -sa; ny = ca; }
        double n_dot_d = nx*dx + ny*dy;
        if(std::fabs(n_dot_d) < GEOMETRY_PRECISION) continue;
        double n_dot_p = nx*px + ny*py;
        double t = -n_dot_p / n_dot_d;
        if(t > 0 && t < GEOMETRY_PRECISION) t = 0;
        double hx = px + t*dx, hy = py + t*dy, hz = pz + t*dz;
        if(hx*ca + hy*sa < -GEOMETRY_PRECISION) continue;
        bool entering = (n_dot_d < 0);
        all_hits[n_all] = {t, siren::math::Vector3D(hx, hy, hz), entering, 1};
        n_all++;
    }

    if(n_all == 0) return {};

    std::sort(all_hits, all_hits + n_all, [](TaggedHit const & a, TaggedHit const & b) {
        return a.distance < b.distance;
    });

    // Determine initial in_wedge state
    bool in_surface = false;
    bool in_wedge = false;
    bool has_wedge_hit = false;
    for(int i = 0; i < n_all; ++i) {
        if(all_hits[i].source == 1) {
            in_wedge = !all_hits[i].entering;
            has_wedge_hit = true;
            break;
        }
    }
    if(!has_wedge_hit) {
        in_wedge = PhiInRange(px, py, start_phi_, delta_phi_);
    }
    bool was_inside = in_surface && in_wedge;

    // CSG intersection walk
    std::vector<Intersection> result;
    for(int i = 0; i < n_all; ++i) {
        if(all_hits[i].source == 0) {
            in_surface = all_hits[i].entering;
        } else {
            in_wedge = all_hits[i].entering;
        }
        bool now_inside = in_surface && in_wedge;
        if(now_inside != was_inside) {
            Intersection isect;
            isect.distance = all_hits[i].distance;
            isect.hierarchy = 0;
            isect.entering = now_inside;
            isect.position = all_hits[i].position;
            result.push_back(isect);
        }
        was_inside = now_inside;
    }
    return result;
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
