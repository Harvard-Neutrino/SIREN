#include "SIREN/injection/DetectorDirected2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"

#include "SIREN/injection/Isotropic2BodyChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/geometry/AABB.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace siren {
namespace injection {

static const double TWO_PI = 2.0 * M_PI;
static const double FOUR_PI = 4.0 * M_PI;

DetectorDirected2BodyChannel::DetectorDirected2BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int daughter_index)
    : target_(std::move(target))
    , daughter_index_(daughter_index)
{}

// Estimate the solid angle subtended by the target geometry as
// seen from a given position.
//
// Uses the bounding sphere of the AABB as a conservative
// estimate.  This overestimates the solid angle (samples some
// directions that miss), but the Density() method checks
// whether each direction actually hits, so the weight is correct.
double DetectorDirected2BodyChannel::TargetSolidAngle(
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double dx = extent.GetX();
    double dy = extent.GetY();
    double dz = extent.GetZ();
    double bounding_radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

    siren::math::Vector3D offset = center - position;
    double dist = offset.magnitude();

    if (dist < bounding_radius) {
        // Inside or near the bounding sphere — full sky
        return FOUR_PI;
    }

    // Solid angle of a sphere: 2*pi*(1 - cos(theta))
    // where sin(theta) = R / d
    double sin_theta = bounding_radius / dist;
    double cos_theta = std::sqrt(1.0 - sin_theta * sin_theta);
    return TWO_PI * (1.0 - cos_theta);
}

// Sample a direction uniformly over the solid angle subtended
// by the target as seen from position.
std::pair<siren::math::Vector3D, double>
DetectorDirected2BodyChannel::SampleTargetDirection(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    siren::math::Vector3D const & position) const
{
    auto aabb = target_->GetWorldBoundingBox();
    siren::math::Vector3D center = (aabb.min_corner + aabb.max_corner) * 0.5;
    siren::math::Vector3D extent = aabb.max_corner - aabb.min_corner;
    double dx = extent.GetX();
    double dy = extent.GetY();
    double dz = extent.GetZ();
    double bounding_radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

    siren::math::Vector3D to_center = center - position;
    double dist = to_center.magnitude();

    double solid_angle;
    siren::math::Vector3D direction;

    if (dist < bounding_radius) {
        // Inside bounding sphere — sample full sky
        double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
        double phi = TWO_PI * random->Uniform(0, 1);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        direction = siren::math::Vector3D(
            sin_theta * std::cos(phi),
            sin_theta * std::sin(phi),
            cos_theta
        );
        solid_angle = FOUR_PI;
    } else {
        // Sample uniformly within the bounding cone around to_center
        double sin_cone = bounding_radius / dist;
        double cos_cone = std::sqrt(1.0 - sin_cone * sin_cone);
        solid_angle = TWO_PI * (1.0 - cos_cone);

        // Sample cos(theta) uniform in [cos_cone, 1] (cone around axis)
        double cos_theta = cos_cone + (1.0 - cos_cone) * random->Uniform(0, 1);
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        double phi = TWO_PI * random->Uniform(0, 1);

        // Local coordinates: z along to_center direction
        siren::math::Vector3D axis = to_center;
        axis.normalize();

        // Build an orthonormal frame
        siren::math::Vector3D perp1, perp2;
        if (std::abs(axis.GetX()) < 0.9) {
            perp1 = siren::math::Vector3D(1, 0, 0);
        } else {
            perp1 = siren::math::Vector3D(0, 1, 0);
        }
        perp2 = siren::math::vector_product(axis, perp1);
        perp2.normalize();
        perp1 = siren::math::vector_product(perp2, axis);
        perp1.normalize();

        direction = axis * cos_theta
                  + perp1 * (sin_theta * std::cos(phi))
                  + perp2 * (sin_theta * std::sin(phi));
        direction.normalize();
    }

    return {direction, solid_angle};
}

bool DetectorDirected2BodyChannel::DirectionHitsTarget(
    siren::math::Vector3D const & position,
    siren::math::Vector3D const & direction) const
{
    auto intersections = target_->Intersections(position, direction);
    for (auto const & isec : intersections) {
        if (isec.distance > 0 && isec.entering) {
            return true;
        }
    }
    return false;
}

void DetectorDirected2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2) {
        throw std::runtime_error(
            "DetectorDirected2BodyChannel requires exactly 2 secondaries");
    }

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    // Parent kinematics
    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];
    double p_parent = std::sqrt(
        px_parent * px_parent + py_parent * py_parent + pz_parent * pz_parent);

    if (p_parent < 1e-15) {
        // Parent at rest — fall back to isotropic
        Isotropic2BodyChannel fallback(daughter_index_);
        fallback.Sample(random, nullptr, record);
        return;
    }

    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;

    // Parent direction unit vector
    siren::math::Vector3D parent_dir(
        px_parent / p_parent, py_parent / p_parent, pz_parent / p_parent);

    // Decay position (interaction vertex of this record)
    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Sample a lab-frame direction toward the target
    auto [lab_dir, solid_angle] = SampleTargetDirection(random, decay_pos);

    // cos(theta_lab) = angle between daughter A and parent direction
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);

    // Solve for rest-frame angle(s)
    auto solutions = SolveLabAngle(
        beta, gamma, p_rest, E_A_rest, m_A, cos_theta_lab);

    int n_valid = (solutions[0].valid ? 1 : 0) + (solutions[1].valid ? 1 : 0);

    if (n_valid == 0) {
        // This lab direction is kinematically forbidden.
        // Fall back to isotropic sampling.
        Isotropic2BodyChannel fallback(daughter_index_);
        fallback.Sample(random, nullptr, record);
        return;
    }

    // Pick a solution (when 2 exist, pick proportional to Jacobian)
    int chosen = 0;
    if (n_valid == 2) {
        double w0 = solutions[0].jacobian;
        double w1 = solutions[1].jacobian;
        double r = random->Uniform(0, 1) * (w0 + w1);
        chosen = (r < w0) ? 0 : 1;
    } else {
        chosen = solutions[0].valid ? 0 : 1;
    }

    double cos_theta_rest = solutions[chosen].cos_theta_rest;
    double p_lab = solutions[chosen].p_lab;

    // Sample azimuthal angle around parent direction
    double phi = TWO_PI * random->Uniform(0, 1);

    // But we already have the full lab direction from the target
    // sampling.  We need to construct rest-frame momentum consistent
    // with this lab direction.
    //
    // The lab direction determines cos_theta_lab w.r.t. parent.
    // The azimuthal angle phi_lab is already set by lab_dir.
    // We need to construct the 4-momentum in the lab frame.

    // Lab-frame 3-momentum of daughter A: magnitude p_lab, direction lab_dir
    double E_A_lab = std::sqrt(p_lab * p_lab + m_A * m_A);

    double px_A = p_lab * lab_dir.GetX();
    double py_A = p_lab * lab_dir.GetY();
    double pz_A = p_lab * lab_dir.GetZ();

    record.secondary_momenta[daughter_index_] = {E_A_lab, px_A, py_A, pz_A};

    // Daughter B by momentum conservation
    double E_B_lab = E_parent - E_A_lab;
    record.secondary_momenta[1 - daughter_index_] = {
        E_B_lab,
        px_parent - px_A,
        py_parent - py_A,
        pz_parent - pz_A
    };
}

double DetectorDirected2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 2) return 0.0;

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    // Parent kinematics
    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];
    double p_parent = std::sqrt(
        px_parent * px_parent + py_parent * py_parent + pz_parent * pz_parent);

    if (p_parent < 1e-15) return 1.0 / FOUR_PI;

    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;

    // Daughter A lab direction
    double px_A = record.secondary_momenta[daughter_index_][1];
    double py_A = record.secondary_momenta[daughter_index_][2];
    double pz_A = record.secondary_momenta[daughter_index_][3];
    double p_A = std::sqrt(px_A * px_A + py_A * py_A + pz_A * pz_A);
    if (p_A < 1e-15) return 0.0;

    siren::math::Vector3D lab_dir(px_A / p_A, py_A / p_A, pz_A / p_A);
    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Check if this direction hits the target
    if (!DirectionHitsTarget(decay_pos, lab_dir)) {
        return 0.0;
    }

    // Solid angle of the target as seen from the decay point
    double solid_angle = TargetSolidAngle(decay_pos);

    // Lab-frame angular density: uniform over the target solid angle
    // g(Omega_lab) = 1 / solid_angle
    double g_lab = 1.0 / solid_angle;

    // cos(theta_lab) between daughter and parent
    siren::math::Vector3D parent_dir(
        px_parent / p_parent, py_parent / p_parent, pz_parent / p_parent);
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);

    // Get the Jacobian |dOmega_lab / dOmega_rest|
    auto solutions = SolveLabAngle(
        beta, gamma, p_rest, E_A_rest, m_A, cos_theta_lab);

    // The generation density in rest-frame solid angle measure is:
    //   p_gen(Omega_rest) = g(Omega_lab) * |dOmega_lab/dOmega_rest|
    //
    // When there are two solutions, we selected with probability
    // proportional to the Jacobian, so the density for the chosen
    // solution is:
    //   p_gen = g(Omega_lab) * jacobian_chosen / (jacobian_0 + jacobian_1)
    //         * (jacobian_0 + jacobian_1)
    //         ... actually:
    //
    // The density at a point is the sum over all solutions that
    // produce this rest-frame angle, weighted by their mapping
    // density.  Since the lab angle maps to 1 or 2 rest-frame
    // angles, and we sample uniformly in lab direction, the density
    // at a rest-frame point Omega_rest that maps to lab direction
    // Omega_lab is:
    //
    //   p_gen(Omega_rest) = g(Omega_lab(Omega_rest)) * |dOmega_lab/dOmega_rest|
    //
    // This is the same formula regardless of how many solutions exist.
    // The multi-solution ambiguity is a property of the inverse map,
    // not the density.

    // Find which solution matches the actual daughter momentum
    double cos_theta_rest_actual;
    {
        // Boost daughter A back to rest frame to find actual cos_theta_rest
        // p_z_rest = gamma * (p_lab * cos_lab - beta * E_lab)
        double E_A_lab = record.secondary_momenta[daughter_index_][0];
        double p_par_lab = p_A * cos_theta_lab;
        double p_par_rest = gamma * (p_par_lab - beta * E_A_lab);
        cos_theta_rest_actual = p_par_rest / p_rest;
    }

    // Find the matching solution and its Jacobian
    for (auto const & sol : solutions) {
        if (sol.valid && std::abs(sol.cos_theta_rest - cos_theta_rest_actual) < 1e-6) {
            return g_lab * sol.jacobian;
        }
    }

    // No matching solution — shouldn't happen for a valid record
    return 0.0;
}

} // namespace injection
} // namespace siren
