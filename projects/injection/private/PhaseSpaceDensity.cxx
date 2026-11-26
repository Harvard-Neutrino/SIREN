#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/utilities/Errors.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace siren {
namespace injection {

using MType = PhaseSpaceMeasure::Type;

// ================================================================ //
//  Measure conversion (density Jacobians)                            //
// ================================================================ //

namespace {

[[noreturn]] void ThrowUnconvertible(
    const char * reason,
    PhaseSpaceMeasure const & from,
    PhaseSpaceMeasure const & to,
    PhaseSpaceTopology topology)
{
    std::ostringstream oss;
    oss << "ConvertDensity: " << reason << " (from "
        << PhaseSpaceMeasureName(from) << " to " << PhaseSpaceMeasureName(to)
        << " in " << PhaseSpaceTopologyName(topology) << " topology; "
        << "from factorization indices [spectator=" << from.spectator
        << ", pair_first=" << from.pair_first
        << ", pair_second=" << from.pair_second << "], "
        << "to factorization indices [spectator=" << to.spectator
        << ", pair_first=" << to.pair_first
        << ", pair_second=" << to.pair_second << "])"
        << " [siren-docs: errors#measure-compat]";
    throw siren::utilities::MeasureCompatibilityError(oss.str());
}

} // anonymous namespace

double ConvertDensity(
    double density,
    PhaseSpaceMeasure const & from,
    PhaseSpaceMeasure const & to,
    PhaseSpaceTopology topology,
    siren::dataclasses::InteractionRecord const & record)
{
    if (from == to) return density;
    // Zero is still a density with a convention. Do not let it silently pass
    // an unsupported conversion, but avoid demanding record kinematics for a
    // supported conversion whose result is necessarily zero.
    if (density == 0.0 &&
        PhaseSpaceDensityConvertible(topology, from, to)) {
        return density;
    }
    if (topology == PhaseSpaceTopology::Decay2Body &&
        from.type == MType::CosThetaRest &&
        (to.type == MType::SolidAngleRest || to.type == MType::SolidAngleLab)) {
        return ConvertDensity(density / (2.0 * M_PI),
            PhaseSpaceMeasure::SolidAngleRest(), to, topology, record);
    }
    namespace J = siren::injection::phase_space_jacobian;

    // ---- Decay2Body: SolidAngleRest <-> SolidAngleLab ----

    if (topology == PhaseSpaceTopology::Decay2Body &&
        from.type == MType::SolidAngleLab &&
        to.type == MType::SolidAngleLab) {
        double rest_density = ConvertDensity(
            density, from, PhaseSpaceMeasure::SolidAngleRest(),
            topology, record);
        return ConvertDensity(
            rest_density, PhaseSpaceMeasure::SolidAngleRest(), to,
            topology, record);
    }

    if (topology == PhaseSpaceTopology::Decay2Body &&
        ((from.type == MType::SolidAngleRest &&
          to.type == MType::SolidAngleLab) ||
         (from.type == MType::SolidAngleLab &&
          to.type == MType::SolidAngleRest)))
    {
        if (record.secondary_masses.size() < 2 ||
            record.secondary_momenta.size() < 2)
            ThrowUnconvertible(
                "Decay2Body rest<->lab requires 2 secondary masses and momenta",
                from, to, topology);

        PhaseSpaceMeasure const & lab_measure =
            from.type == MType::SolidAngleLab ? from : to;
        int daughter_index = lab_measure.spectator;
        if (daughter_index < 0 || daughter_index > 1)
            ThrowUnconvertible(
                "Decay2Body SolidAngleLab daughter index must be 0 or 1",
                from, to, topology);
        int other_index = 1 - daughter_index;

        double M = record.primary_mass;
        double m_A = record.secondary_masses[daughter_index];
        double m_B = record.secondary_masses[other_index];
        double p_rest = TwoBodyRestMomentum(M, m_A, m_B);
        double E_rest = TwoBodyRestEnergy(M, m_A, m_B);

        double E_parent = record.primary_momentum[0];
        double px = record.primary_momentum[1];
        double py = record.primary_momentum[2];
        double pz = record.primary_momentum[3];
        double p_parent = std::sqrt(px*px + py*py + pz*pz);
        // Parent at rest: lab and rest frames coincide.
        if (p_parent < 1e-15) return density;

        if (!(p_rest > 0.0))
            ThrowUnconvertible(
                "record's parent mass is at or below the two-body "
                "threshold for these daughter masses",
                from, to, topology);

        double beta_parent = p_parent / E_parent;
        double gamma_parent = E_parent / M;

        auto const & daughter_momentum =
            record.secondary_momenta[daughter_index];
        double px_A = daughter_momentum[1];
        double py_A = daughter_momentum[2];
        double pz_A = daughter_momentum[3];
        double p_A = std::sqrt(px_A*px_A + py_A*py_A + pz_A*pz_A);
        if (p_A < 1e-15)
            ThrowUnconvertible(
                "daughter momentum is degenerate (|p| ~ 0), so the lab "
                "angle is undefined for this record",
                from, to, topology);

        double cos_theta_lab =
            (px_A*px + py_A*py + pz_A*pz) / (p_A * p_parent);

        double E_A_lab = daughter_momentum[0];
        double p_par_lab = p_A * cos_theta_lab;
        double p_par_rest = gamma_parent * (p_par_lab - beta_parent * E_A_lab);
        double cos_theta_rest_actual = (p_rest > 0) ? p_par_rest / p_rest : 0.0;

        auto solutions = SolveLabAngle(
            beta_parent, gamma_parent, p_rest, E_rest,
            m_A, cos_theta_lab);

        double best_J = 0.0;
        double best_dist = 1e30;
        for (auto const & sol : solutions) {
            if (!sol.valid) continue;
            double dist = std::abs(sol.cos_theta_rest - cos_theta_rest_actual);
            if (dist < best_dist) {
                best_dist = dist;
                best_J = sol.jacobian;
            }
        }
        // An unconvertible point must not pass through unconverted.
        if (best_J <= 0.0 || !std::isfinite(best_J))
            ThrowUnconvertible(
                best_dist > 1e29
                    ? "no lab<->rest solution matches this record (lab "
                      "angle outside the kinematically allowed cone; "
                      "inconsistent record?)"
                    : "lab<->rest Jacobian is singular at this point "
                      "(record sits on the critical-angle fold)",
                from, to, topology);

        if (from.type == MType::SolidAngleRest) {
            return density / best_J;
        } else {
            return density * best_J;
        }
    }

    // ---- Scatter2to2 conversions ----

    if (topology == PhaseSpaceTopology::Scatter2to2) {
        using siren::dataclasses::MeasureHasExplicitAzimuth;
        using siren::dataclasses::MeasureIntegratesAzimuth;
        using siren::dataclasses::MeasureWithExplicitAzimuth;

        if (MeasureHasExplicitAzimuth(from) &&
            MeasureIntegratesAzimuth(to)) {
            ThrowUnconvertible(
                "cannot marginalize an explicit azimuth pointwise; evaluate "
                "the mixture in an explicit-azimuth measure",
                from, to, topology);
        }

        // Lift an azimuth-integrated physical density using its declared
        // uniform conditional p(phi)=1/(2*pi). Opaque/nonuniform omitted
        // azimuths must use Unspecified rather than one of these measures.
        if (MeasureIntegratesAzimuth(from) &&
            MeasureHasExplicitAzimuth(to)) {
            PhaseSpaceMeasure lifted = MeasureWithExplicitAzimuth(from);
            double joint_density = density / (2.0 * M_PI);
            if (lifted == to) return joint_density;
            return ConvertDensity(
                joint_density, lifted, to, topology, record);
        }

        // Bjorken (x,y) <-> (Q2,y). These measures have two non-azimuthal
        // coordinates; they are intentionally distinct from one-dimensional
        // MandelstamQ2.
        bool bjorken_to_q2y =
            (from.type == MType::BjorkenXY &&
             to.type == MType::MandelstamQ2Y) ||
            (from.type == MType::BjorkenXYPhi &&
             to.type == MType::MandelstamQ2YPhi);
        if (bjorken_to_q2y) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end())
                ThrowUnconvertible(
                    "requires 'bjorken_y' (or 'y') in interaction_parameters",
                    from, to, topology);
            return J::BjorkenXYDensityToQ2YDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }
        bool q2y_to_bjorken =
            (from.type == MType::MandelstamQ2Y &&
             to.type == MType::BjorkenXY) ||
            (from.type == MType::MandelstamQ2YPhi &&
             to.type == MType::BjorkenXYPhi);
        if (q2y_to_bjorken) {
            auto it = record.interaction_parameters.find("bjorken_y");
            if (it == record.interaction_parameters.end())
                it = record.interaction_parameters.find("y");
            if (it == record.interaction_parameters.end())
                ThrowUnconvertible(
                    "requires 'bjorken_y' (or 'y') in interaction_parameters",
                    from, to, topology);
            return J::Q2YDensityToBjorkenXYDensity(
                density, it->second, record.target_mass,
                record.primary_momentum[0]);
        }

        // MandelstamQ2 <-> fixed-mass one-dimensional y. Unlike Bjorken
        // (x,y), this Jacobian has no extra factor of y.
        bool y_to_q2 =
            (from.type == MType::FixedMassY &&
             to.type == MType::MandelstamQ2) ||
            (from.type == MType::FixedMassYPhi &&
             to.type == MType::MandelstamQ2Phi);
        if (y_to_q2) {
            return J::FixedMassYDensityToMandelstamQ2Density(
                density, record.target_mass, record.primary_momentum[0]);
        }
        bool q2_to_y =
            (from.type == MType::MandelstamQ2 &&
             to.type == MType::FixedMassY) ||
            (from.type == MType::MandelstamQ2Phi &&
             to.type == MType::FixedMassYPhi);
        if (q2_to_y) {
            return J::MandelstamQ2DensityToFixedMassYDensity(
                density, record.target_mass, record.primary_momentum[0]);
        }

        // SolidAngleRest <-> MandelstamQ2Phi retains azimuth and is therefore
        // an ordinary pointwise Jacobian with no 2*pi factor.
        if (from.type == MType::SolidAngleRest &&
            to.type == MType::MandelstamQ2Phi) {
            if (record.secondary_masses.size() < 2)
                ThrowUnconvertible(
                    "SolidAngleRest<->MandelstamQ2Phi requires 2 secondary masses",
                    from, to, topology);
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::SolidAngleRestDensityToMandelstamQ2PhiDensity(
                density, s, record.primary_mass, record.target_mass,
                record.secondary_masses[0], record.secondary_masses[1]);
        }
        if (from.type == MType::MandelstamQ2Phi &&
            to.type == MType::SolidAngleRest) {
            if (record.secondary_masses.size() < 2)
                ThrowUnconvertible(
                    "SolidAngleRest<->MandelstamQ2Phi requires 2 secondary masses",
                    from, to, topology);
            double E = record.primary_momentum[0];
            double s = record.primary_mass * record.primary_mass
                     + record.target_mass * record.target_mass
                     + 2.0 * record.target_mass * E;
            return J::MandelstamQ2PhiDensityToSolidAngleRestDensity(
                density, s, record.primary_mass, record.target_mass,
                record.secondary_masses[0], record.secondary_masses[1]);
        }

        // Compose the remaining fixed-mass measures through a Q2 coordinate
        // with the same azimuth treatment.
        auto fixed_marginal = [](MType type) {
            return type == MType::MandelstamQ2
                || type == MType::FixedMassY;
        };
        if (fixed_marginal(from.type) && fixed_marginal(to.type)) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2(), topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2(),
                to, topology, record);
        }
        auto fixed_joint = [](MType type) {
            return type == MType::SolidAngleRest
                || type == MType::MandelstamQ2Phi
                || type == MType::FixedMassYPhi;
        };
        if (fixed_joint(from.type) && fixed_joint(to.type)) {
            double d = ConvertDensity(density, from,
                PhaseSpaceMeasure::MandelstamQ2Phi(), topology, record);
            return ConvertDensity(d, PhaseSpaceMeasure::MandelstamQ2Phi(),
                to, topology, record);
        }

        if (from.type == MType::SolidAngleLab || to.type == MType::SolidAngleLab) {
            ThrowUnconvertible(
                "SolidAngleLab boost conversions are only implemented for "
                "Decay2Body topology (the rest<->lab Jacobian is the "
                "parent-rest-frame two-body boost; a topology-correct CM "
                "boost for scattering is not implemented)",
                from, to, topology);
        }
    }

    // ---- Decay3Body / Scatter2to3: 3-body measure conversions ----
    // Uses the measure's index fields to read the correct momenta.

    if (topology == PhaseSpaceTopology::Decay3Body ||
        topology == PhaseSpaceTopology::Scatter2to3) {

        // Different Dalitz pair choices are related by a linear change among
        // the three pair invariants, whose absolute Jacobian is one. Do not
        // route this through the recursive two-body Jacobians below.
        if (from.type == MType::DalitzPair &&
            to.type == MType::DalitzPair) {
            return density;
        }

        // For scattering, the parent mass in the Jacobian is the CM
        // energy sqrt(s), not the beam particle mass.
        double parent_mass = record.primary_mass;
        if (topology == PhaseSpaceTopology::Scatter2to3) {
            double E = record.primary_momentum[0];
            double m_beam = record.primary_mass;
            double m_target = record.target_mass;
            double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E;
            if (s > 0) parent_mass = std::sqrt(s);
        }

        auto compute_s_pair = [&](PhaseSpaceMeasure const & m) -> double {
            int i1 = m.pair_first;
            int i2 = m.pair_second;
            if (i1 < 0 || i2 < 0 ||
                i1 >= static_cast<int>(record.secondary_momenta.size()) ||
                i2 >= static_cast<int>(record.secondary_momenta.size())) {
                return -1.0;
            }
            auto const & p1 = record.secondary_momenta[i1];
            auto const & p2 = record.secondary_momenta[i2];
            double E = p1[0] + p2[0];
            double px = p1[1] + p2[1];
            double py = p1[2] + p2[2];
            double pz = p1[3] + p2[3];
            return E*E - px*px - py*py - pz*pz;
        };

        // Same recursive/helicity type, different indices:
        // cross-factorization conversion through Dalitz.
        if (from.type == to.type && from != to) {
            // Step 1: convert from's Recursive2Body to Dalitz using from's indices
            double s_pair_from = compute_s_pair(from);
            if (s_pair_from < 0)
                ThrowUnconvertible(
                    "invalid 'from' pair invariant mass (bad factorization indices)",
                    from, to, topology);
            double dalitz_density = J::Recursive2BodyDensityToDalitzDensity(
                density, parent_mass,
                record.secondary_masses[from.spectator],
                record.secondary_masses[from.pair_first],
                record.secondary_masses[from.pair_second], s_pair_from);

            // Step 2: convert Dalitz to to's Recursive2Body using to's indices
            double s_pair_to = compute_s_pair(to);
            if (s_pair_to < 0)
                ThrowUnconvertible(
                    "invalid 'to' pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::DalitzDensityToRecursive2BodyDensity(
                dalitz_density, parent_mass,
                record.secondary_masses[to.spectator],
                record.secondary_masses[to.pair_first],
                record.secondary_masses[to.pair_second], s_pair_to);
        }

        // Recursive2Body <-> HelicityAngles: trivial (same measure, different angular reference)
        if ((from.type == MType::Recursive2Body && to.type == MType::HelicityAngles) ||
            (from.type == MType::HelicityAngles && to.type == MType::Recursive2Body)) {
            return density;
        }

        // Recursive2Body <-> DalitzPair
        if (from.type == MType::Recursive2Body && to.type == MType::DalitzPair) {
            double s_pair = compute_s_pair(from);
            if (s_pair < 0)
                ThrowUnconvertible(
                    "invalid pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::Recursive2BodyDensityToDalitzDensity(
                density, parent_mass,
                record.secondary_masses[from.spectator],
                record.secondary_masses[from.pair_first],
                record.secondary_masses[from.pair_second], s_pair);
        }
        if (from.type == MType::DalitzPair && to.type == MType::Recursive2Body) {
            double s_pair = compute_s_pair(to);
            if (s_pair < 0)
                ThrowUnconvertible(
                    "invalid pair invariant mass (bad factorization indices)",
                    from, to, topology);
            return J::DalitzDensityToRecursive2BodyDensity(
                density, parent_mass,
                record.secondary_masses[to.spectator],
                record.secondary_masses[to.pair_first],
                record.secondary_masses[to.pair_second], s_pair);
        }

        // HelicityAngles <-> DalitzPair: compose through Recursive2Body
        if (from.type == MType::HelicityAngles && to.type == MType::DalitzPair) {
            PhaseSpaceMeasure intermediate = PhaseSpaceMeasure::Recursive2Body(
                from.spectator, from.pair_first, from.pair_second);
            return ConvertDensity(density, intermediate, to, topology, record);
        }
        if (from.type == MType::DalitzPair && to.type == MType::HelicityAngles) {
            PhaseSpaceMeasure intermediate = PhaseSpaceMeasure::Recursive2Body(
                to.spectator, to.pair_first, to.pair_second);
            return ConvertDensity(density, from, intermediate, topology, record);
        }
    }

    ThrowUnconvertible(
        "no conversion implemented for this measure pair and topology",
        from, to, topology);
}

} // namespace injection
} // namespace siren
