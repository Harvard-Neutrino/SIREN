#include "SIREN/injection/PhaseSpaceChannel.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/PhaseSpaceJacobian.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/utilities/Errors.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace siren {
namespace injection {

using MType = PhaseSpaceMeasure::Type;

// ================================================================ //
//  Measure conversion (density Jacobians)                            //
// ================================================================ //

namespace {

int MeasurePriority(PhaseSpaceMeasure const & m) {
    switch (m.type) {
        case MType::SolidAngleRest:  return 0;
        case MType::Recursive2Body:  return 0;
        case MType::MandelstamQ2Phi:  return 1;
        case MType::FixedMassYPhi:    return 2;
        case MType::MandelstamQ2YPhi: return 3;
        case MType::BjorkenXYPhi:     return 4;
        case MType::MandelstamQ2:     return 5;
        case MType::FixedMassY:       return 6;
        case MType::MandelstamQ2Y:    return 7;
        case MType::BjorkenXY:        return 8;
        case MType::HelicityAngles:   return 9;
        case MType::DalitzPair:       return 10;
        case MType::SolidAngleLab:    return 11;
        case MType::Unspecified:      return 12;
    }
    return 12;
}

// Signal a conversion ConvertDensity cannot perform, instead of silently
// returning the input density unchanged.  The from==to
// short-circuit handles "no conversion needed"; this throw covers
// "conversion not implemented" and "required inputs missing", which a
// caller forming w = f/g must never silently absorb.
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

// ConvertDensity is a public free function (declared in PhaseSpaceChannel.h) so
// the same measure conversion the mixture uses internally is callable directly
// from Python.  The file-local ThrowUnconvertible / phase_space_jacobian helpers
// remain visible to it within this translation unit.
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
    namespace J = siren::injection::phase_space_jacobian;

    // ---- Decay2Body: SolidAngleRest <-> SolidAngleLab ----

    // Lab-angle measures for different daughters have different boost
    // Jacobians. Convert between them through the rest-frame solid angle,
    // where the two daughter directions differ only by an antipodal map with
    // unit Jacobian.
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

        // A pointwise density conversion cannot integrate a nonuniform
        // azimuth. Mixtures therefore elect an explicit-azimuth measure when
        // any channel supplies one, and only the supported uniform lift below
        // crosses from an integrated measure to an explicit one.
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

        // The rest<->lab Jacobian is the parent-rest-frame two-body boost,
        // which only Decay2Body reaches; no topology-correct CM boost is
        // implemented for scattering, so fail loud instead of silently
        // applying the wrong-frame two-body boost.
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

// ================================================================ //
//  MultiChannelPhaseSpace                                            //
// ================================================================ //

namespace {

// Guard the invariant every sampling/density path relies on: one weight per
// channel, each finite and non-negative, normalized to sum 1.  A mixture
// whose weights were mutated directly (public members) without a Normalize()
// call would otherwise sample from the wrong distribution (the cumulative
// loop silently dumps the leftover mass on the last channel) or mis-scale
// the density.  Throw instead.
void ThrowIfWeightsInconsistent(
    std::vector<std::shared_ptr<PhaseSpaceChannel>> const & channels,
    std::vector<double> const & weights,
    const char * where)
{
    if (channels.size() != weights.size()) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace::" << where
            << ": weights size (" << weights.size()
            << ") does not match channels size (" << channels.size() << ")"
            << " -- call Normalize() or use the validating constructor"
            << " [siren-docs: errors#configuration]";
        throw siren::utilities::ConfigurationError(oss.str());
    }
    for (size_t i = 0; i < weights.size(); ++i) {
        if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
            std::ostringstream oss;
            oss << "MultiChannelPhaseSpace::" << where
                << ": weight " << i << " is " << weights[i]
                << " (must be finite and non-negative)"
                << " -- call Normalize() or use the validating constructor"
                << " [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
    }
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (std::abs(sum - 1.0) > 1e-9) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace::" << where
            << ": weights sum to " << sum << ", not 1"
            << " -- call Normalize() or use the validating constructor"
            << " [siren-docs: errors#configuration]";
        throw siren::utilities::ConfigurationError(oss.str());
    }
}

} // anonymous namespace

MultiChannelPhaseSpace::MultiChannelPhaseSpace(
    std::vector<std::shared_ptr<PhaseSpaceChannel>> channels_,
    std::vector<double> weights_,
    bool allow_incompatible)
    : channels(std::move(channels_)),
      weights(std::move(weights_)),
      allow_incompatible_(allow_incompatible)
{
    Normalize();
    ThrowOnIncompatibility();
}

void MultiChannelPhaseSpace::Normalize() {
    if (channels.empty()) {
        throw siren::utilities::ConfigurationError(
            "MultiChannelPhaseSpace has no channels"
            " [siren-docs: errors#configuration]");
    }
    if (weights.empty()) {
        // Uniform prior: 1/N per channel.
        weights.assign(channels.size(), 1.0 / static_cast<double>(channels.size()));
        return;
    }
    if (weights.size() != channels.size()) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace::Normalize: weights size ("
            << weights.size() << ") does not match channels size ("
            << channels.size() << ") [siren-docs: errors#configuration]";
        throw siren::utilities::ConfigurationError(oss.str());
    }
    for (size_t i = 0; i < weights.size(); ++i) {
        if (!std::isfinite(weights[i]) || weights[i] < 0.0) {
            std::ostringstream oss;
            oss << "MultiChannelPhaseSpace::Normalize: weight " << i << " is "
                << weights[i] << " (must be finite and non-negative)"
                << " [siren-docs: errors#configuration]";
            throw siren::utilities::ConfigurationError(oss.str());
        }
    }
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum <= 0.0) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace::Normalize: weight sum is "
            << sum << " (must be > 0) [siren-docs: errors#configuration]";
        throw siren::utilities::ConfigurationError(oss.str());
    }
    for (double & w : weights) w /= sum;
}

std::size_t MultiChannelPhaseSpace::ConventionFingerprint() const {
    std::size_t seed = channels.size();
    auto combine = [&seed](std::size_t value) {
        seed ^= value + 0x9e3779b97f4a7c15ULL
              + (seed << 6) + (seed >> 2);
    };
    for (auto const & channel : channels) {
        combine(std::hash<PhaseSpaceChannel const *>{}(channel.get()));
        // Pointer identity alone is not enough: a destroyed channel's address
        // can be reused by its replacement, which would serve the old cached
        // measures. Fold in the convention values the cache exists to serve.
        combine(static_cast<std::size_t>(channel->Topology()));
        PhaseSpaceMeasure measure = channel->Measure();
        combine(static_cast<std::size_t>(measure.type));
        combine(static_cast<std::size_t>(measure.spectator));
        combine(static_cast<std::size_t>(measure.pair_first));
        combine(static_cast<std::size_t>(measure.pair_second));
        auto nested = std::dynamic_pointer_cast<NestedMixtureChannel>(channel);
        if (nested && nested->mixture) {
            combine(nested->mixture->ConventionFingerprint());
        }
    }
    return seed;
}

void MultiChannelPhaseSpace::EnsureConventionCache() const {
    std::size_t fingerprint = ConventionFingerprint();
    if (convention_cache_valid_ && fingerprint == convention_fingerprint_) {
        return;
    }

    convention_cache_valid_ = true;
    convention_fingerprint_ = fingerprint;
    cached_topology_error_.clear();
    cached_compatibility_diagnostics_.clear();
    cached_channel_measures_.clear();
    cached_common_topology_ = PhaseSpaceTopology::Unspecified;
    cached_common_measure_ = PhaseSpaceMeasure::Unspecified();
    if (channels.empty()) return;

    std::vector<PhaseSpaceTopology> topologies;
    std::vector<PhaseSpaceMeasure> measures;
    topologies.reserve(channels.size());
    measures.reserve(channels.size());
    for (auto const & channel : channels) {
        topologies.push_back(channel->Topology());
        measures.push_back(channel->Measure());
    }
    cached_channel_measures_ = measures;

    cached_common_topology_ = topologies.front();
    for (std::size_t i = 1; i < channels.size(); ++i) {
        if (topologies[i] == cached_common_topology_) continue;
        std::ostringstream validation;
        validation << "Topology mismatch: channel 0 (" << channels[0]->Name()
                   << ") is " << PhaseSpaceTopologyName(cached_common_topology_)
                   << " but channel " << i << " (" << channels[i]->Name()
                   << ") is " << PhaseSpaceTopologyName(topologies[i])
                   << " [siren-docs: errors#measure-compat]";
        cached_compatibility_diagnostics_.push_back(
            {ChannelDiagnostic::Severity::Fatal, validation.str()});
        if (cached_topology_error_.empty()) {
            std::ostringstream error;
            error << "MultiChannelPhaseSpace topology mismatch: channel 0 ("
                  << channels[0]->Name() << ") is "
                  << PhaseSpaceTopologyName(cached_common_topology_)
                  << " but channel " << i << " (" << channels[i]->Name()
                  << ") is " << PhaseSpaceTopologyName(topologies[i]);
            cached_topology_error_ = error.str();
        }
    }

    // Elect a measure that every channel can reach in the required direction.
    // In particular, one explicit-azimuth channel forces an explicit common
    // measure even if several marginal channels are present: those marginals
    // can be lifted uniformly, while the joint channel cannot be marginalized
    // pointwise. Among viable candidates, retain the majority/priority policy.
    bool found_viable = false;
    int best_count = -1;
    int best_priority = 100;
    bool has_specified = std::any_of(
        measures.begin(), measures.end(),
        [](PhaseSpaceMeasure const & measure) {
            return measure.type != MType::Unspecified;
        });
    for (auto const & candidate : measures) {
        if (has_specified && candidate.type == MType::Unspecified) continue;
        bool viable = std::all_of(
            measures.begin(), measures.end(),
            [&](PhaseSpaceMeasure const & source) {
                return PhaseSpaceDensityConvertible(
                    cached_common_topology_, source, candidate);
            });
        if (found_viable && !viable) continue;
        if (!found_viable && viable) {
            found_viable = true;
            best_count = -1;
            best_priority = 100;
        }
        int count = static_cast<int>(std::count(
            measures.begin(), measures.end(), candidate));
        int priority = MeasurePriority(candidate);
        if (count > best_count ||
            (count == best_count && priority < best_priority)) {
            cached_common_measure_ = candidate;
            best_count = count;
            best_priority = priority;
        }
    }

    // Preserve ValidateChannels' topology-first behavior.
    if (!cached_topology_error_.empty()) return;
    for (std::size_t i = 0; i < channels.size(); ++i) {
        if (measures[i] == cached_common_measure_) continue;
        std::ostringstream diagnostic;
        ChannelDiagnostic::Severity severity;
        if (!PhaseSpaceDensityConvertible(
                cached_common_topology_, measures[i],
                cached_common_measure_)) {
            diagnostic << "Measure incompatibility: channel " << i
                       << " (" << channels[i]->Name() << ") uses "
                       << PhaseSpaceMeasureName(measures[i])
                       << " which is not convertible to "
                       << PhaseSpaceMeasureName(cached_common_measure_)
                       << " within "
                       << PhaseSpaceTopologyName(cached_common_topology_)
                       << " topology"
                       << " [siren-docs: errors#measure-compat]";
            severity = ChannelDiagnostic::Severity::Fatal;
        } else {
            diagnostic << "Channel " << i << " (" << channels[i]->Name()
                       << ") uses " << PhaseSpaceMeasureName(measures[i])
                       << "; will auto-convert to "
                       << PhaseSpaceMeasureName(cached_common_measure_);
            severity = ChannelDiagnostic::Severity::Info;
        }
        cached_compatibility_diagnostics_.push_back(
            {severity, diagnostic.str()});
    }
}

int MultiChannelPhaseSpace::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    ThrowOnIncompatibility();

    if (channels.empty()) {
        throw siren::utilities::ConfigurationError(
            "MultiChannelPhaseSpace has no channels"
            " [siren-docs: errors#configuration]");
    }

    ThrowIfWeightsInconsistent(channels, weights, "Sample");

    double r = random->Uniform(0, 1);
    double cumulative = 0.0;
    int selected = static_cast<int>(channels.size()) - 1;
    for (int i = 0; i < static_cast<int>(channels.size()); ++i) {
        cumulative += weights[i];
        if (r < cumulative) {
            selected = i;
            break;
        }
    }

    channels[selected]->Sample(random, detector_model, record);
    return selected;
}

double MultiChannelPhaseSpace::ComputeContributions(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    std::vector<double> * weighted,
    std::vector<double> * bare,
    std::vector<bool> * directing_active) const
{
    phase_space_detail::ScopedEvaluationCache evaluation_scope;
    ThrowOnIncompatibility();

    if (channels.empty()) {
        throw siren::utilities::ConfigurationError(
            "MultiChannelPhaseSpace has no channels"
            " [siren-docs: errors#configuration]");
    }

    ThrowIfWeightsInconsistent(channels, weights, "ComputeContributions");

    if (!cached_topology_error_.empty()) {
        throw siren::utilities::MeasureCompatibilityError(cached_topology_error_);
    }
    PhaseSpaceTopology common_topo = cached_common_topology_;
    PhaseSpaceMeasure common_meas = cached_common_measure_;

    if (weighted) { weighted->clear(); weighted->reserve(channels.size()); }
    if (bare)     { bare->clear();     bare->reserve(channels.size()); }
    if (directing_active) {
        directing_active->clear();
        directing_active->reserve(channels.size());
    }

    double density = 0.0;
    for (size_t i = 0; i < channels.size(); ++i) {
        PhaseSpaceMeasure ch_meas = cached_channel_measures_[i];
        if (ch_meas != common_meas &&
            ch_meas.type == MType::Unspecified) {
            ThrowUnconvertible(
                "channel measure is unspecified; declare the density measure "
                "explicitly before mixing with specified channels",
                ch_meas, common_meas, common_topo);
        }

        double d;
        if (directing_active) {
            auto evaluation = channels[i]->Evaluate(detector_model, record);
            d = evaluation.density;
            directing_active->push_back(evaluation.directing_active);
        } else {
            d = channels[i]->Density(detector_model, record);
        }
        if (ch_meas != common_meas) {
            d = ConvertDensity(d, ch_meas, common_meas, common_topo, record);
        }
        if (bare) bare->push_back(d);
        double contribution = weights[i] * d;
        if (weighted) weighted->push_back(contribution);
        density += contribution;
    }
    return density;
}

double MultiChannelPhaseSpace::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    return ComputeContributions(detector_model, record, nullptr, nullptr);
}

double MultiChannelPhaseSpace::DensityIn(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceMeasure const & measure) const
{
    double density = ComputeContributions(
        detector_model, record, nullptr, nullptr);
    PhaseSpaceMeasure common = CommonMeasure();
    if (common == measure) return density;
    return ConvertDensity(density, common, measure,
                          CommonTopology(), record);
}

double MultiChannelPhaseSpace::DensityIn(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceConvention const & convention) const
{
    if (CommonTopology() != convention.topology) {
        std::ostringstream oss;
        oss << "Phase-space topology "
            << PhaseSpaceTopologyName(CommonTopology())
            << " does not match the requested density convention "
            << PhaseSpaceTopologyName(convention.topology)
            << " [siren-docs: errors#measure-compat]";
        throw siren::utilities::MeasureCompatibilityError(oss.str());
    }
    return DensityIn(detector_model, record, convention.measure);
}

std::vector<double> MultiChannelPhaseSpace::DensityBreakdown(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    std::vector<double> weighted;
    ComputeContributions(detector_model, record, &weighted, nullptr);
    return weighted;
}

void MultiChannelPhaseSpace::Accumulate(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    double weight, bool discount_fallback, bool recurse)
{
    std::vector<double> bare;
    std::vector<bool> activity;
    double g = ComputeContributions(
        detector_model, record, nullptr, &bare, &activity);
    if (g <= 0.0 || !std::isfinite(g)) return;
    PhaseSpaceTopology denominator_topology = cached_common_topology_;
    PhaseSpaceMeasure denominator_measure = cached_common_measure_;
    CreditAgainst(detector_model, record, weight * weight, g,
                  denominator_topology, denominator_measure,
                  discount_fallback, recurse, &bare, &activity);
}

void MultiChannelPhaseSpace::CreditAgainst(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record,
    double w2, double g_denom,
    PhaseSpaceTopology denominator_topology,
    PhaseSpaceMeasure denominator_measure,
    bool discount_fallback, bool recurse,
    std::vector<double> const * precomputed_bare,
    std::vector<bool> const * precomputed_activity)
{
    EnsureConventionCache();
    if (!cached_topology_error_.empty()) {
        throw siren::utilities::MeasureCompatibilityError(cached_topology_error_);
    }
    PhaseSpaceTopology local_topology = cached_common_topology_;
    PhaseSpaceMeasure local_measure = cached_common_measure_;
    std::vector<double> computed_bare;
    std::vector<bool> computed_activity;
    if (!precomputed_bare || !precomputed_activity) {
        ComputeContributions(
            detector_model, record, nullptr,
            &computed_bare, &computed_activity);
        precomputed_bare = &computed_bare;
        precomputed_activity = &computed_activity;
    }

    if (local_topology != denominator_topology) {
        throw std::runtime_error(
            "Nested Kleiss-Pittau credit topology does not match the outer "
            "density denominator");
    }
    if (local_measure != denominator_measure) {
        if (local_measure.type == MType::Unspecified) {
            ThrowUnconvertible(
                "nested mixture measure is unspecified; declare it before "
                "crediting against a specified outer density",
                local_measure, denominator_measure, denominator_topology);
        }
        computed_bare = *precomputed_bare;
        precomputed_bare = &computed_bare;
        for (double & density : computed_bare) {
            density = ConvertDensity(
                density, local_measure, denominator_measure,
                denominator_topology, record);
        }
    }

    if (kp_accumulator_.size() != channels.size())
        kp_accumulator_.assign(channels.size(), 0.0);
    kp_count_ += 1;

    for (size_t i = 0; i < channels.size(); ++i) {
        bool credit = !discount_fallback || (*precomputed_activity)[i];
        double density = (*precomputed_bare)[i];
        if (credit && density > 0.0 && std::isfinite(density)) {
            kp_accumulator_[i] += w2 * density / g_denom;
        }
        // Recurse into a nested sub-mixture regardless of the outer credit
        // decision: the inner channels are credited their own bare densities
        // against this same outer g (matching the Python inner loop).
        if (recurse) {
            auto nested = std::dynamic_pointer_cast<NestedMixtureChannel>(channels[i]);
            if (nested && nested->mixture) {
                nested->mixture->CreditAgainst(
                    detector_model, record, w2, g_denom,
                    denominator_topology, denominator_measure,
                    discount_fallback, recurse);
            }
        }
    }
}

void MultiChannelPhaseSpace::UpdateWeights(
    std::string const & update_rule, double damping, double min_weight,
    bool recurse, std::string const & failure_mode)
{
    if (update_rule != "sqrt_W" && update_rule != "alpha_sqrt_W") {
        throw std::invalid_argument(
            "MultiChannelPhaseSpace::UpdateWeights: unknown update_rule '"
            + update_rule + "' (expected 'sqrt_W' or 'alpha_sqrt_W')");
    }
    if (failure_mode != "ignore" && failure_mode != "coverage" &&
        failure_mode != "throughput") {
        throw std::invalid_argument(
            "MultiChannelPhaseSpace::UpdateWeights: unknown failure_mode '"
            + failure_mode + "' (expected 'ignore', 'coverage' or 'throughput')");
    }

    size_t n = channels.size();
    if (kp_count_ > 0 && n > 0 && weights.size() == n) {
        // Bare per-channel KP statistic W_i = <w^2 * g_i / g>.
        std::vector<double> W(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            double acc = (i < kp_accumulator_.size()) ? kp_accumulator_[i] : 0.0;
            W[i] = acc / static_cast<double>(kp_count_);
        }

        // Chain failure handling: adjust W_i by the channel's failed-selection
        // fraction f_i = fail/(succ+fail), only when failure-selection data is
        // present (so the single-vertex optimizer is unaffected):
        //   coverage   -- W_i /= (1 - f_i): up-weight lossy channels to keep
        //                 their region sampled (the original behavior).
        //   throughput -- W_i *= (1 - f_i): down-weight lossy channels so the
        //                 sampling tracks the successful contribution.
        //   ignore     -- no adjustment (the success-weighted statistic already
        //                 discounts failures implicitly).
        if (failure_mode != "ignore" && !kp_fail_select_.empty()) {
            for (size_t i = 0; i < n; ++i) {
                double succ = (i < kp_succ_select_.size()) ? kp_succ_select_[i] : 0.0;
                double fail = (i < kp_fail_select_.size()) ? kp_fail_select_[i] : 0.0;
                double tot = succ + fail;
                if (tot > 0.0) {
                    double f_i = fail / tot;
                    if (failure_mode == "coverage") {
                        if (f_i < 1.0) W[i] /= (1.0 - f_i);
                    } else {  // throughput
                        W[i] *= (1.0 - f_i);
                    }
                }
            }
        }

        // Candidate weights: sqrt(W_i) (memoryless) or alpha_i*sqrt(W_i)
        // (canonical).  Mirrors python _kp_update exactly.
        std::vector<double> cand(n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            double root = (W[i] > 0.0) ? std::sqrt(W[i]) : 0.0;
            cand[i] = (update_rule == "alpha_sqrt_W") ? weights[i] * root : root;
        }
        double total = std::accumulate(cand.begin(), cand.end(), 0.0);
        if (total > 0.0) {
            for (size_t i = 0; i < n; ++i) cand[i] /= total;
            // Floor + renormalize on the PRE-damping candidate (as in _kp_update);
            // damping is applied after, with no further floor/renormalize.
            if (min_weight > 0.0) {
                for (size_t i = 0; i < n; ++i) cand[i] = std::max(cand[i], min_weight);
                double t2 = std::accumulate(cand.begin(), cand.end(), 0.0);
                for (size_t i = 0; i < n; ++i) cand[i] /= t2;
            }
            for (size_t i = 0; i < n; ++i) {
                weights[i] = damping * cand[i] + (1.0 - damping) * weights[i];
            }
        }
        // total <= 0: degenerate batch -> keep weights (python returns None).
    }

    if (recurse) {
        for (auto const & channel : channels) {
            auto nested = std::dynamic_pointer_cast<NestedMixtureChannel>(channel);
            if (nested && nested->mixture) {
                nested->mixture->UpdateWeights(update_rule, damping, min_weight, recurse, failure_mode);
            }
        }
    }
    ResetAccumulators(false);
}

void MultiChannelPhaseSpace::ResetAccumulators(bool recurse)
{
    std::fill(kp_accumulator_.begin(), kp_accumulator_.end(), 0.0);
    kp_count_ = 0;
    std::fill(kp_succ_select_.begin(), kp_succ_select_.end(), 0.0);
    std::fill(kp_fail_select_.begin(), kp_fail_select_.end(), 0.0);
    if (recurse) {
        for (auto const & channel : channels) {
            auto nested = std::dynamic_pointer_cast<NestedMixtureChannel>(channel);
            if (nested && nested->mixture) {
                nested->mixture->ResetAccumulators(recurse);
            }
        }
    }
}

void MultiChannelPhaseSpace::AccumulateSelection(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record, bool failed)
{
    std::vector<double> weighted;
    double g = ComputeContributions(detector_model, record, &weighted, nullptr);
    if (g <= 0.0 || !std::isfinite(g)) return;

    if (kp_succ_select_.size() != channels.size())
        kp_succ_select_.assign(channels.size(), 0.0);
    if (kp_fail_select_.size() != channels.size())
        kp_fail_select_.assign(channels.size(), 0.0);

    std::vector<double> & sel = failed ? kp_fail_select_ : kp_succ_select_;
    for (size_t i = 0; i < channels.size(); ++i) {
        double p = weighted[i] / g;  // alpha_i * g_i / g
        if (weighted[i] > 0.0 && std::isfinite(p)) sel[i] += p;
    }
}

PhaseSpaceTopology MultiChannelPhaseSpace::CommonTopology() const {
    EnsureConventionCache();
    if (!cached_topology_error_.empty()) {
        throw std::runtime_error(cached_topology_error_);
    }
    return cached_common_topology_;
}

PhaseSpaceMeasure MultiChannelPhaseSpace::CommonMeasure() const {
    EnsureConventionCache();
    return cached_common_measure_;
}

PhaseSpaceConvention MultiChannelPhaseSpace::CommonConvention() const {
    return {CommonTopology(), CommonMeasure()};
}

std::vector<MultiChannelPhaseSpace::ChannelDiagnostic>
MultiChannelPhaseSpace::ValidateChannelsDetailed() const {
    EnsureConventionCache();
    return cached_compatibility_diagnostics_;
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannels() const {
    EnsureConventionCache();
    std::vector<std::string> messages;
    messages.reserve(cached_compatibility_diagnostics_.size());
    for (auto const & d : cached_compatibility_diagnostics_) {
        messages.push_back(d.message);
    }
    return messages;
}

void MultiChannelPhaseSpace::ThrowOnIncompatibility() const {
    // Populate the cache before the allow_incompatible_ opt-out:
    // ComputeContributions reads the cached conventions right after this
    // returns, so opting out of the throw must not leave the cache empty.
    EnsureConventionCache();
    if (allow_incompatible_) return;
    bool any_fatal = false;
    for (auto const & d : cached_compatibility_diagnostics_) {
        if (d.severity == ChannelDiagnostic::Severity::Fatal) {
            any_fatal = true;
            break;
        }
    }
    if (any_fatal) {
        std::ostringstream oss;
        oss << "MultiChannelPhaseSpace: incompatible channels:";
        for (auto const & d : cached_compatibility_diagnostics_) {
            if (d.severity == ChannelDiagnostic::Severity::Fatal) {
                oss << "\n  - " << d.message;
            }
        }
        throw siren::utilities::MeasureCompatibilityError(oss.str());
    }
}

std::vector<std::string> MultiChannelPhaseSpace::ValidateChannelDensities(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord template_record,
    int samples_per_channel) const
{
    std::vector<std::string> diagnostics;
    auto compat = ValidateChannels();
    diagnostics.insert(diagnostics.end(), compat.begin(), compat.end());

    for (size_t i = 0; i < channels.size(); ++i) {
        int n_zero = 0;
        int n_nan = 0;
        int n_neg = 0;
        int n_unsampleable = 0;
        double max_ratio = 0.0;

        for (int s = 0; s < samples_per_channel; ++s) {
            siren::dataclasses::InteractionRecord record = template_record;
            try {
                channels[i]->Sample(random, detector_model, record);
            } catch (siren::utilities::InjectionFailure const &) {
                // A retryable per-event failure on the synthetic template
                // means the template point lies outside this channel's
                // support (for example a tabulated mapping whose table does
                // not overlap the template's kinematic range), not a broken
                // configuration.  Skip the draw; configuration errors keep
                // propagating.
                ++n_unsampleable;
                continue;
            }

            double d_self = ConvertDensity(
                channels[i]->Density(detector_model, record),
                channels[i]->Measure(), CommonMeasure(),
                CommonTopology(), record);
            if (d_self <= 0 || !std::isfinite(d_self)) continue;

            for (size_t j = 0; j < channels.size(); ++j) {
                if (i == j) continue;
                double d_other = ConvertDensity(
                    channels[j]->Density(detector_model, record),
                    channels[j]->Measure(), CommonMeasure(),
                    CommonTopology(), record);

                if (std::isnan(d_other)) {
                    ++n_nan;
                } else if (d_other < 0) {
                    ++n_neg;
                } else if (d_other == 0) {
                    ++n_zero;
                } else {
                    double ratio = d_self / d_other;
                    if (ratio > max_ratio) max_ratio = ratio;
                }
            }
        }

        if (n_nan > 0) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_nan << "/" << samples_per_channel
                << " samples produced NaN density from other channels";
            diagnostics.push_back(oss.str());
        }
        if (n_neg > 0) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_neg << "/" << samples_per_channel
                << " samples produced negative density from other channels";
            diagnostics.push_back(oss.str());
        }
        if (n_zero > samples_per_channel / 2) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): " << n_zero << "/" << samples_per_channel
                << " samples got zero density from other channels "
                << "(possible measure mismatch or non-overlapping support)";
            diagnostics.push_back(oss.str());
        }
        if (max_ratio > 1e6) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): max density ratio " << max_ratio
                << " (extreme ratio may indicate measure mismatch)";
            diagnostics.push_back(oss.str());
        }
        if (n_unsampleable == samples_per_channel) {
            std::ostringstream oss;
            oss << "Channel " << i << " (" << channels[i]->Name()
                << "): the probe template lies outside the channel's "
                << "support (" << n_unsampleable << "/" << samples_per_channel
                << " draws raised InjectionFailure); density probe skipped "
                << "for this channel";
            diagnostics.push_back(oss.str());
        }
    }

    return diagnostics;
}

// ================================================================ //
//  NestedMixtureChannel                                              //
// ================================================================ //

void NestedMixtureChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    if (!mixture) {
        throw std::runtime_error("NestedMixtureChannel has no inner mixture");
    }
    mixture->Sample(random, detector_model, record);
}

double NestedMixtureChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!mixture) return 0.0;
    return mixture->Density(detector_model, record);
}

std::string NestedMixtureChannel::Name() const {
    return label;
}

PhaseSpaceTopology NestedMixtureChannel::Topology() const {
    if (!mixture) return PhaseSpaceTopology::Unspecified;
    return mixture->CommonTopology();
}

PhaseSpaceMeasure NestedMixtureChannel::Measure() const {
    if (!mixture) return PhaseSpaceMeasure::Unspecified();
    return mixture->CommonMeasure();
}

bool NestedMixtureChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const {
    // The group genuinely directs iff any member does.  A member with no
    // fallback notion (physical / isotropic) returns true by default, so a group
    // containing one is always active; a group of directed channels all in their
    // isotropic fallback returns false, letting the optimizer discount it as one.
    if (!mixture) return true;
    for (auto const & channel : mixture->channels) {
        if (channel && channel->DirectingActive(record)) return true;
    }
    return false;
}

} // namespace injection
} // namespace siren
