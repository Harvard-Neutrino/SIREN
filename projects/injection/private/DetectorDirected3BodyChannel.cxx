#include "SIREN/injection/DetectorDirected3BodyChannel.h"

#include "DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

namespace {

struct FourVector {
    double e;
    siren::math::Vector3D p;
};

FourVector ReadSecondary(siren::dataclasses::InteractionRecord const & record, int index) {
    return FourVector{
        record.secondary_momenta[index][0],
        siren::math::Vector3D(
            record.secondary_momenta[index][1],
            record.secondary_momenta[index][2],
            record.secondary_momenta[index][3])
    };
}

double MassSquared(FourVector const & p4) {
    return p4.e * p4.e - p4.p * p4.p;
}

FourVector Add(FourVector const & a, FourVector const & b) {
    return FourVector{a.e + b.e, a.p + b.p};
}

// Compute the initial-state invariant mass and total 4-momentum.
// For decays: M = primary_mass, P = primary_momentum.
// For scattering (target at rest): M = sqrt(s), P = beam + target.
struct InitialState {
    double M;
    double E, px, py, pz;
};

InitialState GetInitialState(
    siren::dataclasses::InteractionRecord const & record,
    PhaseSpaceTopology topology)
{
    double E = record.primary_momentum[0];
    double px = record.primary_momentum[1];
    double py = record.primary_momentum[2];
    double pz = record.primary_momentum[3];

    if (topology == PhaseSpaceTopology::Scatter2to3) {
        // Target at rest: add target 4-momentum (m_target, 0, 0, 0)
        double m_beam = record.primary_mass;
        double m_target = record.target_mass;
        double E_total = E + m_target;
        double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E;
        if (s <= 0.0) return {0, 0, 0, 0, 0};
        return {std::sqrt(s), E_total, px, py, pz};
    }

    // Decay: parent 4-momentum is the initial state
    return {record.primary_mass, E, px, py, pz};
}

} // anonymous namespace

DetectorDirected3BodyChannel::DetectorDirected3BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int spectator_index,
    int pair_first_index,
    int pair_second_index,
    int directed_pair_index,
    InvariantMassMode mass_mode,
    double resonance_mass,
    double resonance_width,
    double power_law_nu,
    double power_law_offset,
    DetectorDirected2BodyChannel::Mode mode,
    PhaseSpaceTopology topology)
    : target_(std::move(target))
    , spectator_index_(spectator_index)
    , pair_first_index_(pair_first_index)
    , pair_second_index_(pair_second_index)
    , directed_pair_index_(directed_pair_index)
    , mass_mode_(mass_mode)
    , resonance_mass_(resonance_mass)
    , resonance_width_(resonance_width)
    , power_law_nu_(power_law_nu)
    , power_law_offset_(power_law_offset)
    , mode_(mode)
    , topology_(topology)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires a target geometry");
    }
    if (directed_pair_index_ != pair_first_index_ &&
        directed_pair_index_ != pair_second_index_) {
        throw std::runtime_error("directed_pair_index must be one of the pair indices");
    }
    double aabb_volume = detail::AABBVolume(*target_);
    target_volume_ = detail::ComputeGeometryVolume(*target_, aabb_volume);

    if (mode_ == DetectorDirected2BodyChannel::Mode::Volume) {
        if (aabb_volume <= 0.0) {
            throw std::runtime_error("Target bounding box has zero volume");
        }
        if (target_volume_ / aabb_volume < 1e-4) {
            throw std::runtime_error("Target volume is too small relative to its bounding box for sampling to be viable");
        }
    }
}

void DetectorDirected3BodyChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

double DetectorDirected3BodyChannel::SampleInvariantMassSquared(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    double s_min,
    double s_max) const
{
    double r = random->Uniform(0, 1);
    if (mass_mode_ == InvariantMassMode::BreitWigner) {
        if (resonance_mass_ <= 0.0 || resonance_width_ <= 0.0) {
            throw std::runtime_error("BreitWigner invariant-mass mapping requires positive mass and width");
        }
        BreitWignerMapping map(resonance_mass_, resonance_width_, s_min, s_max);
        return map.Forward(r);
    }
    if (mass_mode_ == InvariantMassMode::PowerLaw) {
        PowerLawMapping map(power_law_nu_, power_law_offset_, s_min, s_max);
        return map.Forward(r);
    }
    UniformMapping map(s_min, s_max);
    return map.Forward(r);
}

double DetectorDirected3BodyChannel::InvariantMassDensity(
    double s,
    double s_min,
    double s_max) const
{
    if (s < s_min || s > s_max || s_max <= s_min) return 0.0;
    if (mass_mode_ == InvariantMassMode::BreitWigner) {
        if (resonance_mass_ <= 0.0 || resonance_width_ <= 0.0) return 0.0;
        BreitWignerMapping map(resonance_mass_, resonance_width_, s_min, s_max);
        return map.Density(s);
    }
    if (mass_mode_ == InvariantMassMode::PowerLaw) {
        if (s_min <= power_law_offset_) return 0.0;
        PowerLawMapping map(power_law_nu_, power_law_offset_, s_min, s_max);
        return map.Density(s);
    }
    UniformMapping map(s_min, s_max);
    return map.Density(s);
}

void DetectorDirected3BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 3 ||
        record.secondary_masses.size() != 3 ||
        record.secondary_momenta.size() != 3) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires exactly 3 secondaries");
    }

    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_spectator = record.secondary_masses[spectator_index_];
    double m_first = record.secondary_masses[pair_first_index_];
    double m_second = record.secondary_masses[pair_second_index_];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) {
        throw std::runtime_error("DetectorDirected3BodyChannel has no allowed invariant-mass range");
    }

    // Step 0: sample invariant mass of the pair system
    double s_pair = SampleInvariantMassSquared(random, s_min, s_max);
    double m_pair = std::sqrt(s_pair);

    double E_parent = init.E;
    double px_parent = init.px;
    double py_parent = init.py;
    double pz_parent = init.pz;
    double p_parent = std::sqrt(px_parent*px_parent + py_parent*py_parent + pz_parent*pz_parent);

    // Step 1: P -> pair + spectator (isotropic in rest frame)
    // The pair direction is not biased — it is just a bookkeeping
    // device for factorizing the 3-body phase space.
    double p1_rest = TwoBodyRestMomentum(M, m_pair, m_spectator);
    double E1_rest = TwoBodyRestEnergy(M, m_pair, m_spectator);

    double cos1 = 2.0 * random->Uniform(0, 1) - 1.0;
    double phi1 = 2.0 * M_PI * random->Uniform(0, 1);
    double sin1 = std::sqrt(1.0 - cos1 * cos1);

    double pair_E, pair_px, pair_py, pair_pz;
    if (p_parent < 1e-15) {
        // Parent at rest: lab = rest
        pair_px = p1_rest * sin1 * std::cos(phi1);
        pair_py = p1_rest * sin1 * std::sin(phi1);
        pair_pz = p1_rest * cos1;
        pair_E = E1_rest;
    } else {
        double beta = p_parent / E_parent;
        double gamma = E_parent / M;
        siren::math::Vector3D parent_dir(px_parent/p_parent, py_parent/p_parent, pz_parent/p_parent);

        siren::math::Vector3D perp1, perp2;
        if (std::abs(parent_dir.GetX()) < 0.9)
            perp1 = siren::math::Vector3D(1, 0, 0);
        else
            perp1 = siren::math::Vector3D(0, 1, 0);
        perp2 = siren::math::vector_product(parent_dir, perp1);
        perp2.normalize();
        perp1 = siren::math::vector_product(perp2, parent_dir);
        perp1.normalize();

        double p_rest_par = p1_rest * cos1;
        double p_rest_perp1 = p1_rest * sin1 * std::cos(phi1);
        double p_rest_perp2 = p1_rest * sin1 * std::sin(phi1);

        pair_E = gamma * (E1_rest + beta * p_rest_par);
        double p_lab_par = gamma * (p_rest_par + beta * E1_rest);

        siren::math::Vector3D p_vec =
            parent_dir * p_lab_par + perp1 * p_rest_perp1 + perp2 * p_rest_perp2;
        pair_px = p_vec.GetX();
        pair_py = p_vec.GetY();
        pair_pz = p_vec.GetZ();
    }

    double spec_E = E_parent - pair_E;
    double spec_px = px_parent - pair_px;
    double spec_py = py_parent - pair_py;
    double spec_pz = pz_parent - pair_pz;

    // Step 2: pair -> directed + other (directed toward target)
    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_ : pair_first_index_;
    double m_directed = record.secondary_masses[directed_pair_index_];
    double m_other = record.secondary_masses[other_pair_index];

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    auto step2 = detail::SampleDirectedStep(
        pair_E, pair_px, pair_py, pair_pz,
        m_pair, m_directed, m_other,
        vertex, *target_, target_volume_, mode_, random);

    double dir_E = step2.E_lab;
    double dir_px = step2.p_lab * step2.lab_dir.GetX();
    double dir_py = step2.p_lab * step2.lab_dir.GetY();
    double dir_pz = step2.p_lab * step2.lab_dir.GetZ();

    double oth_E = pair_E - dir_E;
    double oth_px = pair_px - dir_px;
    double oth_py = pair_py - dir_py;
    double oth_pz = pair_pz - dir_pz;

    // Write results
    record.secondary_momenta[spectator_index_] = {spec_E, spec_px, spec_py, spec_pz};
    record.secondary_momenta[directed_pair_index_] = {dir_E, dir_px, dir_py, dir_pz};
    record.secondary_momenta[other_pair_index] = {oth_E, oth_px, oth_py, oth_pz};
    record.interaction_parameters["phase_space_s_pair"] = s_pair;
}

double DetectorDirected3BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 3 ||
        record.secondary_masses.size() != 3 ||
        record.secondary_momenta.size() != 3) {
        return 0.0;
    }

    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_spectator = record.secondary_masses[spectator_index_];
    double m_first = record.secondary_masses[pair_first_index_];
    double m_second = record.secondary_masses[pair_second_index_];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) return 0.0;

    // Reconstruct pair 4-vector
    FourVector first = ReadSecondary(record, pair_first_index_);
    FourVector second = ReadSecondary(record, pair_second_index_);
    FourVector pair = Add(first, second);
    double s_pair = MassSquared(pair);
    if (s_pair <= 0.0) return 0.0;
    double m_pair = std::sqrt(s_pair);

    // Invariant mass density
    double s_density = InvariantMassDensity(s_pair, s_min, s_max);
    if (s_density <= 0.0) return 0.0;

    // Step 1 density: P -> pair + spectator (isotropic)
    static const double INV_FOUR_PI = 1.0 / (4.0 * M_PI);
    double step1_density = INV_FOUR_PI;

    // Step 2 density: pair -> directed + other (directed toward target)
    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_ : pair_first_index_;
    FourVector directed = ReadSecondary(record, directed_pair_index_);

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    double step2_density = detail::DensityDirectedStep(
        pair.e, pair.p.GetX(), pair.p.GetY(), pair.p.GetZ(),
        m_pair,
        record.secondary_masses[directed_pair_index_],
        record.secondary_masses[other_pair_index],
        directed.e, directed.p.GetX(), directed.p.GetY(), directed.p.GetZ(),
        vertex, *target_, target_volume_, mode_);
    if (step2_density <= 0.0) return 0.0;

    return s_density * step1_density * step2_density;
}

} // namespace injection
} // namespace siren
