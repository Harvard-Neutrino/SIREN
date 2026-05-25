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

FourVector ReadPrimary(siren::dataclasses::InteractionRecord const & record) {
    return FourVector{
        record.primary_momentum[0],
        siren::math::Vector3D(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3])
    };
}

FourVector ReadSecondary(siren::dataclasses::InteractionRecord const & record, int index) {
    return FourVector{
        record.secondary_momenta[index][0],
        siren::math::Vector3D(
            record.secondary_momenta[index][1],
            record.secondary_momenta[index][2],
            record.secondary_momenta[index][3])
    };
}

void WriteSecondary(
    siren::dataclasses::InteractionRecord & record,
    int index,
    FourVector const & p4)
{
    record.secondary_momenta[index] = {
        p4.e,
        p4.p.GetX(),
        p4.p.GetY(),
        p4.p.GetZ()
    };
}

double MassSquared(FourVector const & p4) {
    return p4.e * p4.e - p4.p * p4.p;
}

FourVector Add(FourVector const & a, FourVector const & b) {
    return FourVector{a.e + b.e, a.p + b.p};
}

bool DirectedTwoBody(
    FourVector const & parent,
    double parent_mass,
    double directed_mass,
    double other_mass,
    siren::math::Vector3D const & lab_dir,
    std::shared_ptr<siren::utilities::SIREN_random> random,
    FourVector & directed,
    FourVector & other,
    double & jacobian)
{
    double p_rest = TwoBodyRestMomentum(parent_mass, directed_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, directed_mass, other_mass);
    if (p_rest <= 0.0 || E_rest <= 0.0) return false;

    double p_parent = parent.p.magnitude();
    if (p_parent < 1e-15) {
        directed = FourVector{E_rest, lab_dir * p_rest};
        other = FourVector{
            parent.e - directed.e,
            parent.p - directed.p
        };
        jacobian = 1.0;
        return true;
    }

    siren::math::Vector3D parent_dir = parent.p / p_parent;
    double beta = p_parent / parent.e;
    double gamma = parent.e / parent_mass;
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
    cos_theta_lab = std::clamp(cos_theta_lab, -1.0, 1.0);

    auto solutions = SolveLabAngle(
        beta, gamma, p_rest, E_rest, directed_mass, cos_theta_lab);

    double total_j = 0.0;
    for (auto const & sol : solutions) {
        if (sol.valid) total_j += sol.jacobian;
    }
    if (total_j <= 0.0) return false;

    double r = random->Uniform(0, total_j);
    TwoBodyLabSolution chosen = solutions[0];
    double cumulative = 0.0;
    for (auto const & sol : solutions) {
        if (!sol.valid) continue;
        cumulative += sol.jacobian;
        if (r <= cumulative) {
            chosen = sol;
            break;
        }
    }

    double E_lab = std::sqrt(chosen.p_lab * chosen.p_lab
                           + directed_mass * directed_mass);
    directed = FourVector{E_lab, lab_dir * chosen.p_lab};
    other = FourVector{
        parent.e - directed.e,
        parent.p - directed.p
    };
    jacobian = chosen.jacobian;
    return true;
}

double DirectedTwoBodyJacobian(
    FourVector const & parent,
    double parent_mass,
    double directed_mass,
    double other_mass,
    FourVector const & directed)
{
    double p_rest = TwoBodyRestMomentum(parent_mass, directed_mass, other_mass);
    double E_rest = TwoBodyRestEnergy(parent_mass, directed_mass, other_mass);
    if (p_rest <= 0.0 || E_rest <= 0.0) return 0.0;

    double p_parent = parent.p.magnitude();
    double p_directed = directed.p.magnitude();
    if (p_directed <= 0.0) return 0.0;
    if (p_parent < 1e-15) return 1.0;

    siren::math::Vector3D parent_dir = parent.p / p_parent;
    siren::math::Vector3D lab_dir = directed.p / p_directed;
    double beta = p_parent / parent.e;
    double gamma = parent.e / parent_mass;
    double cos_theta_lab = siren::math::scalar_product(lab_dir, parent_dir);
    cos_theta_lab = std::clamp(cos_theta_lab, -1.0, 1.0);

    auto solutions = SolveLabAngle(
        beta, gamma, p_rest, E_rest, directed_mass, cos_theta_lab);

    double p_par_lab = p_directed * cos_theta_lab;
    double p_par_rest = gamma * (p_par_lab - beta * directed.e);
    double cos_theta_rest_actual = p_par_rest / p_rest;

    for (auto const & sol : solutions) {
        if (sol.valid && std::abs(sol.cos_theta_rest - cos_theta_rest_actual) < 1e-6) {
            return sol.jacobian;
        }
    }
    return 0.0;
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
    DetectorDirected2BodyChannel::Mode mode)
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

    double M = record.primary_mass;
    double m_spectator = record.secondary_masses[spectator_index_];
    double m_first = record.secondary_masses[pair_first_index_];
    double m_second = record.secondary_masses[pair_second_index_];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) {
        throw std::runtime_error("DetectorDirected3BodyChannel has no allowed invariant-mass range");
    }

    double s_pair = SampleInvariantMassSquared(random, s_min, s_max);
    double m_pair = std::sqrt(s_pair);

    FourVector parent = ReadPrimary(record);
    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    siren::math::Vector3D pair_dir = detail::SampleDirectedDirection(
        *target_, mode_, random, vertex);

    FourVector pair;
    FourVector spectator;
    double first_jacobian = 0.0;
    if (!DirectedTwoBody(parent, M, m_pair, m_spectator, pair_dir,
                         random, pair, spectator, first_jacobian)) {
        throw std::runtime_error("DetectorDirected3BodyChannel could not solve parent two-body step");
    }

    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_
        : pair_first_index_;

    double m_directed = record.secondary_masses[directed_pair_index_];
    double m_other = record.secondary_masses[other_pair_index];

    siren::math::Vector3D daughter_dir = detail::SampleDirectedDirection(
        *target_, mode_, random, vertex);

    FourVector directed;
    FourVector other;
    double second_jacobian = 0.0;
    if (!DirectedTwoBody(pair, m_pair, m_directed, m_other, daughter_dir,
                         random, directed, other, second_jacobian)) {
        throw std::runtime_error("DetectorDirected3BodyChannel could not solve pair two-body step");
    }

    WriteSecondary(record, spectator_index_, spectator);
    WriteSecondary(record, directed_pair_index_, directed);
    WriteSecondary(record, other_pair_index, other);
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

    double M = record.primary_mass;
    double m_spectator = record.secondary_masses[spectator_index_];
    double m_first = record.secondary_masses[pair_first_index_];
    double m_second = record.secondary_masses[pair_second_index_];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) return 0.0;

    FourVector first = ReadSecondary(record, pair_first_index_);
    FourVector second = ReadSecondary(record, pair_second_index_);
    FourVector pair = Add(first, second);
    double s_pair = MassSquared(pair);
    if (s_pair <= 0.0) return 0.0;
    double m_pair = std::sqrt(s_pair);

    FourVector parent = ReadPrimary(record);
    FourVector spectator = ReadSecondary(record, spectator_index_);
    FourVector directed = ReadSecondary(record, directed_pair_index_);
    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_
        : pair_first_index_;

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    double pair_p = pair.p.magnitude();
    double directed_p = directed.p.magnitude();
    if (pair_p <= 0.0 || directed_p <= 0.0) return 0.0;

    siren::math::Vector3D pair_dir = pair.p / pair_p;
    siren::math::Vector3D directed_dir = directed.p / directed_p;

    double g_pair = detail::AngularDensity(*target_, target_volume_, mode_, vertex, pair_dir);
    double g_directed = detail::AngularDensity(*target_, target_volume_, mode_, vertex, directed_dir);
    if (g_pair <= 0.0 || g_directed <= 0.0) return 0.0;

    double j_pair = DirectedTwoBodyJacobian(parent, M, m_pair, m_spectator, pair);
    double j_directed = DirectedTwoBodyJacobian(
        pair,
        m_pair,
        record.secondary_masses[directed_pair_index_],
        record.secondary_masses[other_pair_index],
        directed);
    if (j_pair <= 0.0 || j_directed <= 0.0) return 0.0;

    double s_density = InvariantMassDensity(s_pair, s_min, s_max);
    return s_density * g_pair * j_pair * g_directed * j_directed;
}

} // namespace injection
} // namespace siren
