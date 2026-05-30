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
        double m_beam = record.primary_mass;
        double m_target = record.target_mass;
        double E_total = E + m_target;
        double s = m_beam * m_beam + m_target * m_target + 2.0 * m_target * E;
        if (s <= 0.0) return {0, 0, 0, 0, 0};
        return {std::sqrt(s), E_total, px, py, pz};
    }

    return {record.primary_mass, E, px, py, pz};
}

// Boost 4-momentum from parent's rest frame to lab frame.
void BoostRestToLab(
    double E_parent, double px_parent, double py_parent, double pz_parent,
    double M_parent,
    double E_rest, double px_rest, double py_rest, double pz_rest,
    double & E_lab, double & px_lab, double & py_lab, double & pz_lab)
{
    double p_parent = std::sqrt(px_parent*px_parent + py_parent*py_parent + pz_parent*pz_parent);
    if (p_parent < 1e-15) {
        E_lab = E_rest; px_lab = px_rest; py_lab = py_rest; pz_lab = pz_rest;
        return;
    }
    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;
    siren::math::Vector3D d(px_parent/p_parent, py_parent/p_parent, pz_parent/p_parent);

    double p_par = px_rest * d.GetX() + py_rest * d.GetY() + pz_rest * d.GetZ();
    double p_lab_par = gamma * (p_par + beta * E_rest);
    E_lab = gamma * (E_rest + beta * p_par);

    px_lab = px_rest + d.GetX() * (p_lab_par - p_par);
    py_lab = py_rest + d.GetY() * (p_lab_par - p_par);
    pz_lab = pz_rest + d.GetZ() * (p_lab_par - p_par);
}

// Boost 4-momentum from lab frame into the rest frame of 'parent'.
void BoostLabToRest(
    double E_parent, double px_parent, double py_parent, double pz_parent,
    double M_parent,
    double E_lab, double px_lab, double py_lab, double pz_lab,
    double & E_rest, double & px_rest, double & py_rest, double & pz_rest)
{
    // Inverse boost: flip the boost direction (negate beta)
    BoostRestToLab(E_parent, -px_parent, -py_parent, -pz_parent,
                   M_parent,
                   E_lab, px_lab, py_lab, pz_lab,
                   E_rest, px_rest, py_rest, pz_rest);
}

} // anonymous namespace

// ------------------------------------------------------------------ //
//  Direct mode constructor                                            //
// ------------------------------------------------------------------ //
DetectorDirected3BodyChannel::DetectorDirected3BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int directed_index,
    InvariantMassMode mass_mode,
    double resonance_mass,
    double resonance_width,
    double power_law_nu,
    double power_law_offset,
    DetectorDirected2BodyChannel::Mode mode,
    PhaseSpaceTopology topology,
    std::vector<double> mass_cdf_nodes,
    std::vector<double> mass_cdf_values)
    : factorization_(Factorization::Direct)
    , target_(std::move(target))
    , directed_index_(directed_index)
    , spectator_index_(-1)
    , pair_first_index_(-1)
    , pair_second_index_(-1)
    , directed_pair_index_(directed_index)
    , mass_mode_(mass_mode)
    , resonance_mass_(resonance_mass)
    , resonance_width_(resonance_width)
    , power_law_nu_(power_law_nu)
    , power_law_offset_(power_law_offset)
    , mass_cdf_nodes_(std::move(mass_cdf_nodes))
    , mass_cdf_values_(std::move(mass_cdf_values))
    , mode_(mode)
    , topology_(topology)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires a target geometry");
    }
    if (directed_index_ < 0 || directed_index_ > 2) {
        throw std::runtime_error("directed_index must be 0, 1, or 2");
    }
    // Infer the other two indices in ascending order
    int k = 0;
    int others[2];
    for (int i = 0; i < 3; ++i) {
        if (i != directed_index_) {
            others[k++] = i;
        }
    }
    other_a_index_ = others[0];
    other_b_index_ = others[1];

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

// ------------------------------------------------------------------ //
//  Recursive mode constructor (backward compatible)                    //
// ------------------------------------------------------------------ //
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
    PhaseSpaceTopology topology,
    std::vector<double> mass_cdf_nodes,
    std::vector<double> mass_cdf_values)
    : factorization_(Factorization::Recursive)
    , target_(std::move(target))
    , directed_index_(directed_pair_index)
    , other_a_index_(-1)
    , other_b_index_(-1)
    , spectator_index_(spectator_index)
    , pair_first_index_(pair_first_index)
    , pair_second_index_(pair_second_index)
    , directed_pair_index_(directed_pair_index)
    , mass_mode_(mass_mode)
    , resonance_mass_(resonance_mass)
    , resonance_width_(resonance_width)
    , power_law_nu_(power_law_nu)
    , power_law_offset_(power_law_offset)
    , mass_cdf_nodes_(std::move(mass_cdf_nodes))
    , mass_cdf_values_(std::move(mass_cdf_values))
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
    if (mass_mode_ == InvariantMassMode::Tabulated) {
        if (mass_cdf_nodes_.size() < 2 ||
            mass_cdf_nodes_.size() != mass_cdf_values_.size()) {
            throw std::runtime_error(
                "Tabulated invariant-mass mapping requires matching "
                "mass_cdf_nodes / mass_cdf_values arrays of length >= 2");
        }
        TabulatedMapping map(mass_cdf_nodes_, mass_cdf_values_, s_min, s_max);
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
    if (mass_mode_ == InvariantMassMode::Tabulated) {
        if (mass_cdf_nodes_.size() < 2 ||
            mass_cdf_nodes_.size() != mass_cdf_values_.size()) {
            return 0.0;
        }
        TabulatedMapping map(mass_cdf_nodes_, mass_cdf_values_, s_min, s_max);
        return map.Density(s);
    }
    UniformMapping map(s_min, s_max);
    return map.Density(s);
}

// ================================================================== //
//  Sample / Density dispatch                                          //
// ================================================================== //

void DetectorDirected3BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 3 ||
        record.secondary_masses.size() != 3 ||
        record.secondary_momenta.size() != 3) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires exactly 3 secondaries");
    }

    if (factorization_ == Factorization::Direct) {
        SampleDirect(random, detector_model, record);
    } else {
        SampleRecursive(random, detector_model, record);
    }
}

double DetectorDirected3BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 3 ||
        record.secondary_masses.size() != 3 ||
        record.secondary_momenta.size() != 3) {
        return 0.0;
    }

    if (factorization_ == Factorization::Direct) {
        return DensityDirect(detector_model, record);
    } else {
        return DensityRecursive(detector_model, record);
    }
}

// ================================================================== //
//  Direct mode: P -> directed + X, X -> other_a + other_b             //
// ================================================================== //

void DetectorDirected3BodyChannel::SampleDirect(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_dir = record.secondary_masses[directed_index_];
    double m_a = record.secondary_masses[other_a_index_];
    double m_b = record.secondary_masses[other_b_index_];

    // Invariant mass range for the complementary system X = (a + b)
    double s_X_min = (m_a + m_b) * (m_a + m_b);
    double s_X_max = (M - m_dir) * (M - m_dir);
    if (s_X_max <= s_X_min) {
        throw std::runtime_error("DetectorDirected3BodyChannel has no allowed invariant-mass range");
    }

    // Step 0: sample complementary invariant mass
    double s_X = SampleInvariantMassSquared(random, s_X_min, s_X_max);
    double M_X = std::sqrt(s_X);

    // Step 1: P -> directed + X as a 2-body decay.
    // Use SampleDirectedStep to bias the directed daughter toward
    // the target, using the PARENT boost.
    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    auto step1 = detail::SampleDirectedStep(
        init.E, init.px, init.py, init.pz,
        M, m_dir, M_X,
        vertex, *target_, target_volume_, mode_, random);

    double dir_E = step1.E_lab;
    double dir_px = step1.p_lab * step1.lab_dir.GetX();
    double dir_py = step1.p_lab * step1.lab_dir.GetY();
    double dir_pz = step1.p_lab * step1.lab_dir.GetZ();

    // X 4-momentum by conservation
    double X_E = init.E - dir_E;
    double X_px = init.px - dir_px;
    double X_py = init.py - dir_py;
    double X_pz = init.pz - dir_pz;

    // Step 2: X -> other_a + other_b (isotropic in X rest frame)
    double p_ab_rest = TwoBodyRestMomentum(M_X, m_a, m_b);
    double E_a_rest = TwoBodyRestEnergy(M_X, m_a, m_b);
    double E_b_rest = M_X - E_a_rest;

    double cos_ab = 2.0 * random->Uniform(0, 1) - 1.0;
    double phi_ab = 2.0 * M_PI * random->Uniform(0, 1);
    double sin_ab = std::sqrt(1.0 - cos_ab * cos_ab);

    double a_px_rest = p_ab_rest * sin_ab * std::cos(phi_ab);
    double a_py_rest = p_ab_rest * sin_ab * std::sin(phi_ab);
    double a_pz_rest = p_ab_rest * cos_ab;

    // Boost a and b from X rest frame to lab frame
    double a_E_lab, a_px_lab, a_py_lab, a_pz_lab;
    BoostRestToLab(X_E, X_px, X_py, X_pz, M_X,
                   E_a_rest, a_px_rest, a_py_rest, a_pz_rest,
                   a_E_lab, a_px_lab, a_py_lab, a_pz_lab);

    double b_E_lab = X_E - a_E_lab;
    double b_px_lab = X_px - a_px_lab;
    double b_py_lab = X_py - a_py_lab;
    double b_pz_lab = X_pz - a_pz_lab;

    // Write results
    record.secondary_momenta[directed_index_] = {dir_E, dir_px, dir_py, dir_pz};
    record.secondary_momenta[other_a_index_] = {a_E_lab, a_px_lab, a_py_lab, a_pz_lab};
    record.secondary_momenta[other_b_index_] = {b_E_lab, b_px_lab, b_py_lab, b_pz_lab};
    record.interaction_parameters["phase_space_s_X"] = s_X;
}

double DetectorDirected3BodyChannel::DensityDirect(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_dir = record.secondary_masses[directed_index_];
    double m_a = record.secondary_masses[other_a_index_];
    double m_b = record.secondary_masses[other_b_index_];

    double s_X_min = (m_a + m_b) * (m_a + m_b);
    double s_X_max = (M - m_dir) * (M - m_dir);
    if (s_X_max <= s_X_min) return 0.0;

    // Reconstruct complementary 4-vector X = other_a + other_b
    FourVector a_fv = ReadSecondary(record, other_a_index_);
    FourVector b_fv = ReadSecondary(record, other_b_index_);
    FourVector X_fv = Add(a_fv, b_fv);
    double s_X = MassSquared(X_fv);
    if (s_X <= 0.0) return 0.0;
    double M_X = std::sqrt(s_X);

    double s_density = InvariantMassDensity(s_X, s_X_min, s_X_max);
    if (s_density <= 0.0) return 0.0;

    // Step 1 density: P -> directed + X, directed biased toward target
    FourVector dir_fv = ReadSecondary(record, directed_index_);

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    double step1_density = detail::DensityDirectedStep(
        init.E, init.px, init.py, init.pz,
        M, m_dir, M_X,
        dir_fv.e, dir_fv.p.GetX(), dir_fv.p.GetY(), dir_fv.p.GetZ(),
        vertex, *target_, target_volume_, mode_);
    if (step1_density <= 0.0) return 0.0;

    // Step 2 density: X -> a + b isotropic
    static const double INV_FOUR_PI = 1.0 / (4.0 * M_PI);

    return s_density * step1_density * INV_FOUR_PI;
}

// ================================================================== //
//  Recursive mode: P -> spectator + pair, pair -> directed + other    //
//  (original implementation, preserved for backward compatibility)    //
// ================================================================== //

void DetectorDirected3BodyChannel::SampleRecursive(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
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

    double s_pair = SampleInvariantMassSquared(random, s_min, s_max);
    double m_pair = std::sqrt(s_pair);

    double E_parent = init.E;
    double px_parent = init.px;
    double py_parent = init.py;
    double pz_parent = init.pz;
    double p_parent = std::sqrt(px_parent*px_parent + py_parent*py_parent + pz_parent*pz_parent);

    // Step 1: P -> pair + spectator (isotropic in rest frame)
    double p1_rest = TwoBodyRestMomentum(M, m_pair, m_spectator);
    double E1_rest = TwoBodyRestEnergy(M, m_pair, m_spectator);

    double cos1 = 2.0 * random->Uniform(0, 1) - 1.0;
    double phi1 = 2.0 * M_PI * random->Uniform(0, 1);
    double sin1 = std::sqrt(1.0 - cos1 * cos1);

    double pair_E, pair_px, pair_py, pair_pz;
    if (p_parent < 1e-15) {
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

    record.secondary_momenta[spectator_index_] = {spec_E, spec_px, spec_py, spec_pz};
    record.secondary_momenta[directed_pair_index_] = {dir_E, dir_px, dir_py, dir_pz};
    record.secondary_momenta[other_pair_index] = {oth_E, oth_px, oth_py, oth_pz};
    record.interaction_parameters["phase_space_s_pair"] = s_pair;
}

double DetectorDirected3BodyChannel::DensityRecursive(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
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

    double s_density = InvariantMassDensity(s_pair, s_min, s_max);
    if (s_density <= 0.0) return 0.0;

    static const double INV_FOUR_PI = 1.0 / (4.0 * M_PI);
    double step1_density = INV_FOUR_PI;

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
