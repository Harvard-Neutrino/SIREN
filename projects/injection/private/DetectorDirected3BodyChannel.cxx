#include "SIREN/injection/DetectorDirected3BodyChannel.h"

#include "DetectorDirectedChannelUtils.h"
#include "InteractionRecordUtils.h"
#include "LorentzBoostUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/injection/InvariantMassMapping.h"
#include "SIREN/injection/TwoBodyKinematics.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

namespace {

using detail::FourVector;
using detail::ReadSecondary;

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
    auto primary = detail::ReadPrimary(record);
    double E = primary.e;
    double px = primary.p.GetX();
    double py = primary.p.GetY();
    double pz = primary.p.GetZ();

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

std::shared_ptr<TabulatedMappingTable const> BuildInvariantMassTable(
    DetectorDirected3BodyChannel::InvariantMassMode mode,
    std::vector<double> nodes,
    std::vector<double> values)
{
    using Mode = DetectorDirected3BodyChannel::InvariantMassMode;
    switch (mode) {
    case Mode::Uniform:
    case Mode::BreitWigner:
    case Mode::PowerLaw:
        return nullptr;
    case Mode::Tabulated:
        return std::make_shared<TabulatedMappingTable>(
            std::move(nodes), std::move(values));
    }
    throw std::invalid_argument(
        "DetectorDirected3BodyChannel received an unknown invariant-mass mode");
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
    std::vector<double> mass_cdf_values,
    double volume)
    : target_(std::move(target))
    , roles_{directed_index, -1, -1, directed_index}
    , sample_implementation_(&DetectorDirected3BodyChannel::SampleDirect)
    , density_implementation_(&DetectorDirected3BodyChannel::DensityDirect)
    , active_implementation_(&DetectorDirected3BodyChannel::DirectingActiveDirect)
    , mass_mode_(mass_mode)
    , resonance_mass_(resonance_mass)
    , resonance_width_(resonance_width)
    , power_law_nu_(power_law_nu)
    , power_law_offset_(power_law_offset)
    , mass_cdf_table_(BuildInvariantMassTable(
          mass_mode, std::move(mass_cdf_nodes), std::move(mass_cdf_values)))
    , mode_(mode)
    , topology_(topology)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires a target geometry");
    }
    if (roles_.directed < 0 || roles_.directed > 2) {
        throw std::runtime_error("directed_index must be 0, 1, or 2");
    }
    // Infer the other two indices in ascending order
    int k = 0;
    int others[2];
    for (int i = 0; i < 3; ++i) {
        if (i != roles_.directed) {
            others[k++] = i;
        }
    }
    roles_.pair_first = others[0];
    roles_.pair_second = others[1];

    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == DetectorDirected2BodyChannel::Mode::Volume, volume);
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
    std::vector<double> mass_cdf_values,
    double volume)
    : target_(std::move(target))
    , roles_{spectator_index, pair_first_index, pair_second_index,
             directed_pair_index}
    , sample_implementation_(&DetectorDirected3BodyChannel::SampleRecursive)
    , density_implementation_(&DetectorDirected3BodyChannel::DensityRecursive)
    , active_implementation_(&DetectorDirected3BodyChannel::DirectingActiveRecursive)
    , mass_mode_(mass_mode)
    , resonance_mass_(resonance_mass)
    , resonance_width_(resonance_width)
    , power_law_nu_(power_law_nu)
    , power_law_offset_(power_law_offset)
    , mass_cdf_table_(BuildInvariantMassTable(
          mass_mode, std::move(mass_cdf_nodes), std::move(mass_cdf_values)))
    , mode_(mode)
    , topology_(topology)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected3BodyChannel requires a target geometry");
    }
    std::array<int, 3> roles{
        roles_.outer, roles_.pair_first, roles_.pair_second};
    std::sort(roles.begin(), roles.end());
    if (roles != std::array<int, 3>{0, 1, 2}) {
        throw std::runtime_error(
            "DetectorDirected3BodyChannel recursive indices must be a "
            "permutation of 0, 1, and 2");
    }
    if (roles_.directed != roles_.pair_first &&
        roles_.directed != roles_.pair_second) {
        throw std::runtime_error("directed_pair_index must be one of the pair indices");
    }
    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == DetectorDirected2BodyChannel::Mode::Volume, volume);
}

void DetectorDirected3BodyChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

bool DetectorDirected3BodyChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 3)) return false;
    return (this->*active_implementation_)(record);
}

bool DetectorDirected3BodyChannel::DirectingActiveDirect(
    siren::dataclasses::InteractionRecord const & record) const
{
    InitialState init = GetInitialState(record, topology_);
    siren::math::Vector3D vertex = detail::ReadVertex(record);
    FourVector first = ReadSecondary(record, roles_.pair_first);
    FourVector second = ReadSecondary(record, roles_.pair_second);
    double pair_mass_squared = MassSquared(Add(first, second));
    if (pair_mass_squared <= 0.0) return false;
    auto geometry = detail::ClassifyDirectedRegime(
        init.E, init.px, init.py, init.pz, init.M,
        record.secondary_masses[roles_.directed],
        std::sqrt(pair_mass_squared), vertex, *target_);
    return detail::IsDirectedStepActive(
        geometry.regime, geometry.inside_geometry);
}

bool DetectorDirected3BodyChannel::DirectingActiveRecursive(
    siren::dataclasses::InteractionRecord const & record) const
{
    siren::math::Vector3D vertex = detail::ReadVertex(record);
    FourVector first = ReadSecondary(record, roles_.pair_first);
    FourVector second = ReadSecondary(record, roles_.pair_second);
    FourVector pair = Add(first, second);
    double pair_mass_squared = MassSquared(pair);
    if (pair_mass_squared <= 0.0) return false;
    int other_pair_index = roles_.PairOther();
    auto geometry = detail::ClassifyDirectedRegime(
        pair.e, pair.p.GetX(), pair.p.GetY(), pair.p.GetZ(),
        std::sqrt(pair_mass_squared),
        record.secondary_masses[roles_.directed],
        record.secondary_masses[other_pair_index], vertex, *target_);
    return detail::IsDirectedStepActive(
        geometry.regime, geometry.inside_geometry);
}

double DetectorDirected3BodyChannel::SampleInvariantMassSquared(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    double s_min,
    double s_max) const
{
    double r = random->Uniform(0, 1);
    switch (mass_mode_) {
    case InvariantMassMode::Uniform: {
        UniformMapping map(s_min, s_max);
        return map.Forward(r);
    }
    case InvariantMassMode::BreitWigner: {
        if (resonance_mass_ <= 0.0 || resonance_width_ <= 0.0) {
            throw std::runtime_error("BreitWigner invariant-mass mapping requires positive mass and width");
        }
        BreitWignerMapping map(resonance_mass_, resonance_width_, s_min, s_max);
        return map.Forward(r);
    }
    case InvariantMassMode::PowerLaw: {
        if (s_min < power_law_offset_) {
            throw std::runtime_error(
                "PowerLaw invariant-mass mapping requires s_min >= power_law_offset");
        }
        PowerLawMapping map(power_law_nu_, power_law_offset_, s_min, s_max);
        return map.Forward(r);
    }
    case InvariantMassMode::Tabulated: {
        if (!mass_cdf_table_) {
            throw std::runtime_error(
                "Tabulated invariant-mass mapping requires matching "
                "mass_cdf_nodes / mass_cdf_values arrays of length >= 2");
        }
        TabulatedMapping map(mass_cdf_table_, s_min, s_max);
        if (!map.HasSupport()) {
            throw siren::utilities::InjectionFailure(
                "Tabulated invariant-mass mapping has no support in the event's "
                "kinematic range");
        }
        return map.Forward(r);
    }
    }
    throw std::logic_error("Unknown invariant-mass mode");
}

double DetectorDirected3BodyChannel::InvariantMassDensity(
    double s,
    double s_min,
    double s_max) const
{
    if (s < s_min || s > s_max || s_max <= s_min) return 0.0;
    switch (mass_mode_) {
    case InvariantMassMode::Uniform: {
        UniformMapping map(s_min, s_max);
        return map.Density(s);
    }
    case InvariantMassMode::BreitWigner: {
        if (resonance_mass_ <= 0.0 || resonance_width_ <= 0.0) return 0.0;
        BreitWignerMapping map(resonance_mass_, resonance_width_, s_min, s_max);
        return map.Density(s);
    }
    case InvariantMassMode::PowerLaw: {
        if (s_min < power_law_offset_) return 0.0;
        PowerLawMapping map(power_law_nu_, power_law_offset_, s_min, s_max);
        return map.Density(s);
    }
    case InvariantMassMode::Tabulated: {
        if (!mass_cdf_table_) return 0.0;
        TabulatedMapping map(mass_cdf_table_, s_min, s_max);
        return map.Density(s);
    }
    }
    throw std::logic_error("Unknown invariant-mass mode");
}

// ================================================================== //
//  Sample / Density dispatch                                          //
// ================================================================== //

void DetectorDirected3BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord & record) const
{
    detail::RequireSecondaryStorage(
        record, 3, "DetectorDirected3BodyChannel");

    (this->*sample_implementation_)(random, detector_model, record);
}

double DetectorDirected3BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const> detector_model,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 3)) {
        return 0.0;
    }

    return (this->*density_implementation_)(detector_model, record);
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
    double m_dir = record.secondary_masses[roles_.directed];
    double m_a = record.secondary_masses[roles_.pair_first];
    double m_b = record.secondary_masses[roles_.pair_second];

    // Invariant mass range for the complementary system X = (a + b)
    double s_X_min = (m_a + m_b) * (m_a + m_b);
    double s_X_max = (M - m_dir) * (M - m_dir);
    if (s_X_max <= s_X_min) {
        throw siren::utilities::InjectionFailure(
            "DetectorDirected3BodyChannel has no allowed invariant-mass range");
    }

    // Step 0: sample complementary invariant mass
    double s_X = SampleInvariantMassSquared(random, s_X_min, s_X_max);
    double M_X = std::sqrt(s_X);

    // Step 1: P -> directed + X as a 2-body decay.
    // Use SampleDirectedStep to bias the directed daughter toward
    // the target, using the PARENT boost.
    siren::math::Vector3D vertex = detail::ReadVertex(record);

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

    double cos_ab = 2.0 * random->Uniform(0, 1) - 1.0;
    double phi_ab = 2.0 * M_PI * random->Uniform(0, 1);
    double sin_ab = std::sqrt(1.0 - cos_ab * cos_ab);

    double a_px_rest = p_ab_rest * sin_ab * std::cos(phi_ab);
    double a_py_rest = p_ab_rest * sin_ab * std::sin(phi_ab);
    double a_pz_rest = p_ab_rest * cos_ab;

    // Boost a and b from X rest frame to lab frame
    auto a_lab = detail::BoostRestFrameToLab(
        X_E, X_px, X_py, X_pz,
        E_a_rest, a_px_rest, a_py_rest, a_pz_rest);
    double a_E_lab = a_lab[0];
    double a_px_lab = a_lab[1];
    double a_py_lab = a_lab[2];
    double a_pz_lab = a_lab[3];

    double b_E_lab = X_E - a_E_lab;
    double b_px_lab = X_px - a_px_lab;
    double b_py_lab = X_py - a_py_lab;
    double b_pz_lab = X_pz - a_pz_lab;

    // Write results
    detail::WriteSecondary(record, roles_.directed, {
        dir_E, siren::math::Vector3D(dir_px, dir_py, dir_pz)});
    detail::WriteSecondary(record, roles_.pair_first, {
        a_E_lab, siren::math::Vector3D(a_px_lab, a_py_lab, a_pz_lab)});
    detail::WriteSecondary(record, roles_.pair_second, {
        b_E_lab, siren::math::Vector3D(b_px_lab, b_py_lab, b_pz_lab)});
    record.interaction_parameters["phase_space_s_X"] = s_X;
}

double DetectorDirected3BodyChannel::DensityDirect(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_dir = record.secondary_masses[roles_.directed];
    double m_a = record.secondary_masses[roles_.pair_first];
    double m_b = record.secondary_masses[roles_.pair_second];

    double s_X_min = (m_a + m_b) * (m_a + m_b);
    double s_X_max = (M - m_dir) * (M - m_dir);
    if (s_X_max <= s_X_min) return 0.0;

    // Reconstruct complementary 4-vector X = other_a + other_b
    FourVector a_fv = ReadSecondary(record, roles_.pair_first);
    FourVector b_fv = ReadSecondary(record, roles_.pair_second);
    FourVector X_fv = Add(a_fv, b_fv);
    double s_X = MassSquared(X_fv);
    if (s_X <= 0.0) return 0.0;
    double M_X = std::sqrt(s_X);

    double s_density = InvariantMassDensity(s_X, s_X_min, s_X_max);
    if (s_density <= 0.0) return 0.0;

    // Step 1 density: P -> directed + X, directed biased toward target
    FourVector dir_fv = ReadSecondary(record, roles_.directed);

    siren::math::Vector3D vertex = detail::ReadVertex(record);

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
// ================================================================== //

void DetectorDirected3BodyChannel::SampleRecursive(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_spectator = record.secondary_masses[roles_.outer];
    double m_first = record.secondary_masses[roles_.pair_first];
    double m_second = record.secondary_masses[roles_.pair_second];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) {
        throw siren::utilities::InjectionFailure(
            "DetectorDirected3BodyChannel has no allowed invariant-mass range");
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

    siren::math::Vector3D pair_rest(
        p1_rest * sin1 * std::cos(phi1),
        p1_rest * sin1 * std::sin(phi1),
        p1_rest * cos1);
    if (p_parent >= 1e-15) {
        siren::math::Vector3D parent_dir(px_parent/p_parent, py_parent/p_parent, pz_parent/p_parent);
        siren::math::Vector3D perp1, perp2;
        detail::SectorPerpFrame(parent_dir, perp1, perp2);
        pair_rest = parent_dir * (p1_rest * cos1)
                  + perp1 * (p1_rest * sin1 * std::cos(phi1))
                  + perp2 * (p1_rest * sin1 * std::sin(phi1));
    }
    auto pair_lab = detail::BoostRestFrameToLab(
        E_parent, px_parent, py_parent, pz_parent,
        E1_rest, pair_rest.GetX(), pair_rest.GetY(), pair_rest.GetZ());
    double pair_E = pair_lab[0];
    double pair_px = pair_lab[1];
    double pair_py = pair_lab[2];
    double pair_pz = pair_lab[3];

    double spec_E = E_parent - pair_E;
    double spec_px = px_parent - pair_px;
    double spec_py = py_parent - pair_py;
    double spec_pz = pz_parent - pair_pz;

    // Step 2: pair -> directed + other (directed toward target)
    int other_pair_index = roles_.PairOther();
    double m_directed = record.secondary_masses[roles_.directed];
    double m_other = record.secondary_masses[other_pair_index];

    siren::math::Vector3D vertex = detail::ReadVertex(record);

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

    detail::WriteSecondary(record, roles_.outer, {
        spec_E, siren::math::Vector3D(spec_px, spec_py, spec_pz)});
    detail::WriteSecondary(record, roles_.directed, {
        dir_E, siren::math::Vector3D(dir_px, dir_py, dir_pz)});
    detail::WriteSecondary(record, other_pair_index, {
        oth_E, siren::math::Vector3D(oth_px, oth_py, oth_pz)});
    record.interaction_parameters["phase_space_s_pair"] = s_pair;
}

double DetectorDirected3BodyChannel::DensityRecursive(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    InitialState init = GetInitialState(record, topology_);
    double M = init.M;
    double m_spectator = record.secondary_masses[roles_.outer];
    double m_first = record.secondary_masses[roles_.pair_first];
    double m_second = record.secondary_masses[roles_.pair_second];
    double s_min = (m_first + m_second) * (m_first + m_second);
    double s_max = (M - m_spectator) * (M - m_spectator);
    if (s_max <= s_min) return 0.0;

    FourVector first = ReadSecondary(record, roles_.pair_first);
    FourVector second = ReadSecondary(record, roles_.pair_second);
    FourVector pair = Add(first, second);
    double s_pair = MassSquared(pair);
    if (s_pair <= 0.0) return 0.0;
    double m_pair = std::sqrt(s_pair);

    double s_density = InvariantMassDensity(s_pair, s_min, s_max);
    if (s_density <= 0.0) return 0.0;

    static const double INV_FOUR_PI = 1.0 / (4.0 * M_PI);
    double step1_density = INV_FOUR_PI;

    int other_pair_index = roles_.PairOther();
    FourVector directed = ReadSecondary(record, roles_.directed);

    siren::math::Vector3D vertex = detail::ReadVertex(record);

    double step2_density = detail::DensityDirectedStep(
        pair.e, pair.p.GetX(), pair.p.GetY(), pair.p.GetZ(),
        m_pair,
        record.secondary_masses[roles_.directed],
        record.secondary_masses[other_pair_index],
        directed.e, directed.p.GetX(), directed.p.GetY(), directed.p.GetZ(),
        vertex, *target_, target_volume_, mode_);
    if (step2_density <= 0.0) return 0.0;

    return s_density * step1_density * step2_density;
}

} // namespace injection
} // namespace siren
