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

    double M = record.primary_mass;
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

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Step 1: P -> pair + spectator (directed toward target)
    auto step1 = detail::SampleDirectedStep(
        record.primary_momentum[0],
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3],
        M, m_pair, m_spectator,
        vertex, *target_, target_volume_, mode_, random);

    // Construct pair and spectator 4-vectors
    double pair_E = step1.E_lab;
    double pair_px = step1.p_lab * step1.lab_dir.GetX();
    double pair_py = step1.p_lab * step1.lab_dir.GetY();
    double pair_pz = step1.p_lab * step1.lab_dir.GetZ();

    double spec_E = record.primary_momentum[0] - pair_E;
    double spec_px = record.primary_momentum[1] - pair_px;
    double spec_py = record.primary_momentum[2] - pair_py;
    double spec_pz = record.primary_momentum[3] - pair_pz;

    // Step 2: pair -> directed + other (directed toward target)
    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_ : pair_first_index_;
    double m_directed = record.secondary_masses[directed_pair_index_];
    double m_other = record.secondary_masses[other_pair_index];

    auto step2 = detail::SampleDirectedStep(
        pair_E, pair_px, pair_py, pair_pz,
        m_pair, m_directed, m_other,
        vertex, *target_, target_volume_, mode_, random);

    // Construct directed and other daughter 4-vectors
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

    double M = record.primary_mass;
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

    siren::math::Vector3D vertex(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    // Step 1 density: P -> pair + spectator
    double step1_density = detail::DensityDirectedStep(
        record.primary_momentum[0],
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3],
        M, m_pair, m_spectator,
        pair.e, pair.p.GetX(), pair.p.GetY(), pair.p.GetZ(),
        vertex, *target_, target_volume_, mode_);
    if (step1_density <= 0.0) return 0.0;

    // Step 2 density: pair -> directed + other
    int other_pair_index = (directed_pair_index_ == pair_first_index_)
        ? pair_second_index_ : pair_first_index_;
    FourVector directed = ReadSecondary(record, directed_pair_index_);

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
