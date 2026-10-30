#include "SIREN/injection/DetectorDirected2BodyChannel.h"

#include "DetectorDirectedChannelUtils.h"
#include "SIREN/dataclasses/RecordKinematics.h"
#include "OnShellDecayKinematics.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/injection/GeometryVolume.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

DetectorDirected2BodyChannel::DetectorDirected2BodyChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    int daughter_index,
    Mode mode,
    double volume)
    : target_(std::move(target))
    , daughter_index_(daughter_index)
    , mode_(mode)
{
    if (!target_) {
        throw std::runtime_error("DetectorDirected2BodyChannel requires a target geometry");
    }
    if (daughter_index_ != 0 && daughter_index_ != 1) {
        throw std::runtime_error(
            "DetectorDirected2BodyChannel daughter_index must be 0 or 1");
    }
    if (mode_ != Mode::Cone && mode_ != Mode::Volume) {
        throw std::invalid_argument("Unknown directed two-body mode");
    }
    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == Mode::Volume, volume);
}


void DetectorDirected2BodyChannel::SetVolume(double volume) {
    // Validated like the constructor argument.
    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == Mode::Volume, volume);
}

bool DetectorDirected2BodyChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!siren::dataclasses::HasSecondaryStorage(record, 2)
        || !detail::OnShellParentValid(record)) {
        return false;
    }
    siren::math::Vector3D decay_pos = siren::dataclasses::InteractionVertex(record);
    auto parent = detail::OnShellParent(record);
    auto geo = detail::ClassifyDirectedRegime(
        parent.Energy(), parent.Momentum().GetX(), parent.Momentum().GetY(), parent.Momentum().GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_);
    return detail::IsDirectedStepActive(geo.regime, geo.inside_geometry);
}

// Regime classification, sampling, and density evaluation are in
// DetectorDirectedChannelUtils.h (detail::ClassifyDirectedRegime,
// detail::SampleDirectedStep, detail::DensityDirectedStep).

// ================================================================ //
//  Sample and Density -- thin adapters over shared step functions  //
// ================================================================ //

void DetectorDirected2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    siren::dataclasses::RequireSecondaryStorage(
        record, 2, "DetectorDirected2BodyChannel");
    detail::RequireOnShellParent(record);

    // Decay frame from the declared mass and three-momentum.
    auto parent = detail::OnShellParent(record);
    siren::math::Vector3D decay_pos = siren::dataclasses::InteractionVertex(record);

    auto result = detail::SampleDirectedStep(
        parent.Energy(), parent.Momentum().GetX(), parent.Momentum().GetY(), parent.Momentum().GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_, target_volume_, mode_, random);

    auto daughter = siren::math::FourMomentum::FromMass(
        record.secondary_masses[daughter_index_], result.lab_dir * result.p_lab);
    siren::dataclasses::SetSecondaryMomentum(record, daughter_index_, daughter);
    siren::dataclasses::SetSecondaryMomentum(record, 1 - daughter_index_, parent - daughter);
}

double DetectorDirected2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!siren::dataclasses::HasSecondaryStorage(record, 2)
        || !detail::OnShellParentValid(record)) {
        return 0.0;
    }

    auto parent = detail::OnShellParent(record);
    auto daughter = siren::dataclasses::SecondaryMomentum(record, daughter_index_);
    siren::math::Vector3D decay_pos = siren::dataclasses::InteractionVertex(record);

    return detail::DensityDirectedStep(
        parent.Energy(), parent.Momentum().GetX(), parent.Momentum().GetY(), parent.Momentum().GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        daughter.Energy(),
        daughter.Momentum().GetX(), daughter.Momentum().GetY(), daughter.Momentum().GetZ(),
        decay_pos, *target_, target_volume_, mode_);
}

} // namespace injection
} // namespace siren
