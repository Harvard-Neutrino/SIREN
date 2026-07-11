#include "SIREN/injection/DetectorDirected2BodyChannel.h"

#include "DetectorDirectedChannelUtils.h"
#include "InteractionRecordUtils.h"

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
    target_volume_ = ResolveDetectorDirectedVolume(
        *target_, mode_ == Mode::Volume, volume);
}

void DetectorDirected2BodyChannel::SetVolume(double volume) {
    target_volume_ = volume;
}

bool DetectorDirected2BodyChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2)) {
        return false;
    }
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);
    auto parent = detail::ReadPrimary(record);
    auto geo = detail::ClassifyDirectedRegime(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
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
    detail::RequireSecondaryStorage(
        record, 2, "DetectorDirected2BodyChannel");

    auto parent = detail::ReadPrimary(record);
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);

    auto result = detail::SampleDirectedStep(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_, target_volume_, mode_, random);

    detail::FourVector daughter{
        result.E_lab,
        result.lab_dir * result.p_lab};
    detail::WriteSecondary(record, daughter_index_, daughter);
    detail::WriteSecondary(record, 1 - daughter_index_, {
        parent.e - daughter.e, parent.p - daughter.p});
}

double DetectorDirected2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2)) {
        return 0.0;
    }

    auto parent = detail::ReadPrimary(record);
    auto daughter = detail::ReadSecondary(record, daughter_index_);
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);

    return detail::DensityDirectedStep(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        daughter.e,
        daughter.p.GetX(), daughter.p.GetY(), daughter.p.GetZ(),
        decay_pos, *target_, target_volume_, mode_);
}

} // namespace injection
} // namespace siren
