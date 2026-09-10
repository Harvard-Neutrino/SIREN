#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"

#include "DetectorDirectedChannelUtils.h"
#include "InteractionRecordUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace siren {
namespace injection {

DetectorDirectedAngularSectorChannel::DetectorDirectedAngularSectorChannel(
    std::shared_ptr<siren::geometry::Geometry const> target,
    double u_lo,
    double u_hi,
    double phi_lo,
    double phi_hi,
    int daughter_index)
    : target_(std::move(target))
    , u_lo_(u_lo)
    , u_hi_(u_hi)
    , phi_lo_(phi_lo)
    , phi_hi_(phi_hi)
    , daughter_index_(daughter_index)
{
    if (!target_) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires a target geometry");
    }
    if (daughter_index_ != 0 && daughter_index_ != 1) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel daughter_index must be 0 or 1");
    }
    if (!(u_hi_ > u_lo_) || u_lo_ < 0.0 || u_hi_ > 1.0) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires 0 <= u_lo < u_hi <= 1");
    }
    double phi_width = phi_hi_ - phi_lo_;
    if (!std::isfinite(phi_lo_) || !std::isfinite(phi_hi_) ||
        !(phi_width > 0.0) ||
        phi_width > detail::kTwoPi + 1e-12) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires a finite phi range "
            "with 0 < phi_hi - phi_lo <= 2*pi");
    }

    // Store a canonical unwrapped interval [start, start + width]. The upper
    // endpoint may exceed 2*pi for a sector that crosses the branch cut;
    // Sample uses that interval directly and Density tests membership modulo
    // 2*pi.
    if (phi_width > detail::kTwoPi) phi_width = detail::kTwoPi;
    phi_lo_ = std::fmod(phi_lo_, detail::kTwoPi);
    if (phi_lo_ < 0.0) phi_lo_ += detail::kTwoPi;
    phi_hi_ = phi_lo_ + phi_width;
}

void DetectorDirectedAngularSectorChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    detail::RequireSecondaryStorage(
        record, 2, "DetectorDirectedAngularSectorChannel");

    auto parent = detail::ReadPrimary(record);
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);

    detail::AngularSectorBin bin{u_lo_, u_hi_, phi_lo_, phi_hi_};

    auto result = detail::SampleAngularSectorStep(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_, bin, random);

    detail::FourVector daughter{
        result.E_lab,
        result.lab_dir * result.p_lab};
    detail::WriteSecondary(record, daughter_index_, daughter);
    detail::WriteSecondary(record, 1 - daughter_index_, {
        parent.e - daughter.e, parent.p - daughter.p});
}

double DetectorDirectedAngularSectorChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2)) {
        return 0.0;
    }

    auto parent = detail::ReadPrimary(record);
    auto daughter = detail::ReadSecondary(record, daughter_index_);
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);

    detail::AngularSectorBin bin{u_lo_, u_hi_, phi_lo_, phi_hi_};

    return detail::DensityAngularSectorStep(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        daughter.e,
        daughter.p.GetX(), daughter.p.GetY(), daughter.p.GetZ(),
        decay_pos, *target_, bin);
}

bool DetectorDirectedAngularSectorChannel::DirectingActive(
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2)) {
        return false;
    }
    auto parent = detail::ReadPrimary(record);
    siren::math::Vector3D decay_pos = detail::ReadVertex(record);
    detail::DirectedGeometry geo = detail::ClassifyDirectedRegime(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_);
    return detail::SectorActive(geo.regime, geo.inside_geometry);
}

} // namespace injection
} // namespace siren
