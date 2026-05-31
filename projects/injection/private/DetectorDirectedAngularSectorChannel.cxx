#include "SIREN/injection/DetectorDirectedAngularSectorChannel.h"

#include "DetectorDirectedChannelUtils.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/geometry/Geometry.h"
#include "SIREN/math/Vector3D.h"
#include "SIREN/utilities/Random.h"

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
    if (!(u_hi_ > u_lo_) || u_lo_ < 0.0 || u_hi_ > 1.0) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires 0 <= u_lo < u_hi <= 1");
    }
    if (!(phi_hi_ > phi_lo_)) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires phi_lo < phi_hi");
    }
}

void DetectorDirectedAngularSectorChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2) {
        throw std::runtime_error(
            "DetectorDirectedAngularSectorChannel requires exactly 2 secondaries");
    }

    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];

    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    detail::AngularSectorBin bin{u_lo_, u_hi_, phi_lo_, phi_hi_};

    auto result = detail::SampleAngularSectorStep(
        E_parent, px_parent, py_parent, pz_parent,
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        decay_pos, *target_, bin, random);

    record.secondary_momenta[daughter_index_] = {
        result.E_lab,
        result.p_lab * result.lab_dir.GetX(),
        result.p_lab * result.lab_dir.GetY(),
        result.p_lab * result.lab_dir.GetZ()
    };
    record.secondary_momenta[1 - daughter_index_] = {
        E_parent - result.E_lab,
        px_parent - result.p_lab * result.lab_dir.GetX(),
        py_parent - result.p_lab * result.lab_dir.GetY(),
        pz_parent - result.p_lab * result.lab_dir.GetZ()
    };
}

double DetectorDirectedAngularSectorChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (record.signature.secondary_types.size() != 2) return 0.0;

    siren::math::Vector3D decay_pos(
        record.interaction_vertex[0],
        record.interaction_vertex[1],
        record.interaction_vertex[2]);

    detail::AngularSectorBin bin{u_lo_, u_hi_, phi_lo_, phi_hi_};

    return detail::DensityAngularSectorStep(
        record.primary_momentum[0],
        record.primary_momentum[1],
        record.primary_momentum[2],
        record.primary_momentum[3],
        record.primary_mass,
        record.secondary_masses[daughter_index_],
        record.secondary_masses[1 - daughter_index_],
        record.secondary_momenta[daughter_index_][0],
        record.secondary_momenta[daughter_index_][1],
        record.secondary_momenta[daughter_index_][2],
        record.secondary_momenta[daughter_index_][3],
        decay_pos, *target_, bin);
}

} // namespace injection
} // namespace siren
