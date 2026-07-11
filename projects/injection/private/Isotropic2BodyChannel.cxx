#include "SIREN/injection/Isotropic2BodyChannel.h"

#include "LorentzBoostUtils.h"
#include "InteractionRecordUtils.h"
#include "SIREN/injection/TwoBodyKinematics.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Errors.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/math/Vector3D.h"

#include <cmath>
#include <stdexcept>

namespace siren {
namespace injection {

static const double TWO_PI = 2.0 * M_PI;
static const double FOUR_PI = 4.0 * M_PI;

namespace {

bool HasTwoBodyDecayPhaseSpace(double parent_mass,
                               double daughter_mass,
                               double other_mass)
{
    return std::isfinite(parent_mass) &&
           std::isfinite(daughter_mass) &&
           std::isfinite(other_mass) &&
           parent_mass > 0.0 &&
           daughter_mass >= 0.0 &&
           other_mass >= 0.0 &&
           parent_mass >= daughter_mass + other_mass;
}

} // namespace

Isotropic2BodyChannel::Isotropic2BodyChannel(int daughter_index)
    : daughter_index_(daughter_index)
{
    if (daughter_index_ != 0 && daughter_index_ != 1) {
        throw std::runtime_error(
            "Isotropic2BodyChannel daughter_index must be 0 or 1");
    }
}

void Isotropic2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    detail::RequireSecondaryStorage(record, 2, "Isotropic2BodyChannel");

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    if (!HasTwoBodyDecayPhaseSpace(M_parent, m_A, m_B)) {
        throw siren::utilities::InjectionFailure(
            "Isotropic2BodyChannel has no kinematically allowed two-body "
            "phase space");
    }

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);

    // Sample isotropic direction in rest frame
    double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = TWO_PI * random->Uniform(0, 1);

    // Rest-frame 3-momentum of daughter A
    double px_rest = p_rest * sin_theta * std::cos(phi);
    double py_rest = p_rest * sin_theta * std::sin(phi);
    double pz_rest = p_rest * cos_theta;

    // Boost to lab frame along parent direction
    auto parent = detail::ReadPrimary(record);
    auto lab = detail::BoostRestFrameToLab(
        parent.e, parent.p.GetX(), parent.p.GetY(), parent.p.GetZ(),
        E_A_rest, px_rest, py_rest, pz_rest);
    detail::FourVector daughter{
        lab[0], siren::math::Vector3D(lab[1], lab[2], lab[3])};
    detail::WriteSecondary(record, daughter_index_, daughter);

    // Daughter B: by momentum conservation
    detail::WriteSecondary(record, 1 - daughter_index_, {
        parent.e - daughter.e, parent.p - daughter.p});
}

double Isotropic2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    if (!detail::HasSecondaryStorage(record, 2)) {
        return 0.0;
    }
    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];
    if (!HasTwoBodyDecayPhaseSpace(M_parent, m_A, m_B)) {
        return 0.0;
    }
    // Isotropic in rest-frame solid angle: density = 1/(4*pi)
    return 1.0 / FOUR_PI;
}

} // namespace injection
} // namespace siren
