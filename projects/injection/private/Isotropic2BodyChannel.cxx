#include "SIREN/injection/Isotropic2BodyChannel.h"
#include "SIREN/injection/TwoBodyKinematics.h"

#include "SIREN/dataclasses/InteractionRecord.h"
#include "SIREN/utilities/Random.h"
#include "SIREN/math/Vector3D.h"

#include <cmath>
#include <stdexcept>

namespace siren {
namespace injection {

static const double TWO_PI = 2.0 * M_PI;
static const double FOUR_PI = 4.0 * M_PI;

Isotropic2BodyChannel::Isotropic2BodyChannel(int daughter_index)
    : daughter_index_(daughter_index) {}

void Isotropic2BodyChannel::Sample(
    std::shared_ptr<siren::utilities::SIREN_random> random,
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord & record) const
{
    if (record.signature.secondary_types.size() != 2) {
        throw std::runtime_error(
            "Isotropic2BodyChannel requires exactly 2 secondary particles");
    }

    double M_parent = record.primary_mass;
    double m_A = record.secondary_masses[daughter_index_];
    double m_B = record.secondary_masses[1 - daughter_index_];

    double p_rest = TwoBodyRestMomentum(M_parent, m_A, m_B);
    double E_A_rest = TwoBodyRestEnergy(M_parent, m_A, m_B);
    double E_B_rest = TwoBodyRestEnergy(M_parent, m_B, m_A);

    // Sample isotropic direction in rest frame
    double cos_theta = 2.0 * random->Uniform(0, 1) - 1.0;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = TWO_PI * random->Uniform(0, 1);

    // Rest-frame 3-momentum of daughter A
    double px_rest = p_rest * sin_theta * std::cos(phi);
    double py_rest = p_rest * sin_theta * std::sin(phi);
    double pz_rest = p_rest * cos_theta;

    // Boost to lab frame along parent direction
    double E_parent = record.primary_momentum[0];
    double px_parent = record.primary_momentum[1];
    double py_parent = record.primary_momentum[2];
    double pz_parent = record.primary_momentum[3];
    double p_parent = std::sqrt(px_parent * px_parent + py_parent * py_parent + pz_parent * pz_parent);

    if (p_parent < 1e-15) {
        // Parent at rest — rest frame IS lab frame
        record.secondary_momenta[daughter_index_] = {E_A_rest, px_rest, py_rest, pz_rest};
        record.secondary_momenta[1 - daughter_index_] = {E_B_rest, -px_rest, -py_rest, -pz_rest};
        return;
    }

    double beta = p_parent / E_parent;
    double gamma = E_parent / M_parent;

    // Parent direction unit vector
    double nx = px_parent / p_parent;
    double ny = py_parent / p_parent;
    double nz = pz_parent / p_parent;

    // Decompose rest-frame momentum into parallel and perpendicular
    // to parent direction
    double p_par = px_rest * nx + py_rest * ny + pz_rest * nz;
    double px_perp = px_rest - p_par * nx;
    double py_perp = py_rest - p_par * ny;
    double pz_perp = pz_rest - p_par * nz;

    // Boost parallel component
    double p_par_lab = gamma * (p_par + beta * E_A_rest);
    double E_A_lab = gamma * (E_A_rest + beta * p_par);

    // Lab-frame momentum of daughter A
    double px_A = p_par_lab * nx + px_perp;
    double py_A = p_par_lab * ny + py_perp;
    double pz_A = p_par_lab * nz + pz_perp;

    record.secondary_momenta[daughter_index_] = {E_A_lab, px_A, py_A, pz_A};

    // Daughter B: by momentum conservation
    double E_B_lab = E_parent - E_A_lab;
    record.secondary_momenta[1 - daughter_index_] = {
        E_B_lab,
        px_parent - px_A,
        py_parent - py_A,
        pz_parent - pz_A
    };
}

double Isotropic2BodyChannel::Density(
    std::shared_ptr<siren::detector::DetectorModel const>,
    siren::dataclasses::InteractionRecord const & record) const
{
    // Isotropic in rest-frame solid angle: density = 1/(4*pi)
    return 1.0 / FOUR_PI;
}

} // namespace injection
} // namespace siren
