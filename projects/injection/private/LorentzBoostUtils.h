#pragma once
#ifndef SIREN_Injection_LorentzBoostUtils_H
#define SIREN_Injection_LorentzBoostUtils_H

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include <rk/geom3.hh>
#include <rk/rk.hh>

namespace siren {
namespace injection {
namespace detail {

inline std::array<double, 4> BoostRestFrameToLab(
    double parent_energy,
    double parent_px, double parent_py, double parent_pz,
    double rest_energy,
    double rest_px, double rest_py, double rest_pz)
{
    double parent_p2 = parent_px * parent_px
                     + parent_py * parent_py
                     + parent_pz * parent_pz;
    if (parent_p2 < 1e-30) {
        return {rest_energy, rest_px, rest_py, rest_pz};
    }

    rk::P4 parent(
        parent_energy, geom3::Vector3(parent_px, parent_py, parent_pz));
    double rest_p2 = rest_px * rest_px + rest_py * rest_py + rest_pz * rest_pz;
    double rest_mass2 = rest_energy * rest_energy - rest_p2;
    double scale = rest_energy * rest_energy + rest_p2 + 1.0;
    if (rest_mass2 < -1e-12 * scale) {
        throw std::runtime_error(
            "Cannot boost a spacelike rest-frame four-momentum");
    }
    if (rest_mass2 < 0.0 && rest_mass2 > -1e-12 * scale) rest_mass2 = 0.0;
    rk::P4 particle(
        geom3::Vector3(rest_px, rest_py, rest_pz),
        std::sqrt(std::max(0.0, rest_mass2)), rest_energy < 0.0);
    rk::P4 lab = particle.boost(parent.labBoost());
    return {lab.e(), lab.px(), lab.py(), lab.pz()};
}

} // namespace detail
} // namespace injection
} // namespace siren

#endif // SIREN_Injection_LorentzBoostUtils_H
