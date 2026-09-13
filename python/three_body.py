"""Three-body decay utilities for M -> l + nu + phi, with massless nu.

Masses, energies and four-momenta (E, px, py, pz) are in GeV. The energy
sampler and probability use a matrix element independent of orientation;
polarization or other angular correlations belong in the model.
The supported masses satisfy m_M > m_l + m_phi > 0, with m_l, m_phi >= 0.

Functions accepting ``physics`` use its ``m_M``, ``m_l``, ``m_phi``,
``E_nu_max``, ``_matel_sq(E_nu, E_phi)`` and ``_E_phi_limits(E_nu)``.
The latter returns (lower, upper) energies or (None, None) outside support.
No inheritance is required. See docs/three_body.md for the measure and
normalization conventions.
"""

import math

import numpy as np
from scipy import integrate as _integrate

__all__ = [
    "dalitz_band", "dalitz_width", "find_max_weight", "sample_energies",
    "build_rest_momenta", "boost_four_vector", "boost_to_rest",
    "final_state_probability",
]


def dalitz_band(m_M, m_l, m_phi, E_nu):
    """Kinematically allowed E_phi range at fixed E_nu (parent rest frame).

    At fixed E_nu the (l, phi) system recoils with momentum E_nu and
    invariant mass squared m_M^2 - 2 m_M E_nu; the band follows from the
    two-body decay of that system boosted back to the parent frame.
    Returns (None, None) outside the Dalitz region.
    """
    E_nu_max = (m_M ** 2 - (m_l + m_phi) ** 2) / (2.0 * m_M)
    if E_nu < 0.0 or E_nu > E_nu_max:
        return None, None

    M_lph2 = m_M ** 2 - 2.0 * m_M * E_nu
    if M_lph2 < (m_l + m_phi) ** 2:
        return None, None

    M_lph = math.sqrt(M_lph2)
    lam = ((M_lph2 - (m_l + m_phi) ** 2)
           * (M_lph2 - (m_l - m_phi) ** 2))
    if lam < 0.0:
        return None, None

    pstar = math.sqrt(lam) / (2.0 * M_lph)
    Estar = (M_lph2 - m_l ** 2 + m_phi ** 2) / (2.0 * M_lph)

    gamma = (m_M - E_nu) / M_lph
    beta_gamma = E_nu / M_lph

    E_phi_hi = gamma * Estar + beta_gamma * pstar
    E_phi_lo = max(gamma * Estar - beta_gamma * pstar, m_phi)

    if m_M - E_nu - E_phi_lo < m_l:
        return None, None
    return E_phi_lo, E_phi_hi


def dalitz_width(physics):
    """Three-body width [GeV]: |M|^2 integrated over the Dalitz region.

    dGamma = |M|^2 dE_nu dE_phi / (64 pi^3 m_M). The integral uses the
    same matrix element as the energy sampler; its relative error also
    affects the normalization of final_state_probability.

    Quadrature requests a relative tolerance of 1e-6 with no absolute
    error floor. A fixed absolute tolerance would relax the accuracy
    requirement when small couplings scale down the integrand.
    """
    d = physics
    prefactor = 1.0 / (64.0 * math.pi ** 3 * d.m_M)

    def integrand(E_phi, E_nu):
        lims = d._E_phi_limits(E_nu)
        if lims[0] is None:
            return 0.0
        if E_phi < lims[0] or E_phi > lims[1]:
            return 0.0
        if d.m_M - E_nu - E_phi < d.m_l:
            return 0.0
        return d._matel_sq(E_nu, E_phi)

    result, _ = _integrate.dblquad(
        integrand,
        0.0, d.E_nu_max,
        lambda E_nu: (d._E_phi_limits(E_nu)[0] or 0.0),
        lambda E_nu: (d._E_phi_limits(E_nu)[1] or 0.0),
        epsabs=0.0,
        epsrel=1e-6,
    )
    return max(result * prefactor, 0.0)


def find_max_weight(physics):
    """Estimate max(|M|^2 * band) on a 400-by-200 grid, with 20% padding.

    This is not a guaranteed bound for an arbitrary matrix element. Supply
    a model-specific bound when narrow peaks are not resolved by the grid.
    """
    d = physics
    best = 0.0
    for e_nu in np.linspace(0.0, d.E_nu_max, 400):
        lims = d._E_phi_limits(e_nu)
        if lims[0] is None:
            continue
        band = lims[1] - lims[0]
        for e_phi in np.linspace(lims[0], lims[1], 200):
            val = d._matel_sq(e_nu, e_phi) * band
            if val > best:
                best = val
    return best * 1.2


def sample_energies(physics, max_weight, random):
    """Rejection-sample (E_nu, E_phi) from |M|^2 on the Dalitz region.

    The proposal draws E_nu uniformly and then E_phi uniformly within its
    E_nu-dependent kinematic band, so the proposal density carries a
    1/band(E_nu) factor. The acceptance weight is |M|^2 * band(E_nu): the
    band factors cancel and the accepted density is proportional to |M|^2
    alone, matching the Dalitz-plane density that FinalStateProbability
    and total_width integrate. max_weight bounds |M|^2 * band over the
    region (see find_max_weight). A nonpositive/nonfinite bound raises
    ValueError. Invalid weights, a bound overrun, or 10,000 rejected proposals
    raise RuntimeError; no substitute event is returned.
    """
    if not math.isfinite(max_weight) or max_weight <= 0:
        raise ValueError("sample_energies requires a finite positive max_weight")
    d = physics
    for _ in range(10000):
        E_nu = random.Uniform(0.0, d.E_nu_max)
        lims = d._E_phi_limits(E_nu)
        if lims[0] is None:
            continue
        E_phi = random.Uniform(lims[0], lims[1])
        if d.m_M - E_nu - E_phi < d.m_l:
            continue
        weight = d._matel_sq(E_nu, E_phi) * (lims[1] - lims[0])
        if not math.isfinite(weight) or weight < 0:
            raise RuntimeError("sample_energies encountered an invalid acceptance weight")
        if weight > max_weight:
            raise RuntimeError("sample_energies acceptance weight exceeds max_weight")
        u = random.Uniform(0.0, max_weight)
        if weight > 0 and u <= weight:
            return E_nu, E_phi
    raise RuntimeError("sample_energies exhausted 10000 proposals")


def build_rest_momenta(physics, E_nu, E_phi, random):
    """Construct the three rest-frame four-momenta at a Haar-uniform
    orientation: the nu direction is isotropic, the phi direction sits at
    the opening angle fixed by momentum balance with a uniform azimuth
    around the nu axis, and the lepton carries the balancing momentum.
    The nu is massless; the phi momentum comes from physics.m_phi.

    Returns (P_nu, P_l, P_phi) as four-vectors in the parent rest frame.
    """
    m_M, m_l = physics.m_M, physics.m_l
    E_l = m_M - E_nu - E_phi

    p_nu = E_nu
    p_phi = math.sqrt(max(E_phi ** 2 - physics.m_phi ** 2, 0.0))

    cos_nu = random.Uniform(-1.0, 1.0)
    phi_nu = random.Uniform(0.0, 2.0 * math.pi)
    sin_nu = math.sqrt(max(1.0 - cos_nu ** 2, 0.0))
    nu_dir = np.array([sin_nu * math.cos(phi_nu),
                       sin_nu * math.sin(phi_nu),
                       cos_nu])

    # Opening angle between the nu and phi momenta from momentum balance:
    # p_l = -(p_nu + p_phi) gives
    # E_l^2 - m_l^2 = p_nu^2 + p_phi^2 + 2 p_nu p_phi cos(theta).
    if p_nu > 0 and p_phi > 0:
        cos_open = (E_l ** 2 - m_l ** 2 - p_nu ** 2 - p_phi ** 2) \
            / (2.0 * p_nu * p_phi)
        cos_open = max(-1.0, min(1.0, cos_open))
    else:
        cos_open = 0.0
    sin_open = math.sqrt(max(1.0 - cos_open ** 2, 0.0))

    perp1 = np.cross(nu_dir, np.array([0.0, 0.0, 1.0]))
    if np.linalg.norm(perp1) < 1e-10:
        perp1 = np.cross(nu_dir, np.array([0.0, 1.0, 0.0]))
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(nu_dir, perp1)

    azimuth = random.Uniform(0.0, 2.0 * math.pi)
    phi_dir = (cos_open * nu_dir
               + sin_open * math.cos(azimuth) * perp1
               + sin_open * math.sin(azimuth) * perp2)

    p_nu_vec = p_nu * nu_dir
    p_phi_vec = p_phi * phi_dir
    p_l_vec = -(p_nu_vec + p_phi_vec)

    return (np.array([E_nu, *p_nu_vec]),
            np.array([E_l, *p_l_vec]),
            np.array([E_phi, *p_phi_vec]))


def boost_four_vector(P_parent, P_rest, M_parent):
    """Boost a rest-frame four-vector to the frame where the parent has
    four-momentum P_parent. Summing the boosted daughters reproduces the
    parent four-momentum exactly."""
    E_parent = P_parent[0]
    p_parent = np.asarray(P_parent[1:], dtype=float)
    p_mag = float(np.linalg.norm(p_parent))
    if p_mag < 1e-12 or M_parent < 1e-12:
        return np.array(P_rest, dtype=float)

    beta = p_mag / E_parent
    gamma = E_parent / M_parent
    beta_hat = p_parent / p_mag

    E_rest = P_rest[0]
    p_rest = np.asarray(P_rest[1:], dtype=float)
    p_par = float(np.dot(p_rest, beta_hat))
    p_perp = p_rest - p_par * beta_hat

    E_lab = gamma * (E_rest + beta * p_par)
    p_par_lab = gamma * (p_par + beta * E_rest)
    return np.array([E_lab, *(p_par_lab * beta_hat + p_perp)])


def boost_to_rest(P_parent, P_lab):
    """Boost a lab-frame four-vector into the parent rest frame.

    Applies a pure boost without rotating the spatial basis, matching
    bsim::calcEnuWgt. BeamDecays' muon density uses this basis to compare
    the daughter momentum with the dk2nu polarization axis.
    """
    P_parent = np.asarray(P_parent, dtype=float)
    P_lab = np.asarray(P_lab, dtype=float)
    p_par = P_parent[1:]
    p_mag = float(np.linalg.norm(p_par))
    if p_mag < 1e-12:
        return P_lab.copy()
    E_par = P_parent[0]
    mass_sq = max(E_par * E_par - p_mag * p_mag, 0.0)
    mass = math.sqrt(mass_sq)
    if mass < 1e-12:
        return P_lab.copy()
    gamma = E_par / mass
    beta = p_par / E_par
    beta_dot_p = float(np.dot(beta, P_lab[1:]))
    partial = P_lab[0] - gamma * beta_dot_p / (gamma + 1.0)
    p_rest = P_lab[1:] - beta * gamma * partial
    E_rest = gamma * (P_lab[0] - beta_dot_p)
    return np.array([E_rest, *p_rest])


def final_state_probability(physics, total_width, pdgid_nu, pdgid_phi,
                             record):
    """Spin-averaged density over the Recursive2Body measure.

    Returns dGamma / (ds_pair dOmega_pair dOmega_sub Gamma_total) with the
    factorization M -> l(spectator) + (nu phi)(pair). The rest-frame
    energies are recovered from the record by inverse boost; the pair
    direction is isotropic for the spin-averaged element, so the density
    depends on (s_pair, cos_theta_sub) through the matrix element and the
    Recursive2Body-to-Dalitz Jacobian.
    """
    if total_width <= 0:
        return 0.0

    d = physics
    m_M, m_l, m_phi = d.m_M, d.m_l, d.m_phi

    P_parent = np.asarray(record.primary_momentum, dtype=float)
    E_parent = P_parent[0]
    p_parent = float(np.linalg.norm(P_parent[1:]))

    idx_nu = idx_phi = -1
    for idx, stype in enumerate(record.signature.secondary_types):
        pid = int(stype)
        if pid == pdgid_nu:
            idx_nu = idx
        elif pid == pdgid_phi:
            idx_phi = idx
    if idx_nu < 0 or idx_phi < 0:
        return 0.0

    P_nu = np.asarray(record.secondary_momenta[idx_nu], dtype=float)
    P_phi = np.asarray(record.secondary_momenta[idx_phi], dtype=float)

    if p_parent < 1e-12:
        E_nu_rf = P_nu[0]
        E_phi_rf = P_phi[0]
    else:
        parent_mass_sq = max(E_parent ** 2 - p_parent ** 2, 0.0)
        parent_mass = math.sqrt(parent_mass_sq) if parent_mass_sq > 0 else m_M
        beta = p_parent / E_parent
        gamma = E_parent / parent_mass
        beta_hat = P_parent[1:] / p_parent
        E_nu_rf = gamma * (P_nu[0] - beta * float(np.dot(P_nu[1:], beta_hat)))
        E_phi_rf = gamma * (P_phi[0]
                            - beta * float(np.dot(P_phi[1:], beta_hat)))

    if E_nu_rf < 0 or E_phi_rf < m_phi:
        return 0.0

    matel_sq = d._matel_sq(E_nu_rf, E_phi_rf)
    if matel_sq <= 0:
        return 0.0

    # Pair invariant mass from the spectator energy:
    # s_pair = (P - p_l)^2 = m_M^2 + m_l^2 - 2 m_M E_l.
    E_l_rf = m_M - E_nu_rf - E_phi_rf
    s_pair = m_M ** 2 + m_l ** 2 - 2.0 * m_M * E_l_rf
    if s_pair <= 0:
        return 0.0
    sqrt_s = math.sqrt(s_pair)

    # Recursive2Body-to-Dalitz Jacobian |ds_13/dcos_sub|
    # = 2 * p_spectator_in_pair_frame * p_sub_in_pair_frame.
    lam_parent = ((m_M ** 2 - (m_l + sqrt_s) ** 2)
                  * (m_M ** 2 - (m_l - sqrt_s) ** 2))
    lam_pair = (s_pair - m_phi ** 2) * (s_pair - m_phi ** 2)
    if lam_parent <= 0 or lam_pair <= 0:
        return 0.0
    p_spec_pair = math.sqrt(lam_parent) / (2.0 * sqrt_s)
    p_sub_pair = math.sqrt(lam_pair) / (2.0 * sqrt_s)
    jacobian = 2.0 * p_spec_pair * p_sub_pair

    # Dalitz density |M|^2 / (256 pi^3 m_M^3), spread over the isotropic
    # pair direction (4 pi) and sub azimuth (2 pi).
    prefactor = jacobian / (2048.0 * math.pi ** 5 * m_M ** 3)

    return matel_sq * prefactor / total_width
