"""
Standard-model beam decays for dk2nu-driven SIREN injection.

Provides the decay channels that dominate conventional neutrino beams:

    MesonTwoBodyLeptonicDecay   pi+- -> mu nu,  K+- -> mu nu
    MuonThreeBodyDecay          mu-  -> e- nubar_e nu_mu   (and mu+ mirror)

Both are built on the siren.DecayModel authoring base, so injection uses
them like any other decay model: the two-body decay samples isotropically
in the parent rest frame (closure by construction for the SolidAngleRest
measure), and the muon decay samples its polarized matrix element over
the Recursive2Body measure through the three-body machinery implemented
in this module. The module is self-contained: it depends only on siren
and numpy/scipy, not on other process resources.

Width semantics: a dk2nu entry records a decay that already happened into
a known channel, so injection from dk2nu rows fixes the channel and the
branching ratio is carried by the row ensemble. total_width() therefore
returns the PARTIAL width of the implemented channel and
FinalStateProbability is normalized within the channel. Do not add these
models to a collection where the relative rates of several decay channels
of the same parent must come from the widths.

Muon polarization: a beam muon from meson two-body decay is fully
polarized, and the neutrino spectra depend on the spin direction. The
spin axis reaches the model through the pol_x/pol_y/pol_z interaction
parameters that the dk2nu reader attaches to each muon row, using the
same convention as bsim::calcEnuWgt in the dk2nu package: the axis is the
direction of the muon's parent boosted into the muon rest frame, which
equals +spin for mu+ and -spin for mu-. A missing or zero axis gives the
unpolarized average.
"""

_SIREN_RESOURCE_MODULE_ONLY = True

import math

import numpy as np
from scipy import integrate as _integrate

import siren
from siren.dataclasses import Particle

HBAR = 6.582119569e-25       # GeV s
G_F = 1.1663788e-5           # GeV^-2
M_E = 0.51099895e-3          # GeV
M_MU = 0.1056583755          # GeV
M_PI = 0.13957039            # GeV
M_K = 0.493677               # GeV
TAU_PI = 2.6033e-8           # s
TAU_K = 1.2380e-8            # s
TAU_MU = 2.1969811e-6        # s
BR_PI_MUNU = 0.999877
BR_K_MUNU = 0.6356


class MesonTwoBodyLeptonicDecay(siren.DecayModel):
    """Two-body leptonic decay of a charged pion or kaon, M -> mu nu.

    The parent is spin zero, so the decay is isotropic in its rest frame
    and the authoring base derives the sampler and the 1/(4 pi)
    FinalStateProbability from the declared SolidAngleRest measure.
    total_width() is the mu-nu partial width, hbar/tau times the PDG
    branching ratio (given-channel semantics; see the module docstring).
    """

    measure = siren.Measure.SolidAngleRest()
    daughter_index = 1  # the neutrino

    _SPECIES = {
        211: (M_PI, TAU_PI, BR_PI_MUNU),
        321: (M_K, TAU_K, BR_K_MUNU),
    }

    def __init__(self, pdgid_meson):
        siren.DecayModel.__init__(self)
        if abs(pdgid_meson) not in self._SPECIES:
            raise ValueError(
                "MesonTwoBodyLeptonicDecay supports charged pions and kaons "
                "(pdg +-211, +-321), got %r" % pdgid_meson)
        self.m_meson, tau, branching = self._SPECIES[abs(pdgid_meson)]
        self.m_lepton = M_MU
        self.pdgid_meson = pdgid_meson

        # pi+ (211) -> mu+ (-13) nu_mu (14); pi- mirrors to mu- nubar_mu.
        sign = 1 if pdgid_meson > 0 else -1
        self.parent = Particle.ParticleType(pdgid_meson)
        self.daughters = (Particle.ParticleType(-13 * sign),
                          Particle.ParticleType(14 * sign))

        self._width = HBAR / tau * branching

    def total_width(self):
        return self._width

    def differential_width(self, record):
        if int(record.signature.primary_type) != self.pdgid_meson:
            return 0.0
        return self._width / (4.0 * math.pi)

    def density_variables(self):
        return ["cos_theta", "phi"]

    def SecondaryMasses(self, secondary_types):
        return [self.m_lepton, 0.0]

    def SecondaryHelicities(self, record):
        # Angular momentum in the two-body decay of a spin-zero parent
        # forces both daughters into the same helicity: the left-handed
        # nu_mu fixes the mu+ to be left handed as well (the origin of the
        # helicity suppression of pi -> e nu), so a positive meson yields
        # -1/2 for both daughters and a negative meson the mirror. The mu+
        # is thereby fully polarized, which is what the muon-decay model's
        # polarization axis encodes.
        h = -0.5 if self.pdgid_meson > 0 else 0.5
        return [h, h]

    def equal(self, other):
        return self is other


# ---------------------------------------------------------------------------
# Three-body rest-frame machinery
#
# The decay M -> l + nu + phi is parameterized by the rest-frame energies
# (E_nu, E_phi) on the Dalitz region plus the orientation of the decay
# plane, with l the spectator, nu a massless neutrino, and phi a daughter
# that may carry mass. Every function is generic over a physics object
# with m_M, m_l, m_phi, E_nu_max, _matel_sq(E_nu, E_phi), and
# _E_phi_limits(E_nu). The muon decay below uses l = electron, nu = the
# muon-flavor neutrino, phi = the electron-flavor neutrino; the
# DarkNewsTables meson-production models reuse the same machinery with a
# massive mediator as phi.
# ---------------------------------------------------------------------------

def _dalitz_band(m_M, m_l, m_phi, E_nu):
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


def _dalitz_width(physics):
    """Three-body width [GeV]: |M|^2 integrated over the Dalitz region.

    dGamma = |M|^2 dE_nu dE_phi / (64 pi^3 m_M). Integrating the
    implemented matrix element numerically keeps the density that
    FinalStateProbability reports normalized against exactly what the
    sampler draws.
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
        epsrel=1e-3,
    )
    return max(result * prefactor, 0.0)


class _MuonDecayPhysics:
    """Spin-averaged mu -> e nu nubar on the Dalitz plane.

    The matrix element for mu- -> e- nubar_e nu_mu,

        |M|^2 = 64 G_F^2 (p_mu . p_nubar_e)(p_e . p_nu_mu)
              = 32 G_F^2 m_mu E_phi (m_mu^2 - m_e^2 - 2 m_mu E_phi),

    depends only on the rest-frame nubar_e energy; the mu+ mirror is
    identical under CP. The polarization enters as a separate orientation
    factor handled by MuonThreeBodyDecay, so this object stays the
    spin-averaged Dalitz-plane density that the energy sampler and the
    width integral consume.
    """

    def __init__(self):
        self.m_M = M_MU
        self.m_l = M_E
        self.m_phi = 0.0
        self.E_nu_max = (M_MU ** 2 - M_E ** 2) / (2.0 * M_MU)
        self.E_phi_max = (M_MU ** 2 - M_E ** 2) / (2.0 * M_MU)

    def _matel_sq(self, E_nu, E_phi):
        return (32.0 * G_F ** 2 * self.m_M * E_phi
                * (self.m_M ** 2 - self.m_l ** 2 - 2.0 * self.m_M * E_phi))

    def _E_phi_limits(self, E_nu):
        """Kinematically allowed E_phi range at fixed E_nu (rest frame)."""
        return _dalitz_band(self.m_M, self.m_l, self.m_phi, E_nu)

    def total_width(self):
        """Three-body width [GeV] from the shared Dalitz integral."""
        return _dalitz_width(self)


def _find_max_weight(physics):
    """Deterministic bound on |M|^2 * band(E_nu), the acceptance weight of
    the rejection sampler in _sample_energies, over a Dalitz-region grid."""
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


def _sample_energies(physics, max_weight, random):
    """Rejection-sample (E_nu, E_phi) from |M|^2 on the Dalitz region.

    The proposal draws E_nu uniformly and then E_phi uniformly within its
    E_nu-dependent kinematic band, so the proposal density carries a
    1/band(E_nu) factor. The acceptance weight is |M|^2 * band(E_nu): the
    band factors cancel and the accepted density is proportional to |M|^2
    alone, matching the Dalitz-plane density that FinalStateProbability
    and total_width integrate. max_weight bounds |M|^2 * band over the
    region (see _find_max_weight).
    """
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
        u = random.Uniform(0.0, max_weight)
        if u <= weight:
            return E_nu, E_phi
    lims = d._E_phi_limits(d.E_nu_max * 0.5)
    if lims[0] is None:
        return d.E_nu_max * 0.5, d.m_phi
    return d.E_nu_max * 0.5, 0.5 * (lims[0] + lims[1])


def _build_rest_momenta(physics, E_nu, E_phi, random):
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


def _boost_four_vector(P_parent, P_rest, M_parent):
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


def _boost_to_rest(P_parent, P_lab):
    """Boost a lab-frame four-vector into the parent rest frame.

    The spatial part follows the same algebra as bsim::calcEnuWgt, so the
    boosted vectors live in the frame whose spatial basis the dk2nu
    polarization axis is expressed in (pure boost, no rotation).
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


def _final_state_probability(physics, total_width, pdgid_nu, pdgid_phi,
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


class MuonThreeBodyDecay(siren.DecayModel):
    """Three-body muon decay mu -> e nu nubar with beam polarization.

    Declared over the Recursive2Body measure with the electron as the
    spectator and the two neutrinos as the pair, so the directed
    three-body channel (channels.toward_3body) can aim the muon-flavor
    neutrino at a detector.

    The density is the tree-level matrix element at fixed muon spin.
    Writing a-hat for the spin axis in the dk2nu convention (the
    parent-of-muon direction boosted into the muon rest frame, equal to
    +spin for mu+ and -spin for mu-), both charges share one form,

        |M|^2 = 64 G_F^2 [(p_mu + m_mu a) . p_e-type](p_e . p_mu-type),

    with a = (0, a-hat) in the muon rest frame: the unpolarized element
    times (E - a-hat . p)/E of the electron-flavor neutrino. The implied
    single-particle spectra are the classic polarized Michel forms,
    x^2 [(3 - 2x) - (1 - 2x) cos theta] for the muon-flavor neutrino and
    12 x^2 (1 - x)(1 - cos theta) for the electron-flavor one with theta
    measured from a-hat, exactly the correction bsim::calcEnuWgt applies.

    The axis arrives per event through the pol_x/pol_y/pol_z interaction
    parameters (detector-frame components, filled by the dk2nu reader);
    its magnitude is the polarization degree, and a missing or zero axis
    gives the unpolarized average. total_width() integrates the
    spin-averaged matrix element numerically -- polarization does not
    change the width -- which keeps FinalStateProbability normalized
    against exactly what sample() draws; it lands within QED radiative
    corrections (about half a percent) of hbar/tau_mu.
    """

    measure = siren.Measure.Recursive2Body(
        spectator=0, pair_first=1, pair_second=2)

    def __init__(self, pdgid_muon):
        siren.DecayModel.__init__(self)
        if abs(pdgid_muon) != 13:
            raise ValueError(
                "MuonThreeBodyDecay supports pdg +-13, got %r" % pdgid_muon)
        # mu- (13) -> e- (11) nu_mu (14) nubar_e (-12); mu+ mirrors.
        sign = 1 if pdgid_muon > 0 else -1
        self.pdgid_muon = pdgid_muon
        self.pdgid_electron = 11 * sign
        self.pdgid_nu_mu = 14 * sign
        self.pdgid_nu_e = -12 * sign

        self.parent = Particle.ParticleType(pdgid_muon)
        self.daughters = (Particle.ParticleType(self.pdgid_electron),
                          Particle.ParticleType(self.pdgid_nu_mu),
                          Particle.ParticleType(self.pdgid_nu_e))

        self._physics = _MuonDecayPhysics()
        self._total = self._physics.total_width()
        self._max_matel = _find_max_weight(self._physics)

    @staticmethod
    def _polarization(parameters):
        """Spin axis from a record's interaction parameters, or None."""
        try:
            axis = np.array([parameters["pol_x"], parameters["pol_y"],
                             parameters["pol_z"]], dtype=float)
        except KeyError:
            return None
        if not np.any(np.abs(axis) > 0):
            return None
        return axis

    @staticmethod
    def _polarization_factor(pol, P_phi_rest):
        """(E - a . p)/E of the electron-flavor neutrino in the rest frame."""
        E = P_phi_rest[0]
        if E <= 0:
            return 0.0
        return max(1.0 - float(np.dot(pol, P_phi_rest[1:])) / E, 0.0)

    def total_width(self):
        return self._total

    def differential_width(self, record):
        if int(record.signature.primary_type) != self.pdgid_muon:
            return 0.0
        density = _final_state_probability(
            self._physics, self._total,
            self.pdgid_nu_mu, self.pdgid_nu_e, record)
        if density <= 0.0:
            return 0.0
        pol = self._polarization(record.interaction_parameters)
        if pol is not None:
            idx = -1
            for i, stype in enumerate(record.signature.secondary_types):
                if int(stype) == self.pdgid_nu_e:
                    idx = i
                    break
            if idx < 0:
                return 0.0
            P_phi_rest = _boost_to_rest(
                record.primary_momentum, record.secondary_momenta[idx])
            density *= self._polarization_factor(pol, P_phi_rest)
        return density * self._total

    def density_variables(self):
        return ["s_pair", "cos_theta_sub"]

    def sample(self, record, random):
        e_nu_rf, e_phi_rf = _sample_energies(
            self._physics, self._max_matel, random)
        pol = self._polarization(record.interaction_parameters)
        parent = np.asarray(record.primary_momentum, dtype=float)
        # The spin correlation reweights the orientation alone: the
        # builder orients the configuration Haar-uniformly, and the factor
        # (E - a . p)/E of the electron-flavor neutrino averages to one
        # over orientations at fixed energies. Accepting the orientation
        # against that factor therefore leaves the (E_nu, E_phi) marginal
        # of the energy sampler untouched and produces the joint polarized
        # density that differential_width integrates.
        bound = 1.0 + (float(np.linalg.norm(pol)) if pol is not None else 0.0)
        for _ in range(1000):
            p_nu, p_e, p_phi = _build_rest_momenta(
                self._physics, e_nu_rf, e_phi_rf, random)
            if pol is None:
                break
            f = self._polarization_factor(pol, p_phi)
            if random.Uniform(0.0, bound) <= f:
                break
        for sec in record.get_secondary_particle_records():
            pid = int(sec.type)
            if pid == self.pdgid_electron:
                sec.four_momentum = _boost_four_vector(parent, p_e, M_MU)
                sec.mass = M_E
            elif pid == self.pdgid_nu_mu:
                sec.four_momentum = _boost_four_vector(parent, p_nu, M_MU)
                sec.mass = 0.0
            elif pid == self.pdgid_nu_e:
                sec.four_momentum = _boost_four_vector(parent, p_phi, M_MU)
                sec.mass = 0.0

    def SecondaryMasses(self, secondary_types):
        return [M_E, 0.0, 0.0]

    def SecondaryHelicities(self, record):
        # Massless-limit helicities: each lepton left handed (-1/2), each
        # antilepton right handed (+1/2). The muon-flavor neutrino's value
        # matters downstream: DarkNews multiplies it into the helicity of
        # the upscattered state.
        return [-0.5 if pid > 0 else 0.5
                for pid in (self.pdgid_electron, self.pdgid_nu_mu,
                            self.pdgid_nu_e)]

    def equal(self, other):
        return self is other
