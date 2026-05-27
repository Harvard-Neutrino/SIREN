"""
Meson three-body decay and scalar/pseudoscalar production for SIREN.

Implements charged meson three-body decays pi/K -> l nu phi from
Dutta, Kim, Thompson, Thornton, Van de Water, PRL 129, 111803 (2022),
using the Carlson-Rislow matrix element (Phys.Rev.D 86, 035013, 2012).

Also provides MesonSimpleDecay (pi -> mu nu) for external-distribution
pipelines, and build_phi_flux() for constructing the scalar mediator
flux at the detector.

All classes are self-contained with no DarkNews imports.
"""

import os
import math
import numpy as np
import scipy.integrate as _integrate

from siren.interactions import Decay as _Decay
from siren import dataclasses
from siren.dataclasses import Particle
from siren.injection import PhaseSpaceConvention as _PhaseSpaceConvention
from siren.injection import PhaseSpaceTopology as _Topology
from siren.injection import PhaseSpaceMeasure as _Measure


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_GF = 1.16638e-5       # Fermi constant [GeV^-2]
_FPI = 0.1307           # pion decay constant [GeV]
_FK = 0.1598            # kaon decay constant [GeV]
_VUD = 0.9737           # CKM |V_ud|
_VUS = 0.2245           # CKM |V_us|
_ALPHA_EM = 1.0 / 137.036

_M_PI = 0.13957039      # charged pion mass [GeV]
_M_K = 0.49368          # charged kaon mass [GeV]
_M_MU = 0.10565837      # muon mass [GeV]
_M_E = 0.000511         # electron mass [GeV]

# PDG IDs
_PDGID_PIPLUS = 211
_PDGID_KPLUS = 321
_PDGID_MUPLUS = -13
_PDGID_EPLUS = -11
_PDGID_NUMU = 14
_PDGID_NUE = 12


def _primary_type(arg):
    if isinstance(arg, dataclasses.InteractionRecord):
        return arg.signature.primary_type
    return arg


# ---------------------------------------------------------------------------
# Kinematic helpers
# ---------------------------------------------------------------------------

def _two_body_p_cm(M, m1, m2):
    arg = (M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2)
    if arg <= 0.0:
        return 0.0
    return math.sqrt(arg) / (2.0 * M)


def _boost_to_lab(P_parent, p_cm, cos_theta, phi, m_daughter):
    E_parent = P_parent[0]
    p_parent = P_parent[1:]
    M_parent_sq = max(E_parent**2 - np.dot(p_parent, p_parent), 0.0)
    M_parent = math.sqrt(M_parent_sq)

    sin_theta = math.sqrt(max(1.0 - cos_theta**2, 0.0))
    E_rf = math.sqrt(p_cm**2 + m_daughter**2)
    p_rf = np.array([
        p_cm * sin_theta * math.cos(phi),
        p_cm * sin_theta * math.sin(phi),
        p_cm * cos_theta,
    ])

    p_mag = np.linalg.norm(p_parent)
    if p_mag < 1e-12 or M_parent < 1e-12:
        return np.array([E_rf, *p_rf])

    beta = p_mag / E_parent
    gamma = E_parent / M_parent
    beta_hat = p_parent / p_mag

    p_par = np.dot(p_rf, beta_hat)
    p_perp = p_rf - p_par * beta_hat

    E_lab = gamma * (E_rf + beta * p_par)
    p_par_lab = gamma * (p_par + beta * E_rf)
    return np.array([E_lab, *(p_par_lab * beta_hat + p_perp)])


# ---------------------------------------------------------------------------
# Meson decay constants and CKM elements
# ---------------------------------------------------------------------------

def _meson_params(m_meson):
    """Return (f_M, V_Mq) for a given meson mass."""
    if abs(m_meson - _M_PI) < 0.01:
        return _FPI, _VUD
    elif abs(m_meson - _M_K) < 0.01:
        return _FK, _VUS
    return _FPI, _VUD


# ===================================================================
#  Carlson-Rislow matrix element for M -> l nu phi
# ===================================================================

def _matel_sq_scalar(E_nu, E_phi, m_M, m_l, m_phi, C2):
    """
    Spin-summed |M|^2 for scalar coupling.
    Carlson & Rislow, Phys.Rev.D 86, 035013 (2012), Eq. 25.
    """
    t = m_M**2 - 2.0 * m_M * E_nu
    u = m_M**2 + m_phi**2 - 2.0 * m_M * E_phi
    D = t - m_l**2
    if D <= 0.0:
        return 0.0
    common = ((t + u - m_phi**2) * t * D
              - (t**2 - m_l**2 * m_M**2) * (t + m_l**2 - m_phi**2))
    mass_term = m_l**2 * t * (m_M**2 - t)
    T_S = 8.0 * (common + 2.0 * mass_term)
    return max(C2 * T_S / D**2, 0.0)


def _matel_sq_pseudo(E_nu, E_phi, m_M, m_l, m_phi, C2):
    """
    Spin-summed |M|^2 for pseudoscalar coupling.
    Same as scalar but with sign flip on the mass term.
    """
    t = m_M**2 - 2.0 * m_M * E_nu
    u = m_M**2 + m_phi**2 - 2.0 * m_M * E_phi
    D = t - m_l**2
    if D <= 0.0:
        return 0.0
    common = ((t + u - m_phi**2) * t * D
              - (t**2 - m_l**2 * m_M**2) * (t + m_l**2 - m_phi**2))
    mass_term = m_l**2 * t * (m_M**2 - t)
    T_P = 8.0 * (common - 2.0 * mass_term)
    return max(C2 * T_P / D**2, 0.0)


# ===================================================================
#  MesonThreeBodyDecay  --  pi/K -> l nu phi
# ===================================================================

class MesonThreeBodyDecay:
    """
    Three-body charged meson decay M -> l nu phi (scalar or pseudoscalar).

    This is a pure-physics class (no SIREN C++ base).  It can be used
    standalone for rate calculations, or wrapped by build_phi_flux() to
    construct the mediator flux at a detector.

    Parameters
    ----------
    m_meson : float
        Parent meson mass [GeV].
    m_lepton : float
        Final-state charged lepton mass [GeV].
    m_mediator : float
        Scalar/pseudoscalar mediator mass [GeV].
    g_mu : float
        Yukawa coupling of the mediator to the muon.
    mediator_type : str
        "scalar" or "pseudoscalar".
    """

    def __init__(self, m_meson, m_lepton, m_mediator, g_mu, mediator_type="scalar"):
        self.m_M = m_meson
        self.m_l = m_lepton
        self.m_phi = m_mediator
        self.g_mu = g_mu
        self.mediator_type = mediator_type

        f_M, V_Mq = _meson_params(m_meson)
        self.f_M = f_M
        self.V_Mq = V_Mq
        self._C2 = (_GF * f_M * V_Mq * g_mu)**2 / 2.0

        if m_meson < m_lepton + m_mediator:
            raise ValueError(
                f"Kinematically forbidden: m_M={m_meson:.4f} < "
                f"m_l+m_phi={m_lepton + m_mediator:.4f} GeV"
            )

        self.E_nu_max = (m_meson**2 - (m_lepton + m_mediator)**2) / (2.0 * m_meson)
        self.E_phi_max = (m_meson**2 + m_mediator**2 - m_lepton**2) / (2.0 * m_meson)

    def _matel_sq(self, E_nu, E_phi):
        if self.mediator_type == "scalar":
            return _matel_sq_scalar(E_nu, E_phi, self.m_M, self.m_l, self.m_phi, self._C2)
        elif self.mediator_type == "pseudoscalar":
            return _matel_sq_pseudo(E_nu, E_phi, self.m_M, self.m_l, self.m_phi, self._C2)
        raise ValueError(f"Unknown mediator_type: {self.mediator_type}")

    def _E_phi_limits(self, E_nu):
        """Kinematically allowed E_phi range for a given E_nu (meson rest frame)."""
        m_M = self.m_M
        m_l = self.m_l
        m_phi = self.m_phi

        if E_nu < 0.0 or E_nu > self.E_nu_max:
            return None, None

        M_lph2 = m_M**2 - 2.0 * m_M * E_nu
        if M_lph2 < (m_l + m_phi)**2:
            return None, None

        M_lph = math.sqrt(M_lph2)
        lam = (M_lph2 - (m_l + m_phi)**2) * (M_lph2 - (m_l - m_phi)**2)
        if lam < 0.0:
            return None, None

        pstar = math.sqrt(lam) / (2.0 * M_lph)
        Estar = (M_lph2 - m_l**2 + m_phi**2) / (2.0 * M_lph)

        gamma = (m_M - E_nu) / M_lph
        beta_gamma = E_nu / M_lph

        E_phi_hi = gamma * Estar + beta_gamma * pstar
        E_phi_lo = max(gamma * Estar - beta_gamma * pstar, m_phi)

        if m_M - E_nu - E_phi_lo < m_l:
            return None, None

        return E_phi_lo, E_phi_hi

    def total_width(self):
        """Total three-body decay width [GeV]."""
        prefactor = 1.0 / (64.0 * math.pi**3 * self.m_M)

        def integrand(E_phi, E_nu):
            lims = self._E_phi_limits(E_nu)
            if lims[0] is None:
                return 0.0
            if E_phi < lims[0] or E_phi > lims[1]:
                return 0.0
            if self.m_M - E_nu - E_phi < self.m_l:
                return 0.0
            return self._matel_sq(E_nu, E_phi)

        result, _ = _integrate.dblquad(
            integrand,
            0.0, self.E_nu_max,
            lambda Enu: (self._E_phi_limits(Enu)[0] or 0.0),
            lambda Enu: (self._E_phi_limits(Enu)[1] or 0.0),
            epsrel=1e-3,
        )
        return max(result * prefactor, 0.0)

    def differential_decay_rate(self, E_phi_vals):
        """
        dGamma/dE_phi integrated over E_nu, evaluated at each E_phi in the array.
        Returns array in units of GeV^{-1} (partial width per GeV of E_phi).
        """
        m_M = self.m_M
        m_l = self.m_l
        m_phi = self.m_phi
        prefactor = 1.0 / (64.0 * math.pi**3 * m_M)

        results = np.zeros_like(E_phi_vals, dtype=float)
        E_phi_min_global = m_phi
        E_phi_max_global = self.E_phi_max

        if E_phi_max_global <= E_phi_min_global:
            return results

        for i, Ep in enumerate(E_phi_vals):
            if Ep <= E_phi_min_global or Ep >= E_phi_max_global:
                continue

            u_val = m_M**2 + m_phi**2 - 2.0 * m_M * Ep
            q2 = u_val
            if q2 <= m_l**2:
                continue

            mq = math.sqrt(q2)
            E_mu_star = (q2 + m_l**2) / (2.0 * mq)
            p_mu_star = math.sqrt(max(E_mu_star**2 - m_l**2, 0.0))

            E_q = m_M - Ep
            p_q = math.sqrt(max(E_q**2 - q2, 0.0))
            gamma_boost = E_q / mq
            beta_gamma = p_q / mq

            E_mu_max = gamma_boost * E_mu_star + beta_gamma * p_mu_star
            E_mu_min = max(gamma_boost * E_mu_star - beta_gamma * p_mu_star, m_l)

            if E_mu_max <= E_mu_min:
                continue

            def integrand(Emu, _Ep=Ep):
                Enu = m_M - _Ep - Emu
                if Enu < 0.0:
                    return 0.0
                return self._matel_sq(Enu, _Ep)

            try:
                val, _ = _integrate.quad(integrand, E_mu_min, E_mu_max,
                                         limit=40, epsrel=1e-3)
            except Exception:
                val = 0.0

            results[i] = max(val, 0.0)

        return prefactor * results


# ===================================================================
#  MesonSimpleDecay  --  pi+ -> mu+ nu_mu  (SM two-body)
# ===================================================================

class MesonSimpleDecay(_Decay):
    """
    SM two-body decay pi+ -> mu+ nu_mu.
    Width: Gamma = G_F^2 f_pi^2 |V_ud|^2 m_pi m_mu^2 (1 - m_mu^2/m_pi^2)^2 / (8 pi)
    """

    def __init__(
        self,
        m_meson=_M_PI,
        m_lepton=_M_MU,
        *,
        pdgid_meson=_PDGID_PIPLUS,
        pdgid_lepton=_PDGID_MUPLUS,
        pdgid_neutrino=_PDGID_NUMU,
        table_dir=None,
    ):
        _Decay.__init__(self)

        self.m_meson = m_meson
        self.m_lepton = m_lepton
        self.m_nu = 0.0
        self.pdgid_meson = pdgid_meson
        self.pdgid_lepton = pdgid_lepton
        self.pdgid_neutrino = pdgid_neutrino

        f_M, V_Mq = _meson_params(m_meson)
        self.f_M = f_M
        self.V_Mq = V_Mq

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)
        self._total_width = self._compute_width()

    def _compute_width(self):
        r = (self.m_lepton / self.m_meson)**2
        return (_GF**2 * self.f_M**2 * self.V_Mq**2
                * self.m_meson * self.m_lepton**2
                * (1.0 - r)**2 / (8.0 * math.pi))

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_meson)
        sig.target_type = Particle.ParticleType.Decay
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_lepton),
            Particle.ParticleType(self.pdgid_neutrino),
        ]
        return [sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if int(primary_type) == self.pdgid_meson:
            return self.GetPossibleSignatures()
        return []

    def TotalDecayWidthAllFinalStates(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def TotalDecayWidth(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def DifferentialDecayWidth(self, record):
        if int(record.signature.primary_type) != self.pdgid_meson:
            return 0.0
        return self._total_width / (4.0 * math.pi)

    def FinalStateProbability(self, record):
        if self._total_width <= 0.0:
            return 0.0
        return self.DifferentialDecayWidth(record) / self._total_width

    def DensityVariables(self):
        return ["cos_theta"]

    def Convention(self):
        return _PhaseSpaceConvention.RestFrameSolidAngle

    def Topology(self):
        return _Topology.Decay2Body

    def Measure(self):
        return _Measure.SolidAngleRest

    def SecondaryMasses(self, secondary_types):
        return [self.m_lepton, self.m_nu]

    def SecondaryHelicities(self, record):
        return [0, 0]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other

    def SampleFinalState(self, record, random):
        P_parent = np.array(record.primary_momentum)
        p_cm = _two_body_p_cm(self.m_meson, self.m_lepton, self.m_nu)

        cos_theta = random.Uniform(-1.0, 1.0)
        phi = random.Uniform(0.0, 2.0 * math.pi)

        P_lep = _boost_to_lab(P_parent, p_cm, cos_theta, phi, self.m_lepton)
        P_nu = _boost_to_lab(P_parent, p_cm, -cos_theta, phi + math.pi, self.m_nu)

        for sec in record.get_secondary_particle_records():
            if int(sec.type) == self.pdgid_lepton:
                sec.four_momentum = P_lep
                sec.mass = self.m_lepton
            elif int(sec.type) == self.pdgid_neutrino:
                sec.four_momentum = P_nu
                sec.mass = self.m_nu
        return


# ===================================================================
#  Shared helpers for rest-frame <-> lab-frame kinematics
# ===================================================================

def _rest_frame_E_phi(E_V_lab, cos_theta_V_lab, E_pi, p_pi, m_M, m_phi):
    """Compute rest-frame E_phi from lab-frame V1 kinematics via inverse boost.

    The pion defines +z.  The boost from lab to rest frame has
    beta = p_pi / E_pi,  gamma = E_pi / m_M.
    E_phi_rf = gamma * (E_V_lab - beta * p_V_lab * cos_theta_V_lab)
    """
    gamma = E_pi / m_M
    beta = p_pi / E_pi if E_pi > 0 else 0.0
    p_V_lab = math.sqrt(max(E_V_lab**2 - m_phi**2, 0.0))
    return gamma * (E_V_lab - beta * p_V_lab * cos_theta_V_lab)


def _dGamma_dEV_dcosV(E_V_lab, cos_theta_V_lab, E_pi, p_pi, decay_obj):
    """Marginal differential width d^2 Gamma / (dE_V d cos_theta_V) in lab frame.

    Integrates the rest-frame matrix element over the neutrino energy E_nu
    at fixed rest-frame E_phi (determined by the lab-frame V1 kinematics
    via inverse boost).

    Returns the differential width in GeV^-1 (per unit E_V per unit cos_theta_V).
    """
    d = decay_obj
    m_M = d.m_M
    m_l = d.m_l
    m_phi = d.m_phi

    E_phi_rf = _rest_frame_E_phi(E_V_lab, cos_theta_V_lab, E_pi, p_pi, m_M, m_phi)

    if E_phi_rf < m_phi or E_phi_rf > d.E_phi_max:
        return 0.0

    u_val = m_M**2 + m_phi**2 - 2.0 * m_M * E_phi_rf
    if u_val <= m_l**2:
        return 0.0

    mq = math.sqrt(u_val)
    E_mu_star = (u_val + m_l**2) / (2.0 * mq)
    p_mu_star = math.sqrt(max(E_mu_star**2 - m_l**2, 0.0))

    E_q = m_M - E_phi_rf
    p_q = math.sqrt(max(E_q**2 - u_val, 0.0))
    gamma_sub = E_q / mq
    beta_gamma_sub = p_q / mq

    E_mu_max = gamma_sub * E_mu_star + beta_gamma_sub * p_mu_star
    E_mu_min = max(gamma_sub * E_mu_star - beta_gamma_sub * p_mu_star, m_l)

    if E_mu_max <= E_mu_min:
        return 0.0

    def integrand(Emu):
        Enu = m_M - E_phi_rf - Emu
        if Enu < 0.0:
            return 0.0
        return d._matel_sq(Enu, E_phi_rf)

    try:
        val, _ = _integrate.quad(integrand, E_mu_min, E_mu_max,
                                 limit=40, epsrel=1e-3)
    except Exception:
        val = 0.0

    prefactor_rf = 1.0 / (64.0 * math.pi**3 * m_M)

    gamma = E_pi / m_M
    beta = p_pi / E_pi if E_pi > 0 else 0.0
    p_V_lab = math.sqrt(max(E_V_lab**2 - m_phi**2, 0.0))
    jacobian = gamma * abs(1.0 - beta * E_V_lab / max(p_V_lab, 1e-30) * cos_theta_V_lab) if p_V_lab > 0 else 1.0

    return max(prefactor_rf * val / max(jacobian, 1e-30), 0.0)


def _sample_rest_frame(decay_obj, max_matel, random):
    """Rejection-sample (E_nu, E_phi) from the matrix element in the rest frame."""
    d = decay_obj
    for _ in range(10000):
        E_nu = random.Uniform(0.0, d.E_nu_max)
        lims = d._E_phi_limits(E_nu)
        if lims[0] is None:
            continue
        E_phi = random.Uniform(lims[0], lims[1])
        if d.m_M - E_nu - E_phi < d.m_l:
            continue
        matel = d._matel_sq(E_nu, E_phi)
        u = random.Uniform(0.0, max_matel)
        if u <= matel:
            return E_nu, E_phi
    lims = d._E_phi_limits(d.E_nu_max * 0.5)
    return d.E_nu_max * 0.5, 0.5 * (lims[0] + lims[1]) if lims[0] is not None else d.m_phi


def _build_lab_momenta(E_nu_rf, E_phi_rf, decay_obj, P_parent, random):
    """Given rest-frame (E_nu, E_phi), construct lab-frame 4-momenta for all 3 daughters."""
    d = decay_obj
    m_M, m_l, m_phi = d.m_M, d.m_l, d.m_phi
    E_l = m_M - E_nu_rf - E_phi_rf

    p_nu = E_nu_rf
    p_l = math.sqrt(max(E_l**2 - m_l**2, 0.0))
    p_phi = math.sqrt(max(E_phi_rf**2 - m_phi**2, 0.0))

    cos_nu = random.Uniform(-1.0, 1.0)
    phi_nu = random.Uniform(0.0, 2.0 * math.pi)

    P_nu_lab = _boost_to_lab(P_parent, p_nu, cos_nu, phi_nu, 0.0)

    if p_nu > 0 and p_phi > 0:
        cos_phi_rf = ((m_M - E_nu_rf)**2 - E_phi_rf**2 - E_l**2 + m_l**2
                      + 2.0 * E_nu_rf * E_phi_rf) / (2.0 * p_nu * p_phi)
        cos_phi_rf = max(-1.0, min(1.0, cos_phi_rf))
    else:
        cos_phi_rf = 0.0

    sin_nu = math.sqrt(max(1.0 - cos_nu**2, 0.0))
    nu_dir = np.array([sin_nu * math.cos(phi_nu), sin_nu * math.sin(phi_nu), cos_nu])

    perp1 = np.cross(nu_dir, np.array([0, 0, 1]))
    if np.linalg.norm(perp1) < 1e-10:
        perp1 = np.cross(nu_dir, np.array([0, 1, 0]))
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(nu_dir, perp1)

    azimuth = random.Uniform(0.0, 2.0 * math.pi)
    sin_phi_rf = math.sqrt(max(1.0 - cos_phi_rf**2, 0.0))
    phi_dir = (cos_phi_rf * nu_dir
               + sin_phi_rf * math.cos(azimuth) * perp1
               + sin_phi_rf * math.sin(azimuth) * perp2)

    P_phi_lab = _boost_to_lab(
        P_parent, p_phi,
        float(phi_dir[2] / max(np.linalg.norm(phi_dir), 1e-12)),
        math.atan2(float(phi_dir[1]), float(phi_dir[0])),
        m_phi,
    )

    P_l_rf_3 = np.array([m_M, 0, 0, 0]) - np.array([E_nu_rf, *(p_nu * nu_dir)]) \
               - np.array([E_phi_rf, *(p_phi * phi_dir)])
    p_l_3 = P_l_rf_3[1:]
    p_l_mag = np.linalg.norm(p_l_3)
    if p_l_mag > 0:
        l_dir = p_l_3 / p_l_mag
        P_l_lab = _boost_to_lab(P_parent, p_l, float(l_dir[2]),
                                math.atan2(float(l_dir[1]), float(l_dir[0])), m_l)
    else:
        P_l_lab = _boost_to_lab(P_parent, p_l, 0.0, 0.0, m_l)

    return P_nu_lab, P_l_lab, P_phi_lab


def _extract_V_lab_angles(P_parent, P_V_lab, m_M):
    """Extract (E_V, cos_theta_V, phi_V) in the frame where +z = pion direction."""
    p_pi_3 = np.array(P_parent[1:])
    p_pi_mag = np.linalg.norm(p_pi_3)
    if p_pi_mag < 1e-12:
        return P_V_lab[0], 0.0, 0.0

    z_hat = p_pi_3 / p_pi_mag
    p_V_3 = np.array(P_V_lab[1:])
    p_V_mag = np.linalg.norm(p_V_3)
    if p_V_mag < 1e-12:
        return P_V_lab[0], 0.0, 0.0

    cos_theta_V = np.dot(p_V_3, z_hat) / p_V_mag

    arb = np.array([0, 1, 0]) if abs(z_hat[1]) < 0.9 else np.array([1, 0, 0])
    x_hat = np.cross(z_hat, arb)
    x_hat /= np.linalg.norm(x_hat)
    y_hat = np.cross(z_hat, x_hat)

    phi_V = math.atan2(np.dot(p_V_3, y_hat), np.dot(p_V_3, x_hat))
    return P_V_lab[0], float(cos_theta_V), phi_V


# ===================================================================
#  MesonThreeBodySIRENDecay -- pi/K -> l nu V1 (physical, SIREN interface)
# ===================================================================

class MesonThreeBodySIRENDecay(_Decay):
    """
    Physical three-body meson decay M -> l nu V1 for SIREN injection.

    DifferentialDecayWidth and FinalStateProbability are both
    differential in lab-frame V1 variables (E_V, cos_theta_V)
    with +z along the pion momentum direction.

    Since this is the physical (unbiased) version,
    FinalStateProbability = DifferentialDecayWidth / TotalDecayWidth.
    """

    def __init__(
        self,
        m_meson=_M_PI,
        m_lepton=_M_MU,
        m_mediator=0.017,
        g_mu=1.0,
        mediator_type="scalar",
        *,
        pdgid_meson=_PDGID_PIPLUS,
        pdgid_lepton=_PDGID_MUPLUS,
        pdgid_neutrino=_PDGID_NUMU,
        pdgid_mediator=5922,
        table_dir=None,
    ):
        _Decay.__init__(self)

        self.m_meson = m_meson
        self.m_lepton = m_lepton
        self.m_mediator = m_mediator
        self.m_nu = 0.0

        self.pdgid_meson = pdgid_meson
        self.pdgid_lepton = pdgid_lepton
        self.pdgid_neutrino = pdgid_neutrino
        self.pdgid_mediator = pdgid_mediator

        self._decay = MesonThreeBodyDecay(
            m_meson, m_lepton, m_mediator, g_mu, mediator_type
        )

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)

        self._total_width = self._decay.total_width()
        self._max_matel = self._find_max_matel()

    def _find_max_matel(self, n_samples=10000):
        d = self._decay
        max_val = 0.0
        for _ in range(n_samples):
            E_nu = np.random.uniform(0, d.E_nu_max)
            lims = d._E_phi_limits(E_nu)
            if lims[0] is None:
                continue
            E_phi = np.random.uniform(lims[0], lims[1])
            if d.m_M - E_nu - E_phi < d.m_l:
                continue
            val = d._matel_sq(E_nu, E_phi)
            if val > max_val:
                max_val = val
        return max_val * 1.2

    # -- SIREN Decay interface --

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_meson)
        sig.target_type = Particle.ParticleType.Decay
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_lepton),
            Particle.ParticleType(self.pdgid_neutrino),
            Particle.ParticleType(self.pdgid_mediator),
        ]
        return [sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if int(primary_type) == self.pdgid_meson:
            return self.GetPossibleSignatures()
        return []

    def TotalDecayWidthAllFinalStates(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def TotalDecayWidth(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def DifferentialDecayWidth(self, record):
        """Physical d^2Gamma / (dE_V d cos_theta_V) evaluated at the
        V1 kinematics in this record, marginalised over nu/lepton DOF."""
        if int(record.signature.primary_type) != self.pdgid_meson:
            return 0.0

        P_parent = np.array(record.primary_momentum)
        E_pi = P_parent[0]
        p_pi = math.sqrt(max(E_pi**2 - self.m_meson**2, 0.0))

        for idx, stype in enumerate(record.signature.secondary_types):
            if int(stype) == self.pdgid_mediator:
                P_V = np.array(record.secondary_momenta[idx])
                break
        else:
            return 0.0

        E_V, cos_V, _ = _extract_V_lab_angles(P_parent, P_V, self.m_meson)
        return _dGamma_dEV_dcosV(E_V, cos_V, E_pi, p_pi, self._decay)

    def FinalStateProbability(self, record):
        """For the physical decay, the generation density equals the
        physical density: p_gen = (dGamma/dE_V dcosV) / Gamma_total."""
        if self._total_width <= 0:
            return 0.0
        return self.DifferentialDecayWidth(record) / self._total_width

    def DensityVariables(self):
        return ["E_V", "cos_theta_V"]

    def Convention(self):
        return _PhaseSpaceConvention.Custom

    def Topology(self):
        return _Topology.Decay3Body

    def Measure(self):
        return _Measure.Unspecified

    def SecondaryMasses(self, secondary_types):
        return [self.m_lepton, self.m_nu, self.m_mediator]

    def SecondaryHelicities(self, record):
        return [0, 0, 0]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other

    def SampleFinalState(self, record, random):
        """Sample from the physical matrix element via rejection sampling."""
        E_nu_rf, E_phi_rf = _sample_rest_frame(self._decay, self._max_matel, random)
        P_parent = np.array(record.primary_momentum)
        P_nu, P_l, P_V = _build_lab_momenta(E_nu_rf, E_phi_rf, self._decay, P_parent, random)

        for sec in record.get_secondary_particle_records():
            pid = int(sec.type)
            if pid == self.pdgid_lepton:
                sec.four_momentum = P_l
                sec.mass = self._decay.m_l
            elif pid == self.pdgid_neutrino:
                sec.four_momentum = P_nu
                sec.mass = 0.0
            elif pid == self.pdgid_mediator:
                sec.four_momentum = P_V
                sec.mass = self._decay.m_phi
        return


# ===================================================================
#  BiasedMesonThreeBodyDecay -- pi/K -> l nu V1 (biased, cone-directed V1)
# ===================================================================

class BiasedMesonThreeBodyDecay(_Decay):
    """
    Biased three-body meson decay M -> l nu V1 for SIREN injection.

    The V1 direction is sampled from an energy-dependent cone pointed
    at the detector center.  The cone opening angle adapts to the
    V1 lab energy: for a given V1 energy, the cone covers the maximum
    angular spread of chi particles from the subsequent V1 -> chi chi
    decay (determined by m_chi and the V1 boost).

    This ensures that no chi-reachable phase space is cut regardless
    of the V1 energy, while maintaining high injection efficiency for
    boosted V1.

    DifferentialDecayWidth returns the same physical rate as the
    unbiased version.  FinalStateProbability returns the biased
    generation density in the same lab-frame variables.
    """

    def __init__(
        self,
        m_meson=_M_PI,
        m_lepton=_M_MU,
        m_mediator=0.017,
        m_chi=0.008,
        g_mu=1.0,
        mediator_type="scalar",
        *,
        detector_position=(0.0, 0.0, 0.0),
        cone_half_angle=None,
        detector_radius=2.0,
        safety_factor=1.5,
        pdgid_meson=_PDGID_PIPLUS,
        pdgid_lepton=_PDGID_MUPLUS,
        pdgid_neutrino=_PDGID_NUMU,
        pdgid_mediator=5922,
        table_dir=None,
    ):
        _Decay.__init__(self)

        self.m_meson = m_meson
        self.m_lepton = m_lepton
        self.m_mediator = m_mediator
        self.m_nu = 0.0

        self.pdgid_meson = pdgid_meson
        self.pdgid_lepton = pdgid_lepton
        self.pdgid_neutrino = pdgid_neutrino
        self.pdgid_mediator = pdgid_mediator

        self.detector_position = np.array(detector_position, dtype=float)
        self.m_chi = m_chi
        self.cone_half_angle = cone_half_angle
        self.detector_radius = detector_radius
        self.safety_factor = safety_factor

        self._decay = MesonThreeBodyDecay(
            m_meson, m_lepton, m_mediator, g_mu, mediator_type
        )

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)

        self._total_width = self._decay.total_width()
        self._max_matel = self._find_max_matel()

        # Precompute chi kinematics for energy-dependent cone
        self._p_cm_chi = _two_body_p_cm(m_mediator, m_chi, m_chi)
        self._E_chi_rf = math.sqrt(self._p_cm_chi**2 + m_chi**2) if m_mediator >= 2 * m_chi else 0.0
        self._beta_chi = self._p_cm_chi / self._E_chi_rf if self._E_chi_rf > 0 else 0.0

    def _find_max_matel(self, n_samples=10000):
        d = self._decay
        max_val = 0.0
        for _ in range(n_samples):
            E_nu = np.random.uniform(0, d.E_nu_max)
            lims = d._E_phi_limits(E_nu)
            if lims[0] is None:
                continue
            E_phi = np.random.uniform(lims[0], lims[1])
            if d.m_M - E_nu - E_phi < d.m_l:
                continue
            val = d._matel_sq(E_nu, E_phi)
            if val > max_val:
                max_val = val
        return max_val * 1.2

    def _chi_theta_max(self, E_V1):
        """Maximum lab-frame chi angle from V1 -> chi chi at given V1 energy."""
        m_V = self.m_mediator
        if E_V1 <= m_V or self._E_chi_rf <= 0:
            return math.pi
        gamma = E_V1 / m_V
        beta_V = math.sqrt(1.0 - 1.0 / gamma**2)
        if self._beta_chi >= beta_V:
            return math.pi
        sin_max = self._p_cm_chi / (gamma * self._E_chi_rf * beta_V)
        return math.asin(min(sin_max, 1.0))

    def _get_cone_params(self, decay_vertex, E_V1_lab=None):
        """Compute cone axis and energy-dependent half-angle."""
        delta = self.detector_position - np.array(decay_vertex)
        dist = np.linalg.norm(delta)
        if dist < 1e-6:
            return np.array([0, 0, 1.0]), math.pi
        axis = delta / dist

        if self.cone_half_angle is not None:
            half_angle = self.cone_half_angle
        else:
            theta_det = math.atan2(self.detector_radius, dist)
            if E_V1_lab is not None and self._E_chi_rf > 0:
                theta_chi = self._chi_theta_max(E_V1_lab)
                half_angle = max(theta_det, theta_chi * self.safety_factor)
            else:
                half_angle = theta_det
            half_angle = min(half_angle, math.pi)

        return axis, half_angle

    def _cone_solid_angle(self, half_angle):
        """Solid angle of a cone with given half-angle."""
        return 2.0 * math.pi * (1.0 - math.cos(half_angle))

    # -- SIREN Decay interface --

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_meson)
        sig.target_type = Particle.ParticleType.Decay
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_lepton),
            Particle.ParticleType(self.pdgid_neutrino),
            Particle.ParticleType(self.pdgid_mediator),
        ]
        return [sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if int(primary_type) == self.pdgid_meson:
            return self.GetPossibleSignatures()
        return []

    def TotalDecayWidthAllFinalStates(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def TotalDecayWidth(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_meson:
            return 0.0
        return self._total_width

    def DifferentialDecayWidth(self, record):
        """Physical d^2Gamma / (dE_V d cos_theta_V) -- same as unbiased."""
        if int(record.signature.primary_type) != self.pdgid_meson:
            return 0.0

        P_parent = np.array(record.primary_momentum)
        E_pi = P_parent[0]
        p_pi = math.sqrt(max(E_pi**2 - self.m_meson**2, 0.0))

        for idx, stype in enumerate(record.signature.secondary_types):
            if int(stype) == self.pdgid_mediator:
                P_V = np.array(record.secondary_momenta[idx])
                break
        else:
            return 0.0

        E_V, cos_V, _ = _extract_V_lab_angles(P_parent, P_V, self.m_meson)
        return _dGamma_dEV_dcosV(E_V, cos_V, E_pi, p_pi, self._decay)

    def FinalStateProbability(self, record):
        """Biased generation density.

        Direction: uniform within the cone = 1 / Omega_cone.
        Energy: proportional to the physical marginal dGamma/dE_V
                (we rejection-sample E_phi_rf from the physical distribution).

        The combined density in (E_V, cos_theta_V) is:
            p_gen = p(E_V) * (1 / Omega_cone)
        where p(E_V) = (1/Gamma_total) * integral_cosV dGamma/(dE_V dcosV) dcosV
        evaluated by integrating over all cos_theta_V (the physical marginal in E_V).
        """
        if self._total_width <= 0:
            return 0.0

        P_parent = np.array(record.primary_momentum)
        E_pi = P_parent[0]
        p_pi = math.sqrt(max(E_pi**2 - self.m_meson**2, 0.0))

        for idx, stype in enumerate(record.signature.secondary_types):
            if int(stype) == self.pdgid_mediator:
                P_V = np.array(record.secondary_momenta[idx])
                break
        else:
            return 0.0

        E_V, cos_V, _ = _extract_V_lab_angles(P_parent, P_V, self.m_meson)

        _, half_angle = self._get_cone_params(record.interaction_vertex, E_V)
        omega = self._cone_solid_angle(half_angle)

        E_phi_rf = _rest_frame_E_phi(E_V, cos_V, E_pi, p_pi, self.m_meson, self.m_mediator)

        d = self._decay
        if E_phi_rf < d.m_phi or E_phi_rf > d.E_phi_max:
            return 0.0

        dGamma_dEphi = d.differential_decay_rate(np.array([E_phi_rf]))[0]
        p_E = dGamma_dEphi / self._total_width if self._total_width > 0 else 0.0

        gamma = E_pi / self.m_meson
        beta = p_pi / E_pi if E_pi > 0 else 0.0
        p_V_lab = math.sqrt(max(E_V**2 - self.m_mediator**2, 0.0))
        jacobian = gamma * abs(1.0 - beta * E_V / max(p_V_lab, 1e-30) * cos_V) if p_V_lab > 0 else 1.0

        p_EV = p_E / max(jacobian, 1e-30)

        return p_EV / omega

    def DensityVariables(self):
        return ["lab_E_V", "lab_cos_theta_V"]

    def Convention(self):
        return _PhaseSpaceConvention.Custom

    def Topology(self):
        return _Topology.Decay3Body

    def Measure(self):
        return _Measure.Unspecified

    def SecondaryMasses(self, secondary_types):
        return [self.m_lepton, self.m_nu, self.m_mediator]

    def SecondaryHelicities(self, record):
        return [0, 0, 0]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other

    def SampleFinalState(self, record, random):
        """Sample V1 direction from cone toward detector, energy from physical distribution."""
        P_parent = np.array(record.primary_momentum)
        E_pi = P_parent[0]
        p_pi = math.sqrt(max(E_pi**2 - self.m_meson**2, 0.0))
        gamma = E_pi / self.m_meson
        beta = p_pi / E_pi if E_pi > 0 else 0.0

        decay_vertex = np.array(record.interaction_vertex)

        E_nu_rf, E_phi_rf = _sample_rest_frame(self._decay, self._max_matel, random)

        p_phi_rf = math.sqrt(max(E_phi_rf**2 - self.m_mediator**2, 0.0))
        E_V_lab_fwd = gamma * (E_phi_rf + beta * p_phi_rf)
        E_V_lab_bwd = gamma * (E_phi_rf - beta * p_phi_rf)
        E_V_lab_avg = (E_V_lab_fwd + E_V_lab_bwd) / 2
        p_V_lab = math.sqrt(max(E_V_lab_avg**2 - self.m_mediator**2, 0.0))

        # Energy-dependent cone: adapts to the chi opening angle at this V1 energy
        cone_axis, half_angle = self._get_cone_params(decay_vertex, E_V_lab_avg)

        cos_cone = random.Uniform(math.cos(half_angle), 1.0)
        phi_cone = random.Uniform(0.0, 2.0 * math.pi)
        sin_cone = math.sqrt(max(1.0 - cos_cone**2, 0.0))

        arb = np.array([0, 1, 0]) if abs(cone_axis[1]) < 0.9 else np.array([1, 0, 0])
        perp1 = np.cross(cone_axis, arb)
        perp1 /= np.linalg.norm(perp1)
        perp2 = np.cross(cone_axis, perp1)

        V_dir = (cos_cone * cone_axis
                 + sin_cone * math.cos(phi_cone) * perp1
                 + sin_cone * math.sin(phi_cone) * perp2)

        cos_theta_V_pion = np.dot(V_dir, P_parent[1:] / max(p_pi, 1e-12))
        E_V_lab = gamma * (E_phi_rf + beta * p_phi_rf * cos_theta_V_pion)
        p_V_lab = math.sqrt(max(E_V_lab**2 - self.m_mediator**2, 0.0))

        P_V_lab = np.array([E_V_lab, *(p_V_lab * V_dir)])

        P_parent_4 = np.array([E_pi, 0, 0, p_pi])
        P_V_rf = np.array([E_phi_rf, *(p_phi_rf * np.array([0, 0, 1]))])
        P_nu_rf = np.array([E_nu_rf, 0, 0, -E_nu_rf])

        E_l = self.m_meson - E_nu_rf - E_phi_rf
        p_l = math.sqrt(max(E_l**2 - self.m_lepton**2, 0.0))
        P_l_rf_3 = -P_nu_rf[1:] - P_V_rf[1:]
        p_l_mag = np.linalg.norm(P_l_rf_3)
        if p_l_mag > 0:
            P_l_rf_3 = P_l_rf_3 / p_l_mag * p_l

        P_nu_lab = _boost_to_lab(P_parent, E_nu_rf, 0.0, 0.0, 0.0)
        P_l_lab = _boost_to_lab(P_parent, p_l,
                                float(P_l_rf_3[2] / max(np.linalg.norm(P_l_rf_3), 1e-12)),
                                math.atan2(float(P_l_rf_3[1]), float(P_l_rf_3[0])),
                                self.m_lepton)

        for sec in record.get_secondary_particle_records():
            pid = int(sec.type)
            if pid == self.pdgid_lepton:
                sec.four_momentum = P_l_lab
                sec.mass = self.m_lepton
            elif pid == self.pdgid_neutrino:
                sec.four_momentum = P_nu_lab
                sec.mass = 0.0
            elif pid == self.pdgid_mediator:
                sec.four_momentum = P_V_lab
                sec.mass = self.m_mediator
        return


# ===================================================================
#  build_phi_flux  --  construct mediator flux at detector
# ===================================================================

def build_phi_flux(
    m_meson,
    m_lepton,
    m_phi,
    g_mu,
    mediator_type="scalar",
    flux_tag="pion_numu",
    min_energy=0.0,
    max_energy=3.0,
    n_bins=50,
    physically_normalized=True,
):
    """
    Construct the scalar/pseudoscalar mediator phi flux at the detector
    by convolving the parent meson spectrum with the three-body
    differential decay rate and boosting to the lab frame.

    Returns a siren.distributions.TabulatedFluxDistribution.
    """
    import siren

    raw_flux = siren.utilities.load_flux(
        "PionKaon",
        tag=flux_tag,
        physically_normalized=physically_normalized,
    )

    try:
        decay = MesonThreeBodyDecay(m_meson, m_lepton, m_phi, g_mu, mediator_type)
    except ValueError:
        energies = [min_energy, max_energy]
        flux_arr = [0.0, 0.0]
        return siren.distributions.TabulatedFluxDistribution(
            min_energy, max_energy, energies, flux_arr, physically_normalized
        )

    E_phi_rf_max = decay.E_phi_max
    if E_phi_rf_max <= m_phi:
        energies = [min_energy, max_energy]
        flux_arr = [0.0, 0.0]
        return siren.distributions.TabulatedFluxDistribution(
            min_energy, max_energy, energies, flux_arr, physically_normalized
        )

    n_rf = 200
    E_phi_rf = np.linspace(m_phi, E_phi_rf_max, n_rf)
    dGamma = decay.differential_decay_rate(E_phi_rf)

    E_phi_out = np.linspace(min_energy, max_energy, n_bins)
    phi_flux = np.zeros(n_bins)
    dEout = (max_energy - min_energy) / (n_bins - 1) if n_bins > 1 else 1.0
    dErf = (E_phi_rf_max - m_phi) / (n_rf - 1) if n_rf > 1 else 1.0

    nu_energies = list(raw_flux.GetEnergyNodes())

    E_nu_rf_2body = (m_meson**2 - m_lepton**2) / (2.0 * m_meson)
    nu_to_meson = m_meson / E_nu_rf_2body if E_nu_rf_2body > 0 else 1.0

    for E_nu in nu_energies:
        E_meson = E_nu * nu_to_meson
        if E_meson < m_meson:
            continue

        meson_flux = raw_flux.EvaluatePDF(E_nu)
        if meson_flux <= 0.0:
            continue

        p_meson = math.sqrt(max(E_meson**2 - m_meson**2, 0.0))
        gamma = E_meson / m_meson
        beta = p_meson / E_meson if E_meson > 0 else 0.0

        for j in range(n_rf):
            Erf = E_phi_rf[j]
            dG = dGamma[j]
            if dG <= 0.0:
                continue

            p_rf = math.sqrt(max(Erf**2 - m_phi**2, 0.0))
            E_lo = gamma * (Erf - beta * p_rf)
            E_hi = gamma * (Erf + beta * p_rf)

            if E_hi <= min_energy or E_lo >= max_energy:
                continue
            dE_lab = max(E_hi - E_lo, 1e-9)
            density = dG * dErf / dE_lab * meson_flux

            for k in range(n_bins):
                Eout = E_phi_out[k]
                bin_lo = Eout - 0.5 * dEout
                bin_hi = Eout + 0.5 * dEout
                if bin_lo < E_hi and bin_hi > E_lo:
                    overlap = min(bin_hi, E_hi) - max(bin_lo, E_lo)
                    phi_flux[k] += density * overlap

    return siren.distributions.TabulatedFluxDistribution(
        min_energy, max_energy, list(E_phi_out), list(phi_flux), physically_normalized
    )
