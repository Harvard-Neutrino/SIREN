"""
Vector-portal dark matter processes for SIREN.

Implements the physics from Dutta, Kim, Thompson, Thornton, Van de Water,
PRL 129, 111803 (2022) [arXiv:2110.11944]:

    Production:  pi/K -> l nu V1  (three-body meson decay)
    Decay:       V1 -> chi chi'   (dark photon to DM pair)
    Scattering:  chi N -> chi' N  (upscattering via t-channel V2)
    Decay:       chi' -> chi V1   (de-excitation)
    Decay:       V1 -> e+ e-      (visible signal)

All classes are self-contained with no DarkNews imports.
Cross-section classes implement the ups_case duck-type interface
expected by PyDarkNewsCrossSection.  Decay classes inherit directly
from siren.interactions.DarkNewsDecay.
"""

import os
import math
import numpy as np
import scipy.integrate as _integrate

from siren.interactions import Decay as _Decay, CrossSection as _CrossSection
from siren import DecayModel as _DecayModel, CrossSectionModel as _CrossSectionModel
from siren import dataclasses
from siren.dataclasses import Particle
from siren.injection import PhaseSpaceConvention as _PhaseSpaceConvention
from siren.injection import PhaseSpaceTopology as _Topology
from siren.injection import PhaseSpaceMeasure as _Measure
from siren.injection import BreitWignerMapping as _BreitWignerMapping

_ALPHA_EM = 1.0 / 137.036
_GEV2_TO_CM2 = 3.8938e-28
_M_ELECTRON = 0.000511  # GeV


# ---------------------------------------------------------------------------
# Lightweight particle / target stubs
# ---------------------------------------------------------------------------

class _Stub:
    def __init__(self, pdgid, mass=0.0, name=""):
        self.pdgid = pdgid
        self.mass = mass
        self.name = name

    def __eq__(self, other):
        if isinstance(other, int):
            return self.pdgid == other
        return self.pdgid == getattr(other, "pdgid", None)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __int__(self):
        return self.pdgid

    def __hash__(self):
        return hash(self.pdgid)

    def __repr__(self):
        return f"_Stub({self.pdgid}, {self.name!r})"


# ---------------------------------------------------------------------------
# Kinematic helpers (inlined from DarkNews phase_space)
# ---------------------------------------------------------------------------

def _Q2max(E, m_ups, M):
    """Maximum Q2 for 2->2 scattering m1 + M -> m3 + M at lab energy E.
    Q2 = -(p1 - p3)^2 = 2*p1cm*p3cm*(1 + cos_theta_cm) at backward scattering."""
    m1 = 0.0  # approximation: m_chi << everything else
    s = m1**2 + M**2 + 2.0 * M * E
    if s <= (m_ups + M)**2:
        return 0.0
    E1cm = (s + m1**2 - M**2) / (2.0 * math.sqrt(s))
    E3cm = (s + m_ups**2 - M**2) / (2.0 * math.sqrt(s))
    p1cm = math.sqrt(max(E1cm**2 - m1**2, 0.0))
    p3cm = math.sqrt(max(E3cm**2 - m_ups**2, 0.0))
    t_min = m1**2 + m_ups**2 - 2.0 * E1cm * E3cm - 2.0 * p1cm * p3cm
    return -t_min


def _Q2min(E, m_ups, M):
    """Minimum Q2 (forward scattering)."""
    m1 = 0.0
    s = m1**2 + M**2 + 2.0 * M * E
    if s <= (m_ups + M)**2:
        return 0.0
    E1cm = (s + m1**2 - M**2) / (2.0 * math.sqrt(s))
    E3cm = (s + m_ups**2 - M**2) / (2.0 * math.sqrt(s))
    p1cm = math.sqrt(max(E1cm**2 - m1**2, 0.0))
    p3cm = math.sqrt(max(E3cm**2 - m_ups**2, 0.0))
    t_max = m1**2 + m_ups**2 - 2.0 * E1cm * E3cm + 2.0 * p1cm * p3cm
    return -t_max


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


def _primary_type(arg):
    if isinstance(arg, dataclasses.InteractionRecord):
        return arg.signature.primary_type
    return arg


def _cos_theta_in_parent_rest(P_parent, P_child):
    P_parent = np.asarray(P_parent, dtype=float)
    P_child = np.asarray(P_child, dtype=float)
    parent_p = P_parent[1:]
    parent_p_mag = np.linalg.norm(parent_p)
    child_p = P_child[1:]

    if parent_p_mag < 1e-12:
        child_p_mag = np.linalg.norm(child_p)
        if child_p_mag < 1e-30:
            return 1.0
        return float(child_p[2] / child_p_mag)

    parent_mass_sq = max(P_parent[0]**2 - parent_p_mag**2, 0.0)
    parent_mass = math.sqrt(parent_mass_sq)
    if parent_mass < 1e-12 or P_parent[0] <= 0.0:
        return 1.0

    beta = parent_p_mag / P_parent[0]
    gamma = P_parent[0] / parent_mass
    beta_hat = parent_p / parent_p_mag

    child_p_parallel = float(np.dot(child_p, beta_hat))
    child_p_perp = child_p - child_p_parallel * beta_hat
    child_p_parallel_rest = gamma * (child_p_parallel - beta * P_child[0])
    child_p_rest = child_p_perp + child_p_parallel_rest * beta_hat
    child_p_rest_mag = np.linalg.norm(child_p_rest)
    if child_p_rest_mag < 1e-30:
        return 1.0
    return float(np.dot(child_p_rest, beta_hat) / child_p_rest_mag)


# ---------------------------------------------------------------------------
# Helm nuclear form factor
# ---------------------------------------------------------------------------

def _helm_F2(Q2, A):
    Q = math.sqrt(max(Q2, 0.0))
    Qfm = Q / 0.197327
    s = 0.9
    r0sq = max((1.2 * A**(1.0 / 3.0))**2 - 5.0 * s**2, 0.0)
    Qr = Qfm * math.sqrt(r0sq)
    if Qr < 1e-6:
        j1_over_Qr = 1.0 / 3.0
    else:
        j1_over_Qr = (math.sin(Qr) - Qr * math.cos(Qr)) / Qr**3
    return (3.0 * j1_over_Qr)**2 * math.exp(-(Qfm * s)**2)


# ===================================================================
#  VectorPortalUpsCase  --  chi N -> chi' N  upscattering
# ===================================================================

class VectorPortalUpsCase:
    """
    Dark-matter upscattering chi N -> chi' N via t-channel dark photon V2.

    Implements the duck-type interface expected by PyDarkNewsCrossSection:
        .nu_projectile, .nu_upscattered, .nuclear_target  (with .pdgid, .mass)
        .MA, .m_ups, .h_upscattered, .Ethreshold, .scattering_regime
        .total_xsec(E), .diff_xsec_Q2(E, Q2)
    """

    def __init__(
        self,
        m_chi,
        m_chi_prime,
        m_V,
        g_D,
        epsilon,
        *,
        pdgid_chi=5917,
        pdgid_chi_prime=5918,
        nuclear_pdgid=1000180400,
        nuclear_mass=37.215,
        nuclear_name="Ar40",
        A=40,
        Z=18,
        scattering_regime="coherent",
    ):
        self.m_chi = m_chi
        self.m_chi_prime = m_chi_prime
        self.m_V = m_V
        self.g_D = g_D
        self.epsilon = epsilon
        self.A = A
        self.Z = Z

        self.nu_projectile = _Stub(pdgid_chi, m_chi, "chi")
        self.nu_upscattered = _Stub(pdgid_chi_prime, m_chi_prime, "chi_prime")
        self.nuclear_target = _Stub(nuclear_pdgid, nuclear_mass, nuclear_name)
        self.nuclear_target.Z = Z
        self.nuclear_target.N = A - Z

        self.MA = nuclear_mass
        self.m_ups = m_chi_prime
        self.h_upscattered = 1
        self.scattering_regime = scattering_regime

        self.Ethreshold = ((m_chi_prime + nuclear_mass)**2
                          - m_chi**2 - nuclear_mass**2) / (2.0 * nuclear_mass)

    def _dsigma_dQ2(self, E, Q2):
        m1 = self.m_chi
        m3 = self.m_chi_prime
        M = self.MA
        mV = self.m_V

        s = m1**2 + M**2 + 2.0 * M * E
        flux_sq = (s - M**2)**2
        if flux_sq <= 0.0:
            return 0.0

        delta_m2 = m3**2 - m1**2
        numerator = 2.0 * M**2 * (2.0 * E * M - Q2 - delta_m2)
        if numerator <= 0.0:
            return 0.0

        propagator = 1.0 / (Q2 + mV**2)**2
        M2 = self.g_D**2 * 4.0 * math.pi * _ALPHA_EM * self.epsilon**2 * numerator * propagator
        F2 = _helm_F2(Q2, self.A)
        Q_eff_sq = (self.Z * self.epsilon)**2 if self.scattering_regime == "coherent" else 1.0

        dsig = M2 * F2 / (16.0 * math.pi * flux_sq)
        return max(0.0, dsig) * _GEV2_TO_CM2

    def diff_xsec_Q2(self, E, Q2):
        return np.array(self._dsigma_dQ2(E, Q2))

    def total_xsec(self, E):
        q2min = _Q2min(E, self.m_ups, self.MA)
        q2max = _Q2max(E, self.m_ups, self.MA)
        if q2max <= q2min:
            return 0.0
        result, _ = _integrate.quad(
            lambda q2: self._dsigma_dQ2(E, q2),
            q2min, q2max,
            limit=80, epsrel=1e-4,
        )
        return max(0.0, result)


# ===================================================================
#  VectorPortalUpscatteringXS  --  chi N -> chi' N
# ===================================================================

class VectorPortalUpscatteringXS(_CrossSectionModel):
    """On-shell vector-portal upscattering chi + N -> chi' + N.

    MandelstamQ2 Scatter2to2: the authoring base derives Topology, Measure,
    DensityVariables, FinalStateProbability (differential / total) and the
    interaction threshold from the total_xs()/differential_xs() hooks; sample()
    is overridden because MandelstamQ2 has no self-contained default channel.
    The HNucleus-aware target/signature methods and the flexible-argument
    TotalCrossSection/DifferentialCrossSection overloads are kept as-is.
    """

    measure = _Measure.MandelstamQ2()

    def __init__(
        self,
        m_chi,
        m_chi_prime,
        m_V2,
        g_D,
        epsilon,
        *,
        pdgid_chi=5917,
        pdgid_chi_prime=5918,
        nuclear_pdgid=1000180400,
        nuclear_mass=37.215,
        nuclear_name="Ar40",
        A=40,
        Z=18,
    ):
        _CrossSectionModel.__init__(self)

        self.m_chi = m_chi
        self.m_chi_prime = m_chi_prime
        self.m_V2 = m_V2
        self.g_D = g_D
        self.epsilon = epsilon
        self.A = A
        self.Z = Z

        self.pdgid_chi = pdgid_chi
        self.pdgid_chi_prime = pdgid_chi_prime
        self.nuclear_pdgid = nuclear_pdgid
        self.m_target = nuclear_mass

        self._ups = VectorPortalUpsCase(
            m_chi=m_chi,
            m_chi_prime=m_chi_prime,
            m_V=m_V2,
            g_D=g_D,
            epsilon=epsilon,
            pdgid_chi=pdgid_chi,
            pdgid_chi_prime=pdgid_chi_prime,
            nuclear_pdgid=nuclear_pdgid,
            nuclear_mass=nuclear_mass,
            nuclear_name=nuclear_name,
            A=A,
            Z=Z,
        )

        self.primary = Particle.ParticleType(pdgid_chi)
        self.target = self.GetPossibleTargets()[0]
        self.finals = tuple(self.GetPossibleSignatures()[0].secondary_types)
        self.threshold = self._ups.Ethreshold

    def GetPossiblePrimaries(self):
        return [Particle.ParticleType(self.pdgid_chi)]

    def GetPossibleTargets(self):
        target_type = Particle.ParticleType(self.nuclear_pdgid)
        if target_type == Particle.ParticleType.PPlus:
            target_type = Particle.ParticleType.HNucleus
        return [target_type]

    def GetPossibleTargetsFromPrimary(self, primary_type):
        if int(primary_type) == self.pdgid_chi:
            return self.GetPossibleTargets()
        return []

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_chi)
        target_type = Particle.ParticleType(self.nuclear_pdgid)
        if target_type == Particle.ParticleType.PPlus:
            target_type = Particle.ParticleType.HNucleus
        sig.target_type = target_type
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_chi_prime),
            target_type,
        ]
        return [sig]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        if int(primary_type) == self.pdgid_chi:
            expected_target = Particle.ParticleType(self.nuclear_pdgid)
            if expected_target == Particle.ParticleType.PPlus:
                expected_target = Particle.ParticleType.HNucleus
            if target_type == expected_target:
                return self.GetPossibleSignatures()
        return []

    def TotalCrossSection(self, arg1, energy=None, target=None):
        if isinstance(arg1, dataclasses.InteractionRecord):
            energy = arg1.primary_momentum[0]
            primary = arg1.signature.primary_type
        else:
            primary = arg1
        if int(primary) != self.pdgid_chi:
            return 0.0
        return self._ups.total_xsec(energy)

    def DifferentialCrossSection(self, arg1, target=None, energy=None, Q2=None):
        if isinstance(arg1, dataclasses.InteractionRecord):
            record = arg1
            primary = np.array(record.primary_momentum, dtype=float)
            chi_prime = np.array(record.secondary_momenta[0], dtype=float)
            m1sq = max(0.0, primary[0]**2 - float(np.dot(primary[1:], primary[1:])))
            m3sq = max(0.0, chi_prime[0]**2 - float(np.dot(chi_prime[1:], chi_prime[1:])))
            p1p3 = primary[0] * chi_prime[0] - float(np.dot(primary[1:], chi_prime[1:]))
            Q2 = max(0.0, -(m1sq + m3sq - 2.0 * p1p3))
            energy = record.primary_momentum[0]
        return float(np.real(self._ups.diff_xsec_Q2(energy, Q2)))

    def total_xs(self, record):
        return self.TotalCrossSection(record)

    def differential_xs(self, record):
        return self.DifferentialCrossSection(record)

    def density_variables(self):
        return ["Q2"]

    def Q2Min(self, interaction):
        return _Q2min(interaction.primary_momentum[0], self.m_chi_prime, self.m_target)

    def Q2Max(self, interaction):
        return _Q2max(interaction.primary_momentum[0], self.m_chi_prime, self.m_target)

    def TargetMass(self, target_type):
        return self.m_target

    def SecondaryMasses(self, secondary_types):
        return [self.m_chi_prime, self.m_target]

    def SecondaryHelicities(self, record):
        return [record.primary_helicity, record.target_helicity]

    def Convention(self):
        return _PhaseSpaceConvention.MandelstamST

    def equal(self, other):
        return self is other

    def _sample_Q2(self, E_chi, random):
        """Rejection-sample Q2 from dsigma/dQ2.

        dsigma/dQ2 falls monotonically with Q2 (the V2 propagator
        dominates), so the envelope maximum sits at q2min.  Sampling from
        this density -- instead of the previous uniform draw -- makes the
        cross section satisfy Sample == Density: an unbiased single-channel
        run then yields the correct propagator-peaked Q2 spectrum, and the
        physical phase-space channel contributes a weight ratio of 1,
        removing the dominant source of event-weight variance.
        """
        q2min = _Q2min(E_chi, self.m_chi_prime, self.m_target)
        q2max = _Q2max(E_chi, self.m_chi_prime, self.m_target)
        if q2max <= q2min:
            return None
        f_max = self._ups._dsigma_dQ2(E_chi, q2min) * 1.5
        if f_max <= 0.0:
            return random.Uniform(q2min, q2max)
        for _ in range(10000):
            cand = random.Uniform(q2min, q2max)
            if random.Uniform(0.0, f_max) <= self._ups._dsigma_dQ2(E_chi, cand):
                return cand
        return random.Uniform(q2min, q2max)

    def sample(self, record, random):
        E_chi = record.primary_momentum[0]
        M = self.m_target
        m_chi = self.m_chi
        m_chi_prime = self.m_chi_prime

        Q2 = self._sample_Q2(E_chi, random)
        if Q2 is None:
            return

        s = m_chi**2 + M**2 + 2.0 * M * E_chi
        sqrt_s = math.sqrt(s)
        E_chi_prime_cm = (s + m_chi_prime**2 - M**2) / (2.0 * sqrt_s)
        E_N_cm = (s - m_chi_prime**2 + M**2) / (2.0 * sqrt_s)
        p_out_cm = math.sqrt(max(E_chi_prime_cm**2 - m_chi_prime**2, 0.0))

        E_in_cm = (s + m_chi**2 - M**2) / (2.0 * sqrt_s)
        p_in_cm = math.sqrt(max(E_in_cm**2 - m_chi**2, 0.0))

        if p_in_cm > 0.0 and p_out_cm > 0.0:
            cos_cm = 1.0 - Q2 / (2.0 * p_in_cm * p_out_cm)
        else:
            cos_cm = 0.0
        cos_cm = max(-1.0, min(1.0, cos_cm))

        gamma_cm = (E_chi + M) / sqrt_s
        beta_cm = math.sqrt(max(E_chi**2 - m_chi**2, 0.0)) / (E_chi + M)

        phi_cm = random.Uniform(0.0, 2.0 * math.pi)
        sin_cm = math.sqrt(max(1.0 - cos_cm**2, 0.0))

        px_cp = p_out_cm * sin_cm * math.cos(phi_cm)
        py_cp = p_out_cm * sin_cm * math.sin(phi_cm)
        pz_cp = p_out_cm * cos_cm

        E_cp_lab = gamma_cm * (E_chi_prime_cm + beta_cm * pz_cp)
        pz_cp_lab = gamma_cm * (pz_cp + beta_cm * E_chi_prime_cm)
        P_chi_prime = np.array([E_cp_lab, px_cp, py_cp, pz_cp_lab])

        E_N_lab = gamma_cm * (E_N_cm - beta_cm * pz_cp)
        pz_N_lab = gamma_cm * (-pz_cp + beta_cm * E_N_cm)
        P_N = np.array([E_N_lab, -px_cp, -py_cp, pz_N_lab])

        p_in_dir = np.array(record.primary_momentum[1:])
        p_in_mag = np.linalg.norm(p_in_dir)
        if p_in_mag > 1e-12:
            z_hat = p_in_dir / p_in_mag
            arb = np.array([0, 1, 0]) if abs(z_hat[1]) < 0.9 else np.array([1, 0, 0])
            x_hat = np.cross(z_hat, arb)
            x_hat /= np.linalg.norm(x_hat)
            y_hat = np.cross(z_hat, x_hat)
            R = np.column_stack([x_hat, y_hat, z_hat])

            P_chi_prime[1:] = R @ P_chi_prime[1:]
            P_N[1:] = R @ P_N[1:]

        for sec in record.get_secondary_particle_records():
            if int(sec.type) == self.pdgid_chi_prime:
                sec.four_momentum = P_chi_prime
                sec.mass = m_chi_prime
            else:
                sec.four_momentum = P_N
                sec.mass = M
        return


# ===================================================================
#  VectorPortalOffShellXS  --  chi N -> chi V1 N  (off-shell chi')
# ===================================================================

class VectorPortalOffShellXS(_CrossSection):
    """
    Off-shell chi' scattering: chi + N -> chi + V1 + N.

    Uses the narrow-width factorization internally:
        dsigma(chi N -> chi V1 N) = dsigma(chi N -> chi' N) * BR(chi' -> chi V1)
    with BR = 1 for the benchmark parameters.

    The cross section is identical to VectorPortalUpsCase, but the
    signature produces 3 secondaries (chi, V1, recoil nucleus) instead
    of 2 (chi', nucleus).  SampleFinalState generates chi' internally,
    decays it to chi + V1 isotropically, and returns 3 secondary momenta.
    """

    def __init__(
        self,
        m_chi,
        m_chi_prime,
        m_V1,
        m_V2,
        g_D,
        epsilon,
        *,
        pdgid_chi=5917,
        pdgid_V1=5922,
        nuclear_pdgid=1000180400,
        nuclear_mass=37.215,
        nuclear_name="Ar40",
        A=40,
        Z=18,
    ):
        _CrossSection.__init__(self)

        self.m_chi = m_chi
        self.m_chi_prime = m_chi_prime
        self.m_V1 = m_V1
        self.m_V2 = m_V2
        self.g_D = g_D
        self.epsilon = epsilon
        self.A = A
        self.Z = Z

        self.pdgid_chi = pdgid_chi
        self.pdgid_V1 = pdgid_V1
        self.nuclear_pdgid = nuclear_pdgid

        self.m_ups = m_chi_prime
        self.m_target = nuclear_mass

        self._ups = VectorPortalUpsCase(
            m_chi=m_chi, m_chi_prime=m_chi_prime, m_V=m_V2,
            g_D=g_D, epsilon=epsilon,
            pdgid_chi=pdgid_chi, pdgid_chi_prime=5918,
            nuclear_pdgid=nuclear_pdgid, nuclear_mass=nuclear_mass,
            nuclear_name=nuclear_name, A=A, Z=Z,
        )

    def GetPossiblePrimaries(self):
        return [Particle.ParticleType(self.pdgid_chi)]

    def GetPossibleTargets(self):
        target_type = Particle.ParticleType(self.nuclear_pdgid)
        if target_type == Particle.ParticleType.PPlus:
            target_type = Particle.ParticleType.HNucleus
        return [target_type]

    def GetPossibleTargetsFromPrimary(self, primary_type):
        if int(primary_type) == self.pdgid_chi:
            return self.GetPossibleTargets()
        return []

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_chi)
        target_type = Particle.ParticleType(self.nuclear_pdgid)
        if target_type == Particle.ParticleType.PPlus:
            target_type = Particle.ParticleType.HNucleus
        sig.target_type = target_type
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_chi),
            Particle.ParticleType(self.pdgid_V1),
            target_type,
        ]
        return [sig]

    def GetPossibleSignaturesFromParents(self, primary_type, target_type):
        if int(primary_type) == self.pdgid_chi:
            expected_target = Particle.ParticleType(self.nuclear_pdgid)
            if expected_target == Particle.ParticleType.PPlus:
                expected_target = Particle.ParticleType.HNucleus
            if target_type == expected_target:
                return self.GetPossibleSignatures()
        return []

    def TotalCrossSection(self, arg1, energy=None, target=None):
        if isinstance(arg1, dataclasses.InteractionRecord):
            energy = arg1.primary_momentum[0]
            primary = arg1.signature.primary_type
        else:
            primary = arg1
        if int(primary) != self.pdgid_chi:
            return 0.0
        return self._ups.total_xsec(energy)

    def DifferentialCrossSection(self, arg1, target=None, energy=None, Q2=None):
        if isinstance(arg1, dataclasses.InteractionRecord):
            record = arg1
            # The momentum transfer to the nucleus is fixed by the (off-shell)
            # chi', which is the chi(0) + V1(1) pair -- not the lone chi.
            p_in = np.array(record.primary_momentum, dtype=float)
            p_pair = (np.array(record.secondary_momenta[0], dtype=float)
                      + np.array(record.secondary_momenta[1], dtype=float))
            q_mu = p_in - p_pair
            Q2 = max(0.0, -(q_mu[0]**2 - float(np.dot(q_mu[1:], q_mu[1:]))))
            energy = record.primary_momentum[0]
        return float(np.real(self._ups.diff_xsec_Q2(energy, Q2)))

    def InteractionThreshold(self, interaction):
        return self._ups.Ethreshold

    def Q2Min(self, interaction):
        return _Q2min(interaction.primary_momentum[0], self.m_chi_prime, self.m_target)

    def Q2Max(self, interaction):
        return _Q2max(interaction.primary_momentum[0], self.m_chi_prime, self.m_target)

    def TargetMass(self, target_type):
        return self.m_target

    def SecondaryMasses(self, secondary_types):
        return [self.m_chi, self.m_V1, self.m_target]

    def SecondaryHelicities(self, record):
        return [record.primary_helicity, 0, record.target_helicity]

    def FinalStateProbability(self, record):
        """Physical density in the Recursive2Body measure for the 2->3 process.

        Factorization: chi N -> N(spectator=2) + {chi(0) V1(1)}(pair), where
        the (chi, V1) pair is the (off-shell) chi'.  In the narrow-width
        approximation the density factorizes into three normalized pieces:

            f = BW(s_pair)                                       [ds_pair]
              * (1/sigma_T) dsigma/dQ2 |dQ2/dcos_pair| / (2 pi)  [dOmega_pair]
              * 1 / (4 pi)                                       [dOmega_sub]

        The chi' production angle carries the physical dsigma/dQ2 (peaked
        forward by the V2 propagator); the chi' -> chi V1 decay is isotropic
        in the chi' rest frame.  SampleFinalState draws from exactly this
        distribution, so the channel density and its sampler are consistent.
        """
        E_chi = record.primary_momentum[0]
        M = self.m_target

        total_xs = self.TotalCrossSection(record)
        if total_xs <= 0:
            return 0.0

        # Pair = chi(0) + V1(1) = the off-shell chi'; its invariant mass is
        # the chi' virtuality s_pair.
        P_pair = (np.array(record.secondary_momenta[0], dtype=float)
                  + np.array(record.secondary_momenta[1], dtype=float))
        s_pair = P_pair[0]**2 - float(np.dot(P_pair[1:], P_pair[1:]))
        if s_pair <= 0.0:
            return 0.0

        diff_xs = self.DifferentialCrossSection(record)
        if diff_xs <= 0:
            return 0.0

        # |dQ2/dcos_pair| in the CM frame, using the actual pair mass so the
        # Jacobian matches the kinematics produced by SampleFinalState.
        s = self.m_chi**2 + M**2 + 2.0 * M * E_chi
        sqrt_s = math.sqrt(s)
        E_in_cm = (s + self.m_chi**2 - M**2) / (2.0 * sqrt_s)
        E_pair_cm = (s + s_pair - M**2) / (2.0 * sqrt_s)
        p_in_cm = math.sqrt(max(E_in_cm**2 - self.m_chi**2, 0.0))
        p_pair_cm = math.sqrt(max(E_pair_cm**2 - s_pair, 0.0))

        dQ2_dcos = 2.0 * p_in_cm * p_pair_cm
        if dQ2_dcos <= 0:
            return 0.0

        # Breit-Wigner in s_pair, normalized over the kinematically allowed
        # [s_min, s_max].  SampleFinalState draws from this SAME mapping
        # object, so the sampled and reported s_pair densities are identical
        # by construction (Contract C1).  The former code normalized the BW
        # over the whole real line (a 1/pi factor) while the sampler used the
        # truncated range (1/(hi-lo)); that made this reported density
        # integrate to (hi-lo)/pi < 1 instead of 1 -- an energy-dependent
        # Sample != Density mismatch, masked today only by the physical
        # channel's small mixture weight.
        bw_map = self._s_pair_mapping(E_chi)
        if bw_map is None:
            return 0.0
        bw = bw_map.Density(s_pair)

        # density per ds_pair * dOmega_pair * dOmega_sub:
        # = (dsigma/dQ2 * |dQ2/dcos|) / (sigma_total * 2*pi * 4*pi) * BW(s_pair)
        return diff_xs * dQ2_dcos * bw / (total_xs * 8.0 * math.pi**2)

    def _chi_prime_width(self):
        """Total decay width of chi' -> chi + V1."""
        m_cp = self.m_chi_prime
        m_chi = self.m_chi
        m_V1 = self.m_V1
        g_D = self._ups.g_D
        p_star = _two_body_p_cm(m_cp, m_chi, m_V1)
        return g_D**2 * p_star**3 / (6.0 * math.pi * m_cp**2)

    def DensityVariables(self):
        return ["s_pair", "cos_theta_sub"]

    def Convention(self):
        return _PhaseSpaceConvention.Recursive2Body

    def Topology(self):
        return _Topology.Scatter2to3

    def Measure(self):
        return _Measure.Recursive2Body(spectator=2, pair_first=0, pair_second=1)

    def equal(self, other):
        return self is other

    def _sample_Q2(self, E_chi, random):
        """Rejection-sample Q2 from dsigma/dQ2.

        dsigma/dQ2 falls monotonically with Q2 (the V2 propagator dominates),
        so the envelope maximum sits at q2min.
        """
        q2min = _Q2min(E_chi, self.m_chi_prime, self.m_target)
        q2max = _Q2max(E_chi, self.m_chi_prime, self.m_target)
        if q2max <= q2min:
            return None
        f_max = self._ups._dsigma_dQ2(E_chi, q2min) * 1.5
        if f_max <= 0.0:
            return random.Uniform(q2min, q2max)
        for _ in range(10000):
            cand = random.Uniform(q2min, q2max)
            if random.Uniform(0.0, f_max) <= self._ups._dsigma_dQ2(E_chi, cand):
                return cand
        return random.Uniform(q2min, q2max)

    def _s_pair_mapping(self, E_chi):
        """Shared Breit-Wigner map for the chi' virtuality s_pair over the
        kinematically allowed range [s_min, s_max].  BOTH SampleFinalState
        (via Forward) and FinalStateProbability (via Density) route through
        this one object, so the sampled and reported s_pair densities are
        identical by construction (Contract C1).  Returns None when the
        allowed range is degenerate.
        """
        m_cp = self.m_chi_prime
        width = self._chi_prime_width()
        s_min = (self.m_chi + self.m_V1)**2
        s = self.m_chi**2 + self.m_target**2 + 2.0 * self.m_target * E_chi
        s_max = (math.sqrt(s) - self.m_target)**2
        if s_max <= s_min or width <= 0.0:
            return None
        return _BreitWignerMapping(m_cp, width, s_min, s_max)

    def _sample_s_pair(self, E_chi, random):
        """Draw s_pair (the chi' virtuality) from the shared Breit-Wigner map."""
        bw_map = self._s_pair_mapping(E_chi)
        if bw_map is None:
            return self.m_chi_prime**2
        return bw_map.Forward(random.Uniform(0.0, 1.0))

    def SampleFinalState(self, record, random):
        """Sample the physical narrow-width distribution: chi' production
        angle from dsigma/dQ2, virtuality s_pair from a Breit-Wigner, and an
        isotropic chi' -> chi + V1 decay in the chi' rest frame.  This matches
        the density reported by FinalStateProbability."""
        E_chi = record.primary_momentum[0]
        M = self.m_target
        m_chi = self.m_chi

        Q2 = self._sample_Q2(E_chi, random)
        if Q2 is None:
            return
        s_pair = self._sample_s_pair(E_chi, random)
        m_pair = math.sqrt(s_pair)

        s = m_chi**2 + M**2 + 2.0 * M * E_chi
        sqrt_s = math.sqrt(s)
        E_pair_cm = (s + s_pair - M**2) / (2.0 * sqrt_s)
        E_N_cm = (s - s_pair + M**2) / (2.0 * sqrt_s)
        p_out_cm = math.sqrt(max(E_pair_cm**2 - s_pair, 0.0))

        E_in_cm = (s + m_chi**2 - M**2) / (2.0 * sqrt_s)
        p_in_cm = math.sqrt(max(E_in_cm**2 - m_chi**2, 0.0))

        cos_cm = 1.0 - Q2 / (2.0 * p_in_cm * p_out_cm) if p_in_cm > 0 and p_out_cm > 0 else 0.0
        cos_cm = max(-1.0, min(1.0, cos_cm))

        gamma_cm = (E_chi + M) / sqrt_s
        beta_cm = math.sqrt(max(E_chi**2 - m_chi**2, 0.0)) / (E_chi + M)

        phi_cm = random.Uniform(0.0, 2.0 * math.pi)
        sin_cm = math.sqrt(max(1.0 - cos_cm**2, 0.0))

        px_cp = p_out_cm * sin_cm * math.cos(phi_cm)
        py_cp = p_out_cm * sin_cm * math.sin(phi_cm)
        pz_cp = p_out_cm * cos_cm

        E_cp_lab = gamma_cm * (E_pair_cm + beta_cm * pz_cp)
        pz_cp_lab = gamma_cm * (pz_cp + beta_cm * E_pair_cm)

        P_chi_prime = np.array([E_cp_lab, px_cp, py_cp, pz_cp_lab])

        p_decay = _two_body_p_cm(m_pair, m_chi, self.m_V1)
        cos_dec = random.Uniform(-1.0, 1.0)
        phi_dec = random.Uniform(0.0, 2.0 * math.pi)

        P_chi_out = _boost_to_lab(P_chi_prime, p_decay, cos_dec, phi_dec, m_chi)
        P_V1 = _boost_to_lab(P_chi_prime, p_decay, -cos_dec, phi_dec + math.pi, self.m_V1)

        E_N_lab = gamma_cm * (E_N_cm - beta_cm * pz_cp)
        pz_N_lab = gamma_cm * (-pz_cp + beta_cm * E_N_cm)
        P_N = np.array([E_N_lab, -px_cp, -py_cp, pz_N_lab])

        p_pi_dir = np.array(record.primary_momentum[1:])
        p_pi_mag = np.linalg.norm(p_pi_dir)
        if p_pi_mag > 1e-12:
            z_hat = p_pi_dir / p_pi_mag
            arb = np.array([0, 1, 0]) if abs(z_hat[1]) < 0.9 else np.array([1, 0, 0])
            x_hat = np.cross(z_hat, arb)
            x_hat /= np.linalg.norm(x_hat)
            y_hat = np.cross(z_hat, x_hat)
            R = np.column_stack([x_hat, y_hat, z_hat])

            for P in [P_chi_out, P_V1, P_N]:
                P[1:] = R @ P[1:]

        secondaries = record.get_secondary_particle_records()
        for sec in secondaries:
            pid = int(sec.type)
            if pid == self.pdgid_chi:
                sec.four_momentum = P_chi_out
                sec.mass = m_chi
            elif pid == self.pdgid_V1:
                sec.four_momentum = P_V1
                sec.mass = self.m_V1
            else:
                sec.four_momentum = P_N
                sec.mass = M
        return


# ===================================================================
#  ChiPrimeDecay  --  chi' -> chi + V1
# ===================================================================

class ChiPrimeDecay(_DecayModel):
    """
    Two-body decay chi' -> chi + V1.
    Width: Gamma = (g_D^2 / 48 pi) m_chi' lambda^{3/2}(1, r_chi^2, r_V^2)

    Isotropic in the chi' rest frame (SolidAngleRest 2-body): the authoring
    base derives the signature methods, the width overload pair, the isotropic
    1/(4 pi) FinalStateProbability, Topology/Measure, and the closure-by-
    construction Isotropic2BodyChannel sampler from total_width() /
    differential_width().
    """

    measure = _Measure.SolidAngleRest()
    daughter_index = 0

    def __init__(
        self,
        m_chi,
        m_chi_prime,
        m_V1,
        g_D,
        *,
        pdgid_chi_prime=5918,
        pdgid_chi=5917,
        pdgid_V1=5922,
        table_dir=None,
    ):
        _DecayModel.__init__(self)
        self.m_chi = m_chi
        self.m_chi_prime = m_chi_prime
        self.m_V1 = m_V1
        self.g_D = g_D

        self.pdgid_chi_prime = pdgid_chi_prime
        self.pdgid_chi = pdgid_chi
        self.pdgid_V1 = pdgid_V1

        self.parent = Particle.ParticleType(pdgid_chi_prime)
        self.daughters = (Particle.ParticleType(pdgid_chi),
                          Particle.ParticleType(pdgid_V1))

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)
        self._total_width = self._compute_width()

    def _compute_width(self):
        p = _two_body_p_cm(self.m_chi_prime, self.m_chi, self.m_V1)
        if p <= 0.0:
            return 0.0
        return self.g_D**2 * p**3 / (6.0 * math.pi * self.m_chi_prime**2)

    def total_width(self):
        return self._total_width

    def differential_width(self, record):
        if int(record.signature.primary_type) != self.pdgid_chi_prime:
            return 0.0
        return self._total_width / (4.0 * math.pi)

    def density_variables(self):
        return ["cos_theta"]

    def Convention(self):
        return _PhaseSpaceConvention.RestFrameSolidAngle

    def SecondaryMasses(self, secondary_types):
        return [self.m_chi, self.m_V1]

    def SecondaryHelicities(self, record):
        return [record.primary_helicity, 0]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other


# ===================================================================
#  DarkPhotonDecay  --  V1 -> e- e+
# ===================================================================

class DarkPhotonDecay(_DecayModel):
    """
    Two-body decay V1 -> e- e+.
    Width: Gamma = (alpha epsilon^2 m_V / 3) sqrt(1 - 4 m_e^2/m_V^2) (1 + 2 m_e^2/m_V^2)
    Angular distribution: dGamma/d(cos theta) ~ 1 + beta^2 cos^2(theta)

    Declared measure SolidAngleRest 2-body, but the rest-frame angular density
    is 1 + beta^2 cos^2(theta), NOT isotropic. differential_width carries that
    shape (the base forms FinalStateProbability = differential / total), and
    sample() is overridden with the matching rejection sampler so Sample and
    Density stay the same distribution.
    """

    measure = _Measure.SolidAngleRest()
    daughter_index = 0

    def __init__(
        self,
        m_V1,
        epsilon,
        *,
        pdgid_V1=5922,
        table_dir=None,
    ):
        _DecayModel.__init__(self)
        self.m_V1 = m_V1
        self.epsilon = epsilon
        self.pdgid_V1 = pdgid_V1

        self.parent = Particle.ParticleType(pdgid_V1)
        self.daughters = (Particle.ParticleType.EMinus,
                          Particle.ParticleType.EPlus)

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)
        self._total_width = self._compute_width()

    def _compute_width(self):
        mV = self.m_V1
        me = _M_ELECTRON
        if mV < 2.0 * me:
            return 0.0
        beta = math.sqrt(max(1.0 - (2.0 * me / mV)**2, 0.0))
        return (_ALPHA_EM * self.epsilon**2 * mV / 3.0) * beta * (1.0 + 2.0 * me**2 / mV**2)

    def total_width(self):
        return self._total_width

    def differential_width(self, record):
        if int(record.signature.primary_type) != self.pdgid_V1:
            return 0.0
        if self._total_width <= 0.0:
            return 0.0

        for idx, stype in enumerate(record.signature.secondary_types):
            if stype == Particle.ParticleType.EMinus:
                cos_theta = _cos_theta_in_parent_rest(
                    record.primary_momentum,
                    record.secondary_momenta[idx],
                )
                break
        else:
            cos_theta = 0.0

        p_cm = _two_body_p_cm(self.m_V1, _M_ELECTRON, _M_ELECTRON)
        e_cm = math.sqrt(p_cm**2 + _M_ELECTRON**2)
        beta_e = p_cm / e_cm if e_cm > 0.0 else 0.0
        norm = 4.0 * math.pi * (1.0 + beta_e**2 / 3.0)
        return self._total_width * (1.0 + beta_e**2 * cos_theta**2) / norm

    def density_variables(self):
        return ["cos_theta"]

    def Convention(self):
        return _PhaseSpaceConvention.RestFrameSolidAngle

    def SecondaryMasses(self, secondary_types):
        return [_M_ELECTRON, _M_ELECTRON]

    def SecondaryHelicities(self, record):
        return [0, 0]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other

    def sample(self, record, random):
        me = _M_ELECTRON
        p_cm = _two_body_p_cm(self.m_V1, me, me)
        beta = p_cm / math.sqrt(p_cm**2 + me**2) if p_cm > 0 else 0.0

        while True:
            cos_theta = random.Uniform(-1.0, 1.0)
            u = random.Uniform(0.0, 1.0 + beta**2)
            if u <= 1.0 + beta**2 * cos_theta**2:
                break

        phi = random.Uniform(0.0, 2.0 * math.pi)

        P_parent = np.array(record.primary_momentum)
        P_eminus = _boost_to_lab(P_parent, p_cm, cos_theta, phi, me)
        P_eplus = _boost_to_lab(P_parent, p_cm, -cos_theta, phi + math.pi, me)

        for sec in record.get_secondary_particle_records():
            if sec.type == Particle.ParticleType.EMinus:
                sec.four_momentum = P_eminus
                sec.mass = me
            elif sec.type == Particle.ParticleType.EPlus:
                sec.four_momentum = P_eplus
                sec.mass = me
        return


# ===================================================================
#  DarkPhotonToChiDecay  --  V1 -> chi chi_bar
# ===================================================================

class DarkPhotonToChiDecay(_DecayModel):
    """
    Two-body decay V1 -> chi chi_bar (dark matter pair production).

    Width: Gamma = (g_D^2 m_V / 12 pi) * beta^3
    where beta = sqrt(1 - 4 m_chi^2 / m_V^2).

    This is the dominant V1 decay when kinematically allowed
    (m_V > 2 m_chi).  The chi and chi_bar are both assigned
    the same PDG ID (chi) - the injector treats them identically
    and the stopping condition prevents re-scattering.

    Isotropic in the V1 rest frame (SolidAngleRest 2-body): the authoring
    base derives the signature methods, the width overload pair, the isotropic
    1/(4 pi) FinalStateProbability, Topology/Measure, and the closure-by-
    construction Isotropic2BodyChannel sampler from total_width() /
    differential_width().
    """

    measure = _Measure.SolidAngleRest()
    daughter_index = 0

    def __init__(
        self,
        m_V1,
        m_chi,
        g_D,
        *,
        pdgid_V1=5922,
        pdgid_chi=5917,
        table_dir=None,
    ):
        _DecayModel.__init__(self)
        self.m_V1 = m_V1
        self.m_chi = m_chi
        self.g_D = g_D
        self.pdgid_V1 = pdgid_V1
        self.pdgid_chi = pdgid_chi

        self.parent = Particle.ParticleType(pdgid_V1)
        self.daughters = (Particle.ParticleType(pdgid_chi),
                          Particle.ParticleType(pdgid_chi))

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)
        self._total_width = self._compute_width()

    def _compute_width(self):
        mV = self.m_V1
        mc = self.m_chi
        if mV < 2.0 * mc:
            return 0.0
        beta = math.sqrt(max(1.0 - (2.0 * mc / mV)**2, 0.0))
        return self.g_D**2 * mV / (12.0 * math.pi) * beta**3

    def total_width(self):
        return self._total_width

    def differential_width(self, record):
        if int(record.signature.primary_type) != self.pdgid_V1:
            return 0.0
        return self._total_width / (4.0 * math.pi)

    def density_variables(self):
        return ["cos_theta"]

    def Convention(self):
        return _PhaseSpaceConvention.RestFrameSolidAngle

    def SecondaryMasses(self, secondary_types):
        return [self.m_chi, self.m_chi]

    def SecondaryHelicities(self, record):
        return [record.primary_helicity, record.primary_helicity]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other


# ===================================================================
#  BiasedDarkPhotonToChiDecay  --  V1 -> chi chi_bar (cone-biased)
# ===================================================================

class BiasedDarkPhotonToChiDecay(_Decay):
    """
    Biased V1 -> chi chi_bar decay with energy-dependent cone.

    The cone opening angle adapts to the V1 energy: for highly boosted
    V1, chi is collimated and a narrow cone suffices.  For slow V1,
    chi spreads widely and a broader cone is needed.

    The cone half-angle is:
        theta_cone = max(theta_det, theta_max_chi * safety_factor)
    where theta_max_chi is the maximum lab-frame chi angle from the
    V1 rest-frame decay kinematics.
    """

    def __init__(
        self,
        m_V1,
        m_chi,
        g_D,
        *,
        detector_position=(0.0, 0.0, 0.0),
        cone_half_angle=None,
        detector_radius=2.0,
        safety_factor=1.5,
        pdgid_V1=5922,
        pdgid_chi=5917,
        table_dir=None,
    ):
        _Decay.__init__(self)
        self.m_V1 = m_V1
        self.m_chi = m_chi
        self.g_D = g_D
        self.pdgid_V1 = pdgid_V1
        self.pdgid_chi = pdgid_chi

        self.detector_position = np.array(detector_position, dtype=float)
        self.cone_half_angle = cone_half_angle
        self.detector_radius = detector_radius
        self.safety_factor = safety_factor

        self.table_dir = table_dir or "."
        os.makedirs(self.table_dir, exist_ok=True)

        mV, mc = m_V1, m_chi
        if mV < 2.0 * mc:
            self._total_width = 0.0
        else:
            beta = math.sqrt(max(1.0 - (2.0 * mc / mV)**2, 0.0))
            self._total_width = g_D**2 * mV / (12.0 * math.pi) * beta**3

        self._p_cm = _two_body_p_cm(m_V1, m_chi, m_chi)
        self._E_chi_rf = math.sqrt(self._p_cm**2 + m_chi**2)
        self._beta_chi = self._p_cm / self._E_chi_rf if self._E_chi_rf > 0 else 0.0

    def _chi_theta_max(self, E_V1):
        """Maximum lab-frame chi opening angle for a V1 with energy E_V1."""
        if E_V1 <= self.m_V1:
            return math.pi
        gamma = E_V1 / self.m_V1
        beta_V = math.sqrt(1.0 - 1.0 / gamma**2)
        if self._beta_chi >= beta_V:
            return math.pi
        # sin(theta_max) = p_cm / (gamma * E_chi_rf * beta_V)
        # (from the standard relativistic decay angle formula)
        sin_max = self._p_cm / (gamma * self._E_chi_rf * beta_V)
        return math.asin(min(sin_max, 1.0))

    def _get_cone_params(self, decay_vertex, E_V1=None):
        delta = self.detector_position - np.array(decay_vertex)
        dist = np.linalg.norm(delta)
        if dist < 1e-6:
            return np.array([0, 0, 1.0]), math.pi
        axis = delta / dist

        if self.cone_half_angle is not None:
            half_angle = self.cone_half_angle
        else:
            theta_det = math.atan2(self.detector_radius, dist)
            if E_V1 is not None:
                theta_chi = self._chi_theta_max(E_V1)
                half_angle = max(theta_det, theta_chi * self.safety_factor)
            else:
                half_angle = theta_det
            half_angle = min(half_angle, math.pi)

        return axis, half_angle

    def GetPossibleSignatures(self):
        sig = dataclasses.InteractionSignature()
        sig.primary_type = Particle.ParticleType(self.pdgid_V1)
        sig.target_type = Particle.ParticleType.Decay
        sig.secondary_types = [
            Particle.ParticleType(self.pdgid_chi),
            Particle.ParticleType(self.pdgid_chi),
        ]
        return [sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        if int(primary_type) == self.pdgid_V1:
            return self.GetPossibleSignatures()
        return []

    def TotalDecayWidthAllFinalStates(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_V1:
            return 0.0
        return self._total_width

    def TotalDecayWidth(self, arg1):
        primary = _primary_type(arg1)
        if int(primary) != self.pdgid_V1:
            return 0.0
        return self._total_width

    def DifferentialDecayWidth(self, record):
        if int(record.signature.primary_type) != self.pdgid_V1:
            return 0.0
        return self._total_width / (4.0 * math.pi)

    def FinalStateProbability(self, record):
        """Biased density: energy-dependent cone -> 1/Omega_cone(E_V1)."""
        E_V1 = record.primary_momentum[0]
        _, half_angle = self._get_cone_params(record.interaction_vertex, E_V1)
        omega = 2.0 * math.pi * (1.0 - math.cos(half_angle))
        return 1.0 / omega if omega > 0 else 0.0

    def DensityVariables(self):
        return ["lab_cone_chi"]

    def Convention(self):
        return _PhaseSpaceConvention.Custom

    def Topology(self):
        return _Topology.Decay2Body

    def Measure(self):
        return _Measure.Unspecified()

    def SecondaryMasses(self, secondary_types):
        return [self.m_chi, self.m_chi]

    def SecondaryHelicities(self, record):
        return [record.primary_helicity, record.primary_helicity]

    def save_to_table(self, table_subdir=None):
        pass

    def equal(self, other):
        return self is other

    def SampleFinalState(self, record, random):
        p_cm = self._p_cm
        P_parent = np.array(record.primary_momentum)
        E_V1 = P_parent[0]
        decay_vertex = np.array(record.interaction_vertex)

        cone_axis, half_angle = self._get_cone_params(decay_vertex, E_V1)

        cos_cone = random.Uniform(math.cos(half_angle), 1.0)
        phi_cone = random.Uniform(0.0, 2.0 * math.pi)
        sin_cone = math.sqrt(max(1.0 - cos_cone**2, 0.0))

        arb = np.array([0, 1, 0]) if abs(cone_axis[1]) < 0.9 else np.array([1, 0, 0])
        perp1 = np.cross(cone_axis, arb)
        perp1 /= np.linalg.norm(perp1)
        perp2 = np.cross(cone_axis, perp1)

        chi_dir_lab = (cos_cone * cone_axis
                       + sin_cone * math.cos(phi_cone) * perp1
                       + sin_cone * math.sin(phi_cone) * perp2)

        cos_theta = float(chi_dir_lab[2])
        phi = math.atan2(float(chi_dir_lab[1]), float(chi_dir_lab[0]))

        P_chi1 = _boost_to_lab(P_parent, p_cm, cos_theta, phi, self.m_chi)
        P_chi2 = _boost_to_lab(P_parent, p_cm, -cos_theta, phi + math.pi, self.m_chi)

        secondaries = record.get_secondary_particle_records()
        secondaries[0].four_momentum = P_chi1
        secondaries[0].mass = self.m_chi
        secondaries[1].four_momentum = P_chi2
        secondaries[1].mass = self.m_chi
        return


# ===================================================================
#  Flux construction
# ===================================================================

def compute_chi_flux(
    m_meson,
    m_lepton,
    m_V1,
    m_chi,
    m_chi_prime,
    g_D,
    epsilon,
    flux_tag,
    min_energy,
    max_energy,
    physically_normalized=True,
):
    """
    Construct the chi (dark matter) flux at the detector by folding:
        neutrino flux -> parent meson energy -> three-body BR(meson -> l nu V1)
        -> V1 -> chi chi' kinematics -> chi energy spectrum.

    Returns a siren.distributions.TabulatedFluxDistribution.
    """
    import siren

    raw_flux = siren.utilities.load_flux(
        "PionKaon",
        tag=flux_tag,
        physically_normalized=physically_normalized,
    )

    available = m_meson - m_lepton
    if available <= m_V1:
        raise RuntimeError(
            "Channel kinematically forbidden: "
            f"m_meson={m_meson*1e3:.1f} MeV, m_lepton={m_lepton*1e3:.1f} MeV, "
            f"m_V1={m_V1*1e3:.1f} MeV"
        )

    E_nu_rf = (m_meson**2 - m_lepton**2) / (2.0 * m_meson)
    nu_to_meson = m_meson / E_nu_rf

    alpha_D = g_D**2 / (4.0 * math.pi)
    x = m_lepton / m_meson
    y = m_V1 / m_meson

    if (1.0 - x - y) <= 0:
        br_ratio = 0.0
    else:
        num = (1.0 - y**2)**2 * (1.0 + 2.0 * y**2)
        den = (1.0 - x**2)**2
        g_ps = num / den if den > 0 else 0.0
        br_ratio = 2.0 * (alpha_D / _ALPHA_EM) * epsilon**2 * g_ps

    E_chi_rf = (m_V1**2 + m_chi**2) / (2.0 * m_V1) if m_V1 > 0 else 0.0
    E_V_rest = (m_meson**2 + m_V1**2 - m_lepton**2) / (2.0 * m_meson)

    nu_energies = list(raw_flux.GetEnergyNodes())
    chi_energies = []
    chi_flux_vals = []

    for E_nu in nu_energies:
        E_meson = E_nu * nu_to_meson
        if E_meson < m_meson:
            continue

        gamma_meson = E_meson / m_meson
        E_V_lab = gamma_meson * E_V_rest

        gamma_V1 = E_V_lab / m_V1 if m_V1 > 0 else 1.0
        E_chi = gamma_V1 * E_chi_rf

        if E_chi < min_energy or E_chi > max_energy:
            continue

        nu_flux_at_E = raw_flux.SamplePDF(E_nu)
        chi_flux_at_E = nu_flux_at_E * br_ratio * 0.5

        chi_energies.append(E_chi)
        chi_flux_vals.append(chi_flux_at_E)

    if not chi_energies:
        chi_energies = [min_energy, max_energy]
        chi_flux_vals = [0.0, 0.0]

    return siren.distributions.TabulatedFluxDistribution(
        min_energy, max_energy, chi_energies, chi_flux_vals, physically_normalized
    )


def compute_chi_flux_from_dk2nu(
    dk2nu_data,
    parent_pdg,
    m_meson,
    m_lepton,
    m_V1,
    m_chi,
    m_chi_prime,
    g_D,
    epsilon,
    min_energy,
    max_energy,
    n_bins=100,
    physically_normalized=True,
):
    """
    Construct the chi flux from dk2nu parent meson data.

    Instead of using tabulated neutrino flux tables, this reads the
    parent meson energies and weights directly from dk2nu simulation
    output, giving a more accurate flux prediction that properly
    accounts for the beamline geometry and focusing.

    Parameters
    ----------
    dk2nu_data : dict
        Output of Dk2nuReader.read_dk2nu().
    parent_pdg : int
        PDG code of parent meson (211 for pi+, 321 for K+, etc.)
    m_meson, m_lepton, m_V1, m_chi, m_chi_prime : float
        Masses in GeV.
    g_D, epsilon : float
        Dark coupling and kinetic mixing.
    min_energy, max_energy : float
        Energy range for the output flux [GeV].
    n_bins : int
        Number of bins in the output flux table.

    Returns
    -------
    siren.distributions.TabulatedFluxDistribution
    """
    import siren

    mask = dk2nu_data["ptype"] == parent_pdg
    E_meson = dk2nu_data["E"][mask]
    weights = dk2nu_data["nimpwt"][mask]

    available = m_meson - m_lepton
    if available <= m_V1 or len(E_meson) == 0:
        energies = [min_energy, max_energy]
        flux_arr = [0.0, 0.0]
        return siren.distributions.TabulatedFluxDistribution(
            min_energy, max_energy, energies, flux_arr, physically_normalized
        )

    alpha_D = g_D**2 / (4.0 * math.pi)
    x = m_lepton / m_meson
    y = m_V1 / m_meson
    if (1.0 - x - y) <= 0:
        br_ratio = 0.0
    else:
        num = (1.0 - y**2)**2 * (1.0 + 2.0 * y**2)
        den = (1.0 - x**2)**2
        g_ps = num / den if den > 0 else 0.0
        br_ratio = 2.0 * (alpha_D / _ALPHA_EM) * epsilon**2 * g_ps

    E_chi_rf = (m_V1**2 + m_chi**2) / (2.0 * m_V1) if m_V1 > 0 else 0.0
    E_V_rest = (m_meson**2 + m_V1**2 - m_lepton**2) / (2.0 * m_meson)

    gamma_mesons = E_meson / m_meson
    E_V_lab = gamma_mesons * E_V_rest
    gamma_V1 = E_V_lab / m_V1 if m_V1 > 0 else np.ones_like(E_V_lab)
    E_chi = gamma_V1 * E_chi_rf

    chi_weights = weights * br_ratio * 0.5

    E_edges = np.linspace(min_energy, max_energy, n_bins + 1)
    hist, _ = np.histogram(E_chi, bins=E_edges, weights=chi_weights)
    dE = E_edges[1] - E_edges[0]
    E_centers = 0.5 * (E_edges[:-1] + E_edges[1:])

    pot = dk2nu_data.get("pot", 0.0)
    if pot > 0:
        flux_vals = hist / (dE * pot)
    else:
        flux_vals = hist / dE

    return siren.distributions.TabulatedFluxDistribution(
        min_energy, max_energy,
        list(E_centers), list(flux_vals),
        physically_normalized,
    )
