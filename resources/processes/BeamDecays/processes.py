"""
Standard-model beam decays for dk2nu-driven SIREN injection.

Provides the decay channels that dominate conventional neutrino beams:

    MesonTwoBodyLeptonicDecay   pi+- -> mu nu,  K+- -> mu nu
    MuonThreeBodyDecay          mu-  -> e- nubar_e nu_mu   (and mu+ mirror)

Both are built on the siren.DecayModel authoring base, so injection uses
them like any other decay model: the two-body decay samples isotropically
in the parent rest frame (closure by construction for the SolidAngleRest
measure), and the muon decay samples its polarized matrix element over
the Recursive2Body measure through siren.three_body. The module depends on
siren and numpy, not on other process resources.

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

import siren
from siren.dataclasses import Particle
from siren import three_body

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
    and sample_isotropic implements the normalized 1/(4 pi) density
    in the declared SolidAngleRest measure.
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

    def sample(self, record, random):
        self.sample_isotropic(record, random)

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
        return three_body.dalitz_band(self.m_M, self.m_l, self.m_phi, E_nu)

    def total_width(self):
        """Three-body width [GeV] from the shared Dalitz integral."""
        return three_body.dalitz_width(self)


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
        self._max_matel = three_body.find_max_weight(self._physics)

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
        density = three_body.final_state_probability(
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
            P_phi_rest = three_body.boost_to_rest(
                record.primary_momentum, record.secondary_momenta[idx])
            density *= self._polarization_factor(pol, P_phi_rest)
        return density * self._total

    def density_variables(self):
        return ["s_pair", "cos_theta_sub"]

    def sample(self, record, random):
        e_nu_rf, e_phi_rf = three_body.sample_energies(
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
            p_nu, p_e, p_phi = three_body.build_rest_momenta(
                self._physics, e_nu_rf, e_phi_rf, random)
            if pol is None:
                break
            f = self._polarization_factor(pol, p_phi)
            if random.Uniform(0.0, bound) <= f:
                break
        for sec in record.get_secondary_particle_records():
            pid = int(sec.type)
            if pid == self.pdgid_electron:
                sec.four_momentum = three_body.boost_four_vector(parent, p_e, M_MU)
                sec.mass = M_E
            elif pid == self.pdgid_nu_mu:
                sec.four_momentum = three_body.boost_four_vector(parent, p_nu, M_MU)
                sec.mass = 0.0
            elif pid == self.pdgid_nu_e:
                sec.four_momentum = three_body.boost_four_vector(parent, p_phi, M_MU)
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
