"""Contract tests for the Dutta-Kim vector-portal chain.

Validates per-vertex channel behavior (momentum conservation, density
positivity, closure) and end-to-end chain assembly for both on-shell
(5-vertex) and off-shell (4-vertex) topologies.
"""

import copy
import math
import os
import sys

import numpy as np
import pytest

siren = pytest.importorskip("siren")
from siren import _util, dataclasses, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

# ------------------------------------------------------------------ #
#  Constants                                                          #
# ------------------------------------------------------------------ #

M_PION = 0.13957039
M_MUON = 0.10565837
M_ELECTRON = 0.000511
M_V1 = 0.017
M_CHI = 0.008
M_CHI_PRIME = 0.050
M_V2 = 0.200
M_ARGON40 = 37.215

G_D = 1.0
EPSILON_1 = 7e-5
EPSILON_2 = 1e-4
G_MU = 1e-3

PT = lambda pdg: dataclasses.Particle.ParticleType(pdg)

# ------------------------------------------------------------------ #
#  Module-scoped fixtures                                             #
# ------------------------------------------------------------------ #

_PROC_DIR = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")


@pytest.fixture(scope="module")
def vp():
    return _util.load_module("chain_VP", os.path.join(_PROC_DIR, "VectorPortal.py"))


@pytest.fixture(scope="module")
def meson():
    return _util.load_module("chain_Meson", os.path.join(_PROC_DIR, "MesonProduction.py"))


@pytest.fixture(scope="module")
def target():
    p = siren.geometry.Placement(siren.math.Vector3D(0, 0, 50))
    return siren.geometry.Box(p, 4.0, 4.0, 4.0)


@pytest.fixture(scope="module")
def rng():
    return siren.utilities.SIREN_random(12345)


# ------------------------------------------------------------------ #
#  Helpers                                                            #
# ------------------------------------------------------------------ #

def _decay_record(sig, parent_mass, energy, sec_masses, vertex=(0, 0, -100)):
    """Build a template record for a decay vertex."""
    rec = dataclasses.InteractionRecord()
    rec.signature = sig
    rec.primary_mass = parent_mass
    pz = math.sqrt(max(energy**2 - parent_mass**2, 0.0))
    rec.primary_momentum = [energy, 0, 0, pz]
    rec.primary_initial_position = [0, 0, 0]
    rec.interaction_vertex = list(vertex)
    rec.secondary_masses = list(sec_masses)
    rec.secondary_momenta = [[0, 0, 0, 0] for _ in sec_masses]
    rec.secondary_helicities = [0] * len(sec_masses)
    return rec


def _scatter_record(sig, primary_mass, energy, target_mass, sec_masses,
                    vertex=(0, 0, -100)):
    """Build a template record for a scattering vertex."""
    rec = _decay_record(sig, primary_mass, energy, sec_masses, vertex)
    rec.target_mass = target_mass
    return rec


def _four_momentum_sum(record):
    """Sum of secondary 4-momenta."""
    total = np.zeros(4)
    for p in record.secondary_momenta:
        total += np.array(p)
    return total


def _check_4mom_conservation(record, tol=1e-6):
    """Assert 4-momentum is conserved at this vertex."""
    p_in = np.array(record.primary_momentum)
    if record.target_mass > 0:
        p_in = p_in + np.array([record.target_mass, 0, 0, 0])
    p_out = _four_momentum_sum(record)
    np.testing.assert_allclose(p_out, p_in, atol=tol,
                               err_msg="4-momentum not conserved")


def _closure_test(mc, ref_channel, template, rng, N=500, tol_lo=0.7, tol_hi=1.3):
    """Run a closure test: mean(ref_density / mc_density) ~ 1.0."""
    weights = []
    for _ in range(N):
        r = copy.deepcopy(template)
        mc.Sample(rng, None, r)
        d_ref = ref_channel.Density(None, r)
        d_mc = mc.Density(None, r)
        if d_mc > 0 and d_ref > 0:
            weights.append(d_ref / d_mc)
    assert len(weights) > N * 0.3, (
        f"Too few valid weights: {len(weights)}/{N}")
    mean_w = sum(weights) / len(weights)
    assert tol_lo < mean_w < tol_hi, (
        f"Closure failed: mean={mean_w:.4f}, expected ~1.0")


def _mc(channels, weights):
    m = injection.MultiChannelPhaseSpace()
    m.channels = channels
    m.weights = weights
    return m


# ------------------------------------------------------------------ #
#  Vertex 0: pi+ -> mu+ nu V1                                        #
# ------------------------------------------------------------------ #

class TestVertex0_PionDecay:

    @pytest.fixture()
    def pion_decay(self, meson):
        return meson.MesonThreeBodySIRENDecay(
            m_mediator=M_V1, g_mu=G_MU,
            pdgid_meson=211, pdgid_lepton=-13,
            pdgid_neutrino=14, pdgid_mediator=5922)

    @pytest.fixture()
    def pion_sig(self, pion_decay):
        return pion_decay.GetPossibleSignatures()[0]

    @pytest.fixture()
    def pion_record(self, pion_sig):
        return _decay_record(pion_sig, M_PION, 2.0,
                             [M_MUON, 0.0, M_V1])

    def test_direct_3body_samples_conserve_momentum(
            self, pion_record, target, rng):
        ch = injection.DetectorDirected3BodyChannel(
            target, directed_index=2,
            mass_mode=injection.InvariantMassMode.Uniform,
            topology=injection.PhaseSpaceTopology.Decay3Body)
        for _ in range(100):
            r = copy.deepcopy(pion_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-6)

    def test_direct_3body_density_positive(
            self, pion_record, target, rng):
        ch = injection.DetectorDirected3BodyChannel(
            target, directed_index=2,
            mass_mode=injection.InvariantMassMode.Uniform,
            topology=injection.PhaseSpaceTopology.Decay3Body)
        for _ in range(20):
            r = copy.deepcopy(pion_record)
            ch.Sample(rng, None, r)
            assert ch.Density(None, r) > 0

    def test_direct_3body_v1_points_toward_target(
            self, pion_record, target, rng):
        ch = injection.DetectorDirected3BodyChannel(
            target, directed_index=2,
            mass_mode=injection.InvariantMassMode.Uniform,
            topology=injection.PhaseSpaceTopology.Decay3Body)
        hits = 0
        N = 50
        for _ in range(N):
            r = copy.deepcopy(pion_record)
            ch.Sample(rng, None, r)
            p_v1 = np.array(r.secondary_momenta[2][1:])
            mag = np.linalg.norm(p_v1)
            if mag > 0:
                cos_z = p_v1[2] / mag
                if cos_z > 0:
                    hits += 1
        assert hits > N * 0.5, (
            f"V1 should mostly point toward target (+z); got {hits}/{N}")

    def test_direct_vs_physical_closure(
            self, pion_decay, pion_sig, pion_record, target, rng):
        phys = injection.PhysicalDecayChannel(pion_decay, pion_sig)
        directed = injection.DetectorDirected3BodyChannel(
            target, directed_index=2,
            mass_mode=injection.InvariantMassMode.Uniform,
            topology=injection.PhaseSpaceTopology.Decay3Body)
        mc = _mc([phys, directed], [0.02, 0.98])
        _closure_test(mc, phys, pion_record, rng, N=2000,
                      tol_lo=0.5, tol_hi=1.5)


# ------------------------------------------------------------------ #
#  Vertex 1: V1(5922) -> chi chi                                     #
# ------------------------------------------------------------------ #

class TestVertex1_V1Decay:

    @pytest.fixture()
    def v1_decay(self, vp):
        return vp.DarkPhotonToChiDecay(M_V1, M_CHI, G_D,
                                       pdgid_V1=5922, pdgid_chi=5917)

    @pytest.fixture()
    def v1_sig(self, v1_decay):
        return v1_decay.GetPossibleSignatures()[0]

    @pytest.fixture()
    def v1_record(self, v1_sig):
        return _decay_record(v1_sig, M_V1, 0.5, [M_CHI, M_CHI])

    def test_directed_2body_chi_conserves_momentum(
            self, v1_record, target, rng):
        ch = injection.DetectorDirected2BodyChannel(target, 0)
        for _ in range(100):
            r = copy.deepcopy(v1_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-6)

    def test_directed_2body_chi_points_toward_target(
            self, v1_record, target, rng):
        ch = injection.DetectorDirected2BodyChannel(target, 0)
        hits = 0
        N = 50
        for _ in range(N):
            r = copy.deepcopy(v1_record)
            ch.Sample(rng, None, r)
            p_chi = np.array(r.secondary_momenta[0][1:])
            mag = np.linalg.norm(p_chi)
            if mag > 0 and p_chi[2] / mag > 0:
                hits += 1
        assert hits > N * 0.5

    def test_directed_2body_closure(
            self, v1_decay, v1_sig, v1_record, target, rng):
        phys = injection.PhysicalDecayChannel(v1_decay, v1_sig)
        directed = injection.DetectorDirected2BodyChannel(target, 0)
        mc = _mc([phys, directed], [0.5, 0.5])
        _closure_test(mc, phys, v1_record, rng, N=500,
                      tol_lo=0.7, tol_hi=1.3)


# ------------------------------------------------------------------ #
#  Vertex 2 on-shell: chi + Ar -> chi' + Ar                          #
# ------------------------------------------------------------------ #

class TestVertex2_ChiUpscatter:

    @pytest.fixture()
    def upscatter(self, vp):
        return vp.VectorPortalUpscatteringXS(
            M_CHI, M_CHI_PRIME, M_V2, G_D, EPSILON_2,
            pdgid_chi=5917, pdgid_chi_prime=5918,
            nuclear_pdgid=1000180400, nuclear_mass=M_ARGON40, A=40, Z=18)

    @pytest.fixture()
    def chi_sig(self, upscatter):
        return upscatter.GetPossibleSignatures()[0]

    @pytest.fixture()
    def chi_record(self, upscatter, chi_sig):
        sec_masses = upscatter.SecondaryMasses(chi_sig.secondary_types)
        return _scatter_record(chi_sig, M_CHI, 1.0, M_ARGON40, sec_masses)

    def test_directed_scattering_conserves_momentum(
            self, chi_record, target, rng):
        ch = injection.DetectorDirectedScatteringChannel(
            target, directed_index=0,
            variable=injection.ScatteringVariable.Q2)
        for _ in range(50):
            r = copy.deepcopy(chi_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-4)

    def test_directed_scattering_closure(
            self, upscatter, chi_sig, chi_record, target, rng):
        phys = injection.PhysicalCrossSectionChannel(upscatter, chi_sig)
        directed = injection.DetectorDirectedScatteringChannel(
            target, directed_index=0,
            variable=injection.ScatteringVariable.Q2)
        mc = _mc([phys, directed], [0.5, 0.5])
        _closure_test(mc, phys, chi_record, rng, N=500,
                      tol_lo=0.5, tol_hi=1.5)


# ------------------------------------------------------------------ #
#  Vertex 2 off-shell: chi + Ar -> chi + V1_sig + Ar                 #
# ------------------------------------------------------------------ #

class TestVertex2_OffshellScatter:

    @pytest.fixture()
    def offshell_xs(self, vp):
        return vp.VectorPortalOffShellXS(
            M_CHI, M_CHI_PRIME, M_V1, M_V2, G_D, EPSILON_2,
            pdgid_V1=5923)

    @pytest.fixture()
    def offshell_sig(self, offshell_xs):
        return offshell_xs.GetPossibleSignatures()[0]

    @pytest.fixture()
    def offshell_record(self, offshell_xs, offshell_sig):
        sec_masses = offshell_xs.SecondaryMasses(offshell_sig.secondary_types)
        return _scatter_record(offshell_sig, M_CHI, 1.0, M_ARGON40,
                               sec_masses)

    def test_offshell_3body_conserves_momentum(
            self, offshell_record, target, rng):
        _P_STAR = injection.TwoBodyRestMomentum(M_CHI_PRIME, M_CHI, M_V1)
        CHI_PRIME_WIDTH = G_D**2 * _P_STAR**3 / (6.0 * math.pi * M_CHI_PRIME**2)
        ch = injection.DetectorDirected3BodyChannel(
            target,
            spectator_index=2, pair_first_index=0, pair_second_index=1,
            directed_pair_index=0,
            mass_mode=injection.InvariantMassMode.BreitWigner,
            resonance_mass=M_CHI_PRIME,
            resonance_width=CHI_PRIME_WIDTH,
            topology=injection.PhaseSpaceTopology.Scatter2to3)
        for _ in range(50):
            r = copy.deepcopy(offshell_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-4)

    def test_offshell_3body_closure(
            self, offshell_xs, offshell_sig, offshell_record, target, rng):
        _P_STAR = injection.TwoBodyRestMomentum(M_CHI_PRIME, M_CHI, M_V1)
        CHI_PRIME_WIDTH = G_D**2 * _P_STAR**3 / (6.0 * math.pi * M_CHI_PRIME**2)
        phys = injection.PhysicalCrossSectionChannel(offshell_xs, offshell_sig)
        directed = injection.DetectorDirected3BodyChannel(
            target,
            spectator_index=2, pair_first_index=0, pair_second_index=1,
            directed_pair_index=0,
            mass_mode=injection.InvariantMassMode.BreitWigner,
            resonance_mass=M_CHI_PRIME,
            resonance_width=CHI_PRIME_WIDTH,
            topology=injection.PhaseSpaceTopology.Scatter2to3)
        mc = _mc([phys, directed], [0.5, 0.5])
        _closure_test(mc, phys, offshell_record, rng, N=1000,
                      tol_lo=0.5, tol_hi=1.5)


# ------------------------------------------------------------------ #
#  Vertex 3: chi'(5918) -> chi(5917) + V1_sig(5923)                  #
# ------------------------------------------------------------------ #

class TestVertex3_ChiPrimeDecay:

    @pytest.fixture()
    def chi_prime_decay(self, vp):
        return vp.ChiPrimeDecay(M_CHI, M_CHI_PRIME, M_V1, G_D,
                                 pdgid_chi_prime=5918, pdgid_chi=5917,
                                 pdgid_V1=5923)

    @pytest.fixture()
    def chip_sig(self, chi_prime_decay):
        return chi_prime_decay.GetPossibleSignatures()[0]

    @pytest.fixture()
    def chip_record(self, chip_sig):
        return _decay_record(chip_sig, M_CHI_PRIME, 0.12,
                             [M_CHI, M_V1])

    def test_directed_2body_v1sig_conserves_momentum(
            self, chip_record, target, rng):
        ch = injection.DetectorDirected2BodyChannel(target, 1)
        for _ in range(100):
            r = copy.deepcopy(chip_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-6)

    def test_directed_2body_closure(
            self, chi_prime_decay, chip_sig, chip_record, target, rng):
        phys = injection.PhysicalDecayChannel(chi_prime_decay, chip_sig)
        directed = injection.DetectorDirected2BodyChannel(target, 1)
        mc = _mc([phys, directed], [0.5, 0.5])
        _closure_test(mc, phys, chip_record, rng, N=500,
                      tol_lo=0.7, tol_hi=1.3)


# ------------------------------------------------------------------ #
#  Vertex 4: V1_sig(5923) -> e+ e-                                   #
# ------------------------------------------------------------------ #

class TestVertex4_V1SignalDecay:

    @pytest.fixture()
    def visible_decay(self, vp):
        return vp.DarkPhotonDecay(M_V1, EPSILON_1, pdgid_V1=5923)

    @pytest.fixture()
    def vis_sig(self, visible_decay):
        return visible_decay.GetPossibleSignatures()[0]

    @pytest.fixture()
    def vis_record(self, vis_sig):
        return _decay_record(vis_sig, M_V1, 0.5,
                             [M_ELECTRON, M_ELECTRON])

    def test_isotropic_2body_conserves_momentum(
            self, vis_record, rng):
        ch = injection.Isotropic2BodyChannel(0)
        for _ in range(100):
            r = copy.deepcopy(vis_record)
            ch.Sample(rng, None, r)
            _check_4mom_conservation(r, tol=1e-6)

    def test_isotropic_closure(
            self, visible_decay, vis_sig, vis_record, rng):
        phys = injection.PhysicalDecayChannel(visible_decay, vis_sig)
        iso = injection.Isotropic2BodyChannel(0)
        mc = _mc([phys, iso], [0.5, 0.5])
        _closure_test(mc, iso, vis_record, rng, N=500,
                      tol_lo=0.8, tol_hi=1.2)


# ------------------------------------------------------------------ #
#  On-shell chain assembly                                            #
# ------------------------------------------------------------------ #

_EXAMPLE4_DIR = os.path.join(
    _util.resource_package_dir(), os.pardir, "resources", "examples", "example4")


@pytest.fixture(scope="module")
def chain_module():
    sys.path.insert(0, os.path.abspath(_EXAMPLE4_DIR))
    mod = _util.load_module(
        "chain_DuttaKim",
        os.path.join(_EXAMPLE4_DIR, "DuttaKim_SBND_full_chain.py"))
    return mod


def _build_onshell_events(chain_module, n_events=12, seed=99):
    """Build the small single-target on-shell chain; return (events, weighter)."""
    detector_model = siren.utilities.load_detector("SBN", detector="SBND")
    fiducial = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
        4.0, 4.0, 4.0)

    models = chain_module.build_onshell_models()
    primary, secondaries = chain_module.build_vertices(models, fiducial)

    injector = Injector(detector=detector_model, primary=primary,
                        secondaries=secondaries, events=n_events, seed=seed)
    events = injector.generate(n_events, on_shortfall="warn")
    weighter = Weighter(injector, primary_physical=primary.physical)
    return events, weighter


class TestOnShellChain:

    def test_chain_models_load(self, chain_module):
        models = chain_module.build_onshell_models()
        sec = models["secondary_interactions"]
        assert PT(5922) in sec
        assert PT(5917) in sec
        assert PT(5918) in sec
        assert PT(5923) in sec
        assert models["pion_decay"] is not None

    def test_chain_vertices_compile(self, chain_module):
        detector_model = siren.utilities.load_detector("SBN", detector="SBND")
        fiducial = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            4.0, 4.0, 4.0)
        models = chain_module.build_onshell_models()
        primary, secondaries = chain_module.build_vertices(models, fiducial)
        assert [v._resolved_particle for v in secondaries] == [
            PT(5922), PT(5917), PT(5918), PT(5923)]
        for is_primary, v in [(True, primary)] + [(False, s)
                                                  for s in secondaries]:
            proc = v.compile(is_primary=is_primary, detector=detector_model)
            assert proc is not None

    def test_chain_end_to_end(self, chain_module):
        events, weighter = _build_onshell_events(chain_module, n_events=10,
                                                 seed=99)
        assert len(events) > 0, "No events generated"
        for ev in events:
            assert len(list(ev.tree)) >= 2
            w = weighter(ev)
            assert math.isfinite(w), f"Non-finite weight: {w}"
            assert w > 0, f"Non-positive weight: {w}"

    def test_breakdown_sum_invariant(self, chain_module):
        """Product of per-vertex physical/generation factors equals the total
        weight up to one global 1/N_gen normalization, identical across events."""
        events, weighter = _build_onshell_events(chain_module, n_events=12, seed=99)
        assert len(events) > 0, "No events generated"

        offsets = []
        for ev in events:
            bd = weighter.explain(ev)
            if not (math.isfinite(bd.total) and bd.total > 0):
                continue
            if not all(v.is_ok for v in bd.vertices):
                continue
            log_prod = sum(math.log(v.physical / v.generation)
                           for v in bd.vertices)
            offsets.append(log_prod - math.log(bd.total))

        assert len(offsets) >= 2, "Need >= 2 clean events to check constancy"
        spread = max(offsets) - min(offsets)
        assert spread < 1e-9, (
            f"per-vertex/total offset is not constant (spread={spread:.3e}); "
            "breakdown() is missing an event-varying factor")
        # The constant is the global 1/N_gen normalization (offset = ln N_gen).
        assert offsets[0] > 0.0


# ------------------------------------------------------------------ #
#  Off-shell chain assembly                                           #
# ------------------------------------------------------------------ #

def _build_offshell_vertices(chain_module, vp, fiducial):
    """Off-shell (4-vertex, virtual chi') chain as spec-form vertices."""
    from siren import channels, expand

    models = chain_module.build_onshell_models()
    pion_decay = models["pion_decay"]
    v1_to_chi = models["models"]["v1_to_chi"]
    offshell_xs = vp.VectorPortalOffShellXS(
        M_CHI, M_CHI_PRIME, M_V1, M_V2, G_D, EPSILON_2, pdgid_V1=5923)
    visible_decay = models["models"]["visible_decay"]

    sx = channels.PairMass.tabulated(*chain_module.build_sX_cdf(pion_decay))
    p_star = injection.TwoBodyRestMomentum(M_CHI_PRIME, M_CHI, M_V1)
    chip_width = G_D ** 2 * p_star ** 3 / (6.0 * math.pi * M_CHI_PRIME ** 2)
    sv = distributions.SecondaryPhysicalVertexDistribution

    primary = siren.Vertex(
        PT(211), pion_decay,
        distributions=[
            distributions.PrimaryMass(M_PION),
            distributions.Monoenergetic(2.0),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
            distributions.PointSourcePositionDistribution([0, 0, 0], 300.0),
        ],
        physical=[
            distributions.Monoenergetic(2.0),
            distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        ],
        kinematics=0.98 * channels.toward_3body(2, fiducial, strategy="direct",
                                                pair_mass=sx)
                   + 0.02 * channels.physical(),
        expand=(expand.child("V1_prod"),))

    v1 = siren.Vertex(
        "V1_prod", v1_to_chi, position=sv(),
        kinematics=0.98 * channels.toward(0, fiducial)
                   + 0.02 * channels.physical(),
        expand=(expand.child("chi", index=0),))

    chi = siren.Vertex(
        "chi", offshell_xs, position=sv(),
        # 2->3 scatter: direct V1_sig (pair index 1) at the detector with a
        # Breit-Wigner pair mass on the virtual chi' resonance
        kinematics=0.98 * channels.toward_3body(
            1, fiducial, spectator=2, strategy="recursive",
            pair_mass=channels.PairMass.breit_wigner(M_CHI_PRIME, chip_width))
                   + 0.02 * channels.physical(),
        expand=(expand.child("V1_sig"),))

    v1s = siren.Vertex(
        "V1_sig", visible_decay, position=sv(),
        kinematics=0.5 * channels.physical() + 0.5 * channels.isotropic(0),
        expand=(expand.depth_below(0),))

    return primary, (v1, chi, v1s)


class TestOffShellChain:

    def test_chain_vertices_build(self, chain_module, vp):
        fiducial = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            4.0, 4.0, 4.0)
        primary, secondaries = _build_offshell_vertices(chain_module, vp,
                                                        fiducial)
        particles = [v._resolved_particle for v in secondaries]
        assert particles == [PT(5922), PT(5917), PT(5923)]
        assert PT(5918) not in particles

    def test_chain_vertices_compile(self, chain_module, vp):
        detector_model = siren.utilities.load_detector("SBN", detector="SBND")
        fiducial = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            4.0, 4.0, 4.0)
        primary, secondaries = _build_offshell_vertices(chain_module, vp,
                                                        fiducial)
        for is_primary, v in [(True, primary)] + [(False, s)
                                                  for s in secondaries]:
            proc = v.compile(is_primary=is_primary, detector=detector_model)
            assert proc is not None

    def test_chain_end_to_end(self, chain_module, vp):
        detector_model = siren.utilities.load_detector("SBN", detector="SBND")
        fiducial = siren.geometry.Box(
            siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)),
            4.0, 4.0, 4.0)
        primary, secondaries = _build_offshell_vertices(chain_module, vp,
                                                        fiducial)

        injector = Injector(detector=detector_model, primary=primary,
                            secondaries=secondaries, events=10, seed=77)
        events = injector.generate(10, on_shortfall="warn")
        weighter = Weighter(injector, primary_physical=primary.physical)

        assert len(events) > 0, "No events generated"
        for ev in events:
            assert len(list(ev.tree)) >= 2
            w = weighter(ev)
            assert math.isfinite(w), f"Non-finite weight: {w}"
            assert w > 0, f"Non-positive weight: {w}"
