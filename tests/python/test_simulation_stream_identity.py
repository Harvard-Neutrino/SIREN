"""BLOCKER 2: the facade default mixture is [physical@0.01, directed@0.99].

The compiled default biasing must reproduce the legacy channel order and the
0.01/0.99 split exactly, so a fixed-seed event stream is unchanged. Channel
selection walks cumulative weights in stored order, so both order and value
are load-bearing.
"""
import pytest

import siren
from siren.Simulation import Simulation
from siren import interactions

PT = siren.dataclasses.ParticleType


def _two_body_signature():
    sig = siren.dataclasses.InteractionSignature()
    sig.primary_type = PT.N4
    sig.target_type = PT.Decay
    sig.secondary_types = [PT.NuLight, PT.Gamma]
    return sig


def _fid():
    return siren.geometry.Box(widths=[2.0, 2.0, 2.0], center=[0.0, 0.0, 5.0])


def test_default_builder_is_physical_first_point_99():
    sig = _two_body_signature()
    xs = interactions.DummyCrossSection()
    mc = Simulation._build_phase_space_for_signature(
        _fid(), sig, xs, PT.NuLight, fraction=0.99)

    assert len(mc.channels) == 2
    # Physical channel FIRST, directed SECOND (the legacy order).
    assert mc.channels[0].Name().startswith("Physical")
    assert mc.channels[1].Name().startswith("DetectorDirected")
    assert list(mc.weights) == pytest.approx([0.01, 0.99], abs=1e-12)


def test_default_fraction_is_point_99_not_point_98():
    """Omitting fraction uses 0.99, reproducing the legacy 0.01/0.99 split."""
    sig = _two_body_signature()
    xs = interactions.DummyCrossSection()
    mc_default = Simulation._build_phase_space_for_signature(
        _fid(), sig, xs, PT.NuLight)
    assert list(mc_default.weights) == pytest.approx([0.01, 0.99], abs=1e-12)

    # A 0.98 fraction is a different split -- proving the default is not 0.98.
    mc_98 = Simulation._build_phase_space_for_signature(
        _fid(), sig, xs, PT.NuLight, fraction=0.98)
    assert list(mc_98.weights) == pytest.approx([0.02, 0.98], abs=1e-12)


def test_default_matches_hand_built_legacy_mixture():
    """The facade default equals a hand-built [physical@0.01, directed@0.99]."""
    sig = _two_body_signature()
    xs = interactions.DummyCrossSection()
    fid = _fid()

    facade = Simulation._build_phase_space_for_signature(
        fid, sig, xs, PT.NuLight, fraction=0.99)

    # Legacy-style hand build: physical appended first, directed second.
    legacy = siren.injection.MultiChannelPhaseSpace()
    legacy.channels = [
        siren.injection.PhysicalCrossSectionChannel(xs, sig),
        siren.injection.DetectorDirected2BodyChannel(fid, 0),
    ]
    legacy.weights = [0.01, 0.99]

    assert len(facade.channels) == len(legacy.channels)
    for a, b in zip(facade.channels, legacy.channels):
        assert a.Name() == b.Name()
    assert list(facade.weights) == pytest.approx(list(legacy.weights), abs=1e-12)


class TestFacadeRunDeterminism:
    """A fixed-seed facade run reproduces its own event stream exactly."""

    def _sim(self):
        return siren.Simulation(
            events=8,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
            process="CC",
            seed=20240607,
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()),
        )

    def test_fixed_seed_stream_identity(self):
        r1 = self._sim().run()
        r2 = self._sim().run()

        assert len(r1) == len(r2)
        for (e1, w1), (e2, w2) in zip(r1, r2):
            rec1 = e1.tree[0].record
            rec2 = e2.tree[0].record
            p1 = list(rec1.primary_momentum)
            p2 = list(rec2.primary_momentum)
            assert p1 == pytest.approx(p2, rel=0, abs=0)
            v1 = list(rec1.interaction_vertex)
            v2 = list(rec2.interaction_vertex)
            assert v1 == pytest.approx(v2, rel=0, abs=0)
            # The DIS final state is what the interaction actually samples;
            # compare it too so a divergence in secondary kinematics (not just
            # the primary/vertex/weight) is caught. Stream identity is bit-exact,
            # so a direct == on the nested four-momenta is the right comparison.
            sm1 = [list(p) for p in rec1.secondary_momenta]
            sm2 = [list(p) for p in rec2.secondary_momenta]
            assert sm1 == sm2
            assert list(rec1.secondary_masses) == list(rec2.secondary_masses)
            assert w1 == pytest.approx(w2, rel=0, abs=0)


def test_bias_targets_default_fraction_is_point_99():
    """The public bias_targets= path defaults the directed split to 0.99.

    A user omitting fraction on bias_targets must get the legacy 0.01/0.99
    split, not the Directed class default (0.98, intentionally different).
    """
    import inspect
    sig = inspect.signature(Simulation._resolve_bias_targets)
    assert sig.parameters["fraction"].default == 0.99

    # Drive the real path: a secondary cross section whose signature is biased
    # toward a target with fraction omitted -> [physical@0.01, directed@0.99].
    xs = interactions.DummyCrossSection()
    daughter = None
    biased_sig = None
    for s in xs.GetPossibleSignatures():
        if s.primary_type == PT.NuMu:
            biased_sig = s
            daughter = s.secondary_types[0]
            break
    assert biased_sig is not None

    sim = object.__new__(Simulation)
    sim._secondary_phase_spaces = {}
    sim._secondary_processes = {PT.NuMu: [xs]}
    sim._resolve_bias_targets(_fid(), daughter, None)

    assert biased_sig in sim._secondary_phase_spaces
    mc = sim._secondary_phase_spaces[biased_sig]
    assert list(mc.weights) == pytest.approx([0.01, 0.99], abs=1e-12)


class TestFacadeForwarding:
    """The facade exposes the Layer-2 objects and forwards to them."""

    def _sim(self):
        return siren.Simulation(
            events=4,
            detector="IceCube",
            primary="NuMu",
            interactions="CSMSDISSplines",
            target="Nucleon",
            process="CC",
            seed=4242,
            energy=siren.dist.PowerLaw(2, 1e3, 1e6),
            direction=siren.dist.IsotropicDirection(),
            position=siren.dist.ColumnDepth(
                600, 600.0, siren.distributions.LeptonDepthFunction()),
        )

    def test_injector_and_weighter_are_layer2_objects(self):
        sim = self._sim()
        assert isinstance(sim.injector, siren.injection.Injector)
        assert isinstance(sim.weighter, siren.injection.Weighter)

    def test_describe_returns_string(self):
        sim = self._sim()
        text = sim.describe()
        assert isinstance(text, str)
        assert "NuMu" in text

    def test_run_optimize_plan_completes(self):
        from siren.tune import Plan
        sim = self._sim()
        results = sim.run(optimize=Plan(rounds=2, events=20))
        assert len(results) == 4

    def test_results_explain_returns_breakdown(self):
        from siren.report import WeightBreakdown
        results = self._sim().run()
        breakdown = results.explain(0)
        assert isinstance(breakdown, WeightBreakdown)
        assert hasattr(breakdown, "total")
        assert hasattr(breakdown, "vertices")

    def test_two_real_same_config_runs_merge(self):
        """Two independently built same-config Simulations pool without error.

        A fresh DetectorModel object per run means config identity must key on
        stable values, not object identity.
        """
        r1 = self._sim().run()
        r2 = self._sim().run()
        merged = siren.Results.merge([r1, r2])
        assert len(merged) == len(r1) + len(r2)
        assert merged.injected == r1.injected + r2.injected
