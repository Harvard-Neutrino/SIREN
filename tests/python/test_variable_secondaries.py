"""A sampler may set the number and types of a vertex's secondaries.

The signature a model declares is the channel: the injector selects it and the
weighter evaluates it. What the model realises is stored separately on the
record. This test wraps the data-free DummyCrossSection (nu + target -> nu +
Nucleon) in a Python model that adds zero to three extra muon neutrinos to the
realised list of each event, borrowing their energy from the recoil system.

It checks that

* the signature stays the declared channel while the realised list varies;
* every realised muon neutrino, and only those, starts a secondary vertex;
* the weight is consistent: each daughter vertex contributes its own factor
  (the probability that its forced interaction happens), independent of how
  many siblings it has, so the full tree's weight factorises over daughters;
  and the root weighs the same as the same root with only the channel's two
  secondaries;
* the realised list and the weight survive pickle, the native archive and a
  HepMC3 round trip.
"""
import os
import pickle

import numpy as np
import pytest

siren = pytest.importorskip("siren")

N_INJECT = 60
MAX_EXTRA = 3
EXTRA_ENERGY_FRACTION = 0.05


def _detector_model():
    from siren import detector
    from siren import _util

    try:
        dm = detector.DetectorModel()
        det_dir = _util.get_detector_model_path("CCM")
        dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
        dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    except (ValueError, RuntimeError) as e:
        pytest.skip(f"CCM detector model unavailable: {e}")
    return dm


def _variable_model(interactions, dc):
    NuMu = dc.Particle.ParticleType.NuMu

    class ExtraNeutrinos(interactions.CrossSection):
        """DummyCrossSection with 0..MAX_EXTRA extra muon neutrinos."""

        def __init__(self):
            super().__init__()
            self._d = interactions.DummyCrossSection()

        def equal(self, other):
            return other is self

        def GetPossiblePrimaries(self):
            return self._d.GetPossiblePrimaries()

        def GetPossibleTargets(self):
            return self._d.GetPossibleTargets()

        def GetPossibleTargetsFromPrimary(self, p):
            return self._d.GetPossibleTargetsFromPrimary(p)

        def GetPossibleSignatures(self):
            return self._d.GetPossibleSignatures()

        def GetPossibleSignaturesFromParents(self, p, t):
            return self._d.GetPossibleSignaturesFromParents(p, t)

        def InteractionThreshold(self, r):
            return self._d.InteractionThreshold(r)

        def DensityVariables(self):
            return self._d.DensityVariables()

        def SampleFinalState(self, r, rand):
            self._d.SampleFinalState(r, rand)
            channel = list(r.secondary_types)
            assert channel == list(r.signature.secondary_types) and channel[0] == NuMu
            lepton = r.get_secondary_particle_record(0)
            hadrons = r.get_secondary_particle_record(1)
            p_lep, h_lep = list(lepton.four_momentum), lepton.helicity
            p_had, h_had = np.array(hadrons.four_momentum), hadrons.helicity
            # Replacing the list rebuilds the particle records, so handles
            # obtained before it must not be used afterwards (as with
            # references into a std::vector that is resized).
            del lepton, hadrons
            n_extra = min(int(rand.Uniform(0, MAX_EXTRA + 1)), MAX_EXTRA)
            if n_extra == 0:
                return
            r.secondary_types = channel + [NuMu] * n_extra
            direction = p_had[1:] / np.linalg.norm(p_had[1:])
            e_extra = EXTRA_ENERGY_FRACTION * p_had[0]
            p_extra = np.concatenate([[e_extra], e_extra * direction])
            p_rest = p_had - n_extra * p_extra
            m_rest = np.sqrt(p_rest[0] ** 2 - np.dot(p_rest[1:], p_rest[1:]))
            new_lepton = r.get_secondary_particle_record(0)
            new_lepton.four_momentum = p_lep
            new_lepton.mass = 0.0
            new_lepton.helicity = h_lep
            rest = r.get_secondary_particle_record(1)
            rest.four_momentum = list(p_rest)
            rest.mass = float(m_rest)
            rest.helicity = h_had
            for i in range(2, 2 + n_extra):
                extra = r.get_secondary_particle_record(i)
                extra.four_momentum = list(p_extra)
                extra.mass = 0.0
                extra.helicity = -1.0

        def TotalCrossSection(self, r):
            return self._d.TotalCrossSection(r)

        def TotalCrossSectionAllFinalStates(self, r):
            return self._d.TotalCrossSection(r)

        def DifferentialCrossSection(self, r):
            return self._d.TotalCrossSection(r)

        def FinalStateProbability(self, r):
            return 1.0

    return ExtraNeutrinos()


@pytest.fixture(scope="module")
def generated():
    from siren import dataclasses as dc
    from siren import injection
    from siren import interactions
    from siren import distributions
    from siren import math as smath
    from siren import utilities

    NuMu = dc.Particle.ParticleType.NuMu
    dm = _detector_model()

    xs = _variable_model(interactions, dc)
    int_col = interactions.InteractionCollection(NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]
    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    sec_col = interactions.InteractionCollection(NuMu, [interactions.DummyCrossSection()])
    sec_inj = injection.SecondaryInjectionProcess()
    sec_inj.primary_type = NuMu
    sec_inj.interactions = sec_col
    sec_inj.distributions = [distributions.SecondaryPhysicalVertexDistribution()]
    sec_phys = injection.PhysicalProcess()
    sec_phys.primary_type = NuMu
    sec_phys.interactions = sec_col
    sec_phys.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    rand = utilities.SIREN_random(20261112)
    inj = injection._Injector(N_INJECT, dm, primary_inj, [sec_inj], rand)
    # One generation of secondaries: expand the root's particles only.
    inj.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
    weighter = injection._Weighter([inj], dm, primary_phys, [sec_phys])

    events = []
    for _ in range(N_INJECT):
        try:
            ev = inj.GenerateEvent()
        except RuntimeError:
            break
        if len(ev.tree) > 0:
            events.append(ev)
    if len(events) < N_INJECT // 2:
        pytest.skip("too few non-empty events generated")
    return dict(events=events, weighter=weighter, keepalive=(xs, inj, dm), dc=dc)


def _root(tree):
    roots = [d for d in tree.tree if d.is_root()]
    assert len(roots) == 1
    return roots[0]


def _subtree(dc, root, children=()):
    tree = dc.InteractionTree()
    parent = tree.add_entry(dc.InteractionRecord(root), None)
    for child in children:
        tree.add_entry(dc.InteractionRecord(child), parent)
    return tree


def _root_only(dc, record):
    return _subtree(dc, record)


def _fields(r):
    return (r.signature.primary_type, r.signature.target_type,
            list(r.signature.secondary_types), list(r.secondary_types),
            r.has_realised_secondary_types, repr(r.primary_id),
            list(r.primary_initial_position), r.primary_initial_time, r.primary_mass,
            list(r.primary_momentum), r.primary_helicity, repr(r.target_id),
            r.target_mass, r.target_helicity, list(r.interaction_vertex),
            r.interaction_time, [repr(i) for i in r.secondary_ids],
            list(r.secondary_masses), [list(m) for m in r.secondary_momenta],
            list(r.secondary_helicities), list(r.secondary_times),
            dict(r.interaction_parameters))


def test_signature_is_the_channel_and_the_realised_list_varies(generated):
    dc = generated["dc"]
    NuMu = dc.Particle.ParticleType.NuMu
    lengths = set()
    for ev in generated["events"]:
        rec = _root(ev).record
        channel = list(rec.signature.secondary_types)
        assert len(channel) == 2 and channel[0] == NuMu
        realised = list(rec.secondary_types)
        assert realised[:2] == channel
        assert all(t == NuMu for t in realised[2:])
        assert rec.has_realised_secondary_types == (len(realised) > 2)
        n = len(realised)
        for field in (rec.secondary_ids, rec.secondary_masses, rec.secondary_momenta,
                      rec.secondary_helicities, rec.secondary_times):
            assert len(field) == n
        lengths.add(n)
    # Every multiplicity from the channel's two to the maximum occurs.
    assert lengths == set(range(2, 3 + MAX_EXTRA))


def test_every_realised_neutrino_feeds_the_secondary_process(generated):
    dc = generated["dc"]
    NuMu = dc.Particle.ParticleType.NuMu
    for ev in generated["events"]:
        root = _root(ev)
        realised = list(root.record.secondary_types)
        ids = [repr(i) for i in root.record.secondary_ids]
        numu_ids = {ids[i] for i, t in enumerate(realised) if t == NuMu}
        children = [d for d in ev.tree if not d.is_root()]
        assert len(children) == len(numu_ids)
        assert {repr(d.record.primary_id) for d in children} == numu_ids
        for d in children:
            assert d.record.signature.primary_type == NuMu
            i = ids.index(repr(d.record.primary_id))
            assert list(d.record.primary_momentum) == list(root.record.secondary_momenta[i])


def test_weight_is_independent_of_the_realised_final_state(generated):
    dc = generated["dc"]
    weighter = generated["weighter"]
    for ev in generated["events"]:
        w = weighter.EventWeight(ev)
        assert np.isfinite(w) and w > 0
        root = _root(ev).record
        children = [d.record for d in ev.tree if not d.is_root()]
        w_root = weighter.EventWeight(_root_only(dc, root))
        assert np.isfinite(w_root) and w_root > 0
        # Each daughter contributes a factor of its own, whatever its siblings.
        factors = [weighter.EventWeight(_subtree(dc, root, [c])) / w_root for c in children]
        assert all(np.isfinite(f) and 0 < f < 1 for f in factors)
        assert w == pytest.approx(w_root * np.prod(factors), rel=1e-9)
        # The root with the channel's two secondaries weighs the same.
        fixed = dc.InteractionRecord(root)
        fixed.secondary_types = list(root.signature.secondary_types)
        assert not fixed.has_realised_secondary_types
        for name in ("secondary_ids", "secondary_masses", "secondary_momenta",
                     "secondary_helicities", "secondary_times"):
            setattr(fixed, name, list(getattr(root, name))[:2])
        assert weighter.EventWeight(_root_only(dc, fixed)) == pytest.approx(w_root, rel=1e-12)


def _assert_same_events(before, after, weighter, rel):
    assert len(after) == len(before)
    for a, b in zip(before, after):
        assert len(a.tree) == len(b.tree)
        ra, rb = _root(a).record, _root(b).record
        assert list(rb.secondary_types) == list(ra.secondary_types)
        assert rb.signature == ra.signature
        assert rb.has_realised_secondary_types == ra.has_realised_secondary_types
        assert weighter.EventWeight(b) == pytest.approx(weighter.EventWeight(a), rel=rel)


def test_pickle_round_trip(generated):
    events, weighter = generated["events"], generated["weighter"]
    restored = [pickle.loads(pickle.dumps(ev)) for ev in events]
    _assert_same_events(events, restored, weighter, rel=0)
    for a, b in zip(events, restored):
        for da, db in zip(a.tree, b.tree):
            # repr is exact for floats and lets NaN (DummyCrossSection gives some
            # recoil masses as NaN) compare equal to itself.
            assert repr(_fields(da.record)) == repr(_fields(db.record))


def test_native_archive_round_trip(generated, tmp_path):
    dc = generated["dc"]
    events, weighter = generated["events"], generated["weighter"]
    path = str(tmp_path / "events.siren_events")
    dc.SaveInteractionTrees(events, path)
    restored = dc.LoadInteractionTrees(path)
    _assert_same_events(events, restored, weighter, rel=0)


def test_hepmc3_round_trip(generated, tmp_path):
    events, weighter = generated["events"], generated["weighter"]
    from siren import hepmc3 as sio

    out = str(tmp_path / "variable.hepmc3")
    try:
        sio.SaveInteractionTreesAsHepMC3(events, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise
    restored = sio.LoadInteractionTreesFromHepMC3(out)
    _assert_same_events(events, restored, weighter, rel=1e-9)
