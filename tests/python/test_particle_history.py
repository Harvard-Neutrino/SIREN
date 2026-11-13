"""Per-particle status and parent links inside a vertex.

A generator's event has history: a struck nucleon, an intermediate
resonance, and the final particles it decayed into. Each secondary now has a
HepMC status (1 final; anything else is intermediate, in the generator's
convention) and an optional parent, the index of the secondary it came from.

The model below wraps DummyCrossSection (nu + target -> nu + Nucleon) and
reports the cascade

    0  NuMu      final
    1  Nucleon   status 21 (struck, before the cascade)
    2  NuMu      status 2  (intermediate), parent 1
    3  NuMu      final, parent 2
    4  Nucleon   final, parent 2

NuMu has a registered secondary process, so the test can check that the two
final neutrinos start secondary vertices and the intermediate one does not.
Statuses, parents and the attachment of daughter vertices must survive
pickle and HepMC3, where the history is written as internal HepMC vertices.
"""
import os
import pickle
import re

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions as dist
from siren import injection as inj
from siren import interactions as xs
from siren import utilities
from siren.math import Vector3D

NuMu = dc.Particle.ParticleType.NuMu
STATUSES = [1, 21, 2, 1, 1]
PARENTS = [-1, -1, 1, 2, 2]


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


class Cascade(xs.CrossSection):
    def __init__(self):
        super().__init__()
        self._d = xs.DummyCrossSection()

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
        lepton, recoil = (r.get_secondary_particle_record(i) for i in (0, 1))
        p_lep, p_rec = list(lepton.four_momentum), np.array(recoil.four_momentum)
        del lepton, recoil
        nucleon = r.secondary_types[1]
        r.secondary_types = [NuMu, nucleon, NuMu, NuMu, nucleon]
        # The intermediate neutrino carries a tenth of the recoil and splits
        # it evenly; the final nucleon keeps the rest.
        direction = p_rec[1:] / np.linalg.norm(p_rec[1:])
        e_mid = 0.1 * p_rec[0]
        p_mid = np.concatenate([[e_mid], e_mid * direction])
        momenta = [p_lep, list(p_rec), list(p_mid), list(0.5 * p_mid), list(p_rec - 0.5 * p_mid)]
        for i, p in enumerate(momenta):
            s = r.get_secondary_particle_record(i)
            s.four_momentum = [float(x) for x in p]
            s.mass = 0.0
            s.helicity = 0.0
        for i, (status, parent) in enumerate(zip(STATUSES, PARENTS)):
            if status != 1:
                r.set_secondary_status(i, status)
            if parent != -1:
                r.set_secondary_parent(i, parent)

    def TotalCrossSection(self, r):
        return self._d.TotalCrossSection(r)

    def TotalCrossSectionAllFinalStates(self, r):
        return self._d.TotalCrossSection(r)

    def DifferentialCrossSection(self, r):
        return self._d.TotalCrossSection(r)

    def FinalStateProbability(self, r):
        return 1.0


@pytest.fixture(scope="module")
def generated():
    dm = _detector_model()
    model = Cascade()
    col = xs.InteractionCollection(NuMu, [model])
    primary = inj.PrimaryInjectionProcess()
    primary.primary_type = NuMu
    primary.interactions = col
    primary.distributions = [
        dist.PrimaryMass(0), dist.PowerLaw(2.0, 0.5, 5.0), dist.IsotropicDirection(),
        dist.PointSourcePositionDistribution(Vector3D(0, 0, 0), 25.0)]
    physical = inj.PhysicalProcess()
    physical.primary_type = NuMu
    physical.interactions = col
    physical.distributions = [dist.PrimaryMass(0), dist.IsotropicDirection()]
    sec_col = xs.InteractionCollection(NuMu, [xs.DummyCrossSection()])
    sec = inj.SecondaryInjectionProcess()
    sec.primary_type = NuMu
    sec.interactions = sec_col
    sec.distributions = [dist.SecondaryPhysicalVertexDistribution()]
    sec_phys = inj.PhysicalProcess()
    sec_phys.primary_type = NuMu
    sec_phys.interactions = sec_col
    sec_phys.distributions = [dist.SecondaryPhysicalVertexDistribution()]
    generator = inj._Injector(40, dm, primary, [sec], utilities.SIREN_random(31))
    generator.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
    weighter = inj._Weighter([generator], dm, physical, [sec_phys])
    events = [generator.GenerateEvent() for _ in range(40)]
    events = [e for e in events if len(e.tree)]
    assert len(events) > 20
    return dict(events=events, weighter=weighter, keepalive=(model, generator, dm))


def _root(event):
    return next(d for d in event.tree if d.is_root())


def test_history_is_stored(generated):
    for e in generated["events"]:
        r = _root(e).record
        assert list(r.secondary_statuses) == STATUSES
        assert list(r.secondary_parents) == PARENTS
        assert [r.is_final_secondary(i) for i in range(5)] == [s == 1 for s in STATUSES]


def test_only_final_particles_start_secondary_vertices(generated):
    for e in generated["events"]:
        root = _root(e).record
        ids = [repr(i) for i in root.secondary_ids]
        children = [repr(d.record.primary_id) for d in e.tree if not d.is_root()]
        assert sorted(children) == sorted([ids[0], ids[3]])


def _state(e):
    out = []
    for d in e.tree:
        r = d.record
        out.append((d.is_root(), [int(t) for t in r.secondary_types], list(r.secondary_statuses),
                    list(r.secondary_parents), repr(r.primary_id)))
    return sorted(out)


def test_pickle_round_trip(generated):
    for e in generated["events"]:
        assert _state(pickle.loads(pickle.dumps(e))) == _state(e)


def test_hepmc3_round_trip(generated, tmp_path):
    events, weighter = generated["events"], generated["weighter"]
    from siren import hepmc3 as sio

    out = str(tmp_path / "cascade.hepmc3")
    try:
        sio.SaveInteractionTreesAsHepMC3(events, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise
    back = sio.LoadInteractionTreesFromHepMC3(out)
    assert len(back) == len(events)
    for a, b in zip(events, back):
        assert _state(b) == _state(a)
        assert weighter.EventWeight(b) == pytest.approx(weighter.EventWeight(a), rel=1e-9)
    text = open(out).read()
    # HepMC particle lines: "P id production_vertex pid px py pz e m status".
    statuses = {int(line.split()[-1]) for line in text.splitlines() if line.startswith("P ")}
    assert {21, 2}.issubset(statuses)
    assert "siren.intra_vertex" in text
