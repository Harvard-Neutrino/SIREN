"""Integer and string vertex attributes, and generator process IDs.

A generator reports more than numbers on a continuous scale: a process ID, a
channel name, a count. A sampler can now store these on the vertex as integer
and string attributes next to the double interaction_parameters. They survive
pickle, the native archive and HepMC3, where they appear as siren.iparam.*
and siren.sparam.* vertex attributes. A root vertex that carries the reserved
integer "process_id" (and optionally "process_name") is registered in the
NuHepMC process list under that ID and gives the event's signal_process_id;
other events keep SIREN's signature-keyed IDs from 700 upward.
"""
import os
import pickle
import re

import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import distributions as dist
from siren import injection as inj
from siren import interactions as xs
from siren import utilities
from siren.math import Vector3D

NuMu = dc.Particle.ParticleType.NuMu
PROCESSES = {1001: "toy QE", 1002: "toy RES"}


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


class Tagged(xs.CrossSection):
    """DummyCrossSection that tags each vertex like a generator would.

    Events below 2 GeV are left without a process ID, so the file mixes
    generator-named and signature-keyed processes.
    """

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
        energy = r.primary_momentum[0]
        ints = {"n_nucleons": 1 + int(3 * rand.Uniform(0, 1)), "large": 2 ** 40 + 1}
        strings = {"channel": "Delta(1232) -> p pi+" if energy > 3 else "QE"}
        if energy >= 2:
            pid = 1002 if energy > 3 else 1001
            ints["process_id"] = pid
            strings["process_name"] = PROCESSES[pid]
        r.interaction_int_parameters = ints
        r.interaction_string_parameters = strings

    def TotalCrossSection(self, r):
        return self._d.TotalCrossSection(r)

    def TotalCrossSectionAllFinalStates(self, r):
        return self._d.TotalCrossSection(r)

    def DifferentialCrossSection(self, r):
        return self._d.DifferentialCrossSection(r)

    def FinalStateProbability(self, r):
        return self._d.FinalStateProbability(r)


@pytest.fixture(scope="module")
def events():
    dm = _detector_model()
    model = Tagged()
    process = inj.PrimaryInjectionProcess()
    process.primary_type = NuMu
    process.interactions = xs.InteractionCollection(NuMu, [model])
    process.distributions = [
        dist.PrimaryMass(0), dist.PowerLaw(1.0, 0.5, 5.0), dist.IsotropicDirection(),
        dist.PointSourcePositionDistribution(Vector3D(0, 0, 0), 25.0)]
    generator = inj._Injector(200, dm, process, utilities.SIREN_random(11))
    out = [generator.GenerateEvent() for _ in range(200)]
    out = [e for e in out if len(e.tree)]
    assert len(out) > 100
    return out, (model, generator, dm)


def _attrs(event):
    r = event.tree[0].record
    return dict(r.interaction_int_parameters), dict(r.interaction_string_parameters)


def test_attributes_are_stored(events):
    evs, _ = events
    with_pid = [e for e in evs if "process_id" in _attrs(e)[0]]
    assert 0 < len(with_pid) < len(evs)
    for e in evs:
        ints, strings = _attrs(e)
        assert ints["large"] == 2 ** 40 + 1
        assert 1 <= ints["n_nucleons"] <= 3
        assert strings["channel"] in ("QE", "Delta(1232) -> p pi+")
        if "process_id" in ints:
            assert strings["process_name"] == PROCESSES[ints["process_id"]]


def test_pickle_and_native_archive(events, tmp_path):
    evs, _ = events
    assert [_attrs(pickle.loads(pickle.dumps(e))) for e in evs] == [_attrs(e) for e in evs]
    path = str(tmp_path / "tagged.siren_events")
    dc.SaveInteractionTrees(evs, path)
    assert [_attrs(e) for e in dc.LoadInteractionTrees(path)] == [_attrs(e) for e in evs]


def test_hepmc3_keeps_attributes_and_registers_generator_processes(events, tmp_path):
    evs, _ = events
    from siren import hepmc3 as sio

    out = str(tmp_path / "tagged.hepmc3")
    try:
        sio.SaveInteractionTreesAsHepMC3(evs, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise
    back = sio.LoadInteractionTreesFromHepMC3(out)
    assert [_attrs(e) for e in back] == [_attrs(e) for e in evs]

    text = open(out).read()
    registry = re.search(r"NuHepMC\.ProcessIDs\s+(.*)", text)
    assert registry is not None
    ids = [int(x) for x in re.findall(r"-?\d+", registry.group(1))]
    assert 1001 in ids and 1002 in ids
    assert any(i >= 700 and i not in PROCESSES for i in ids)
    assert re.search(r"NuHepMC\.ProcessInfo\[1001\]\.Name\s+toy QE", text)
    signal = [int(x) for x in re.findall(r"signal_process_id\s+(\d+)", text)]
    assert len(signal) == len(evs)
    for e, sid in zip(evs, signal):
        pid = _attrs(e)[0].get("process_id")
        if pid is not None:
            assert sid == pid
        else:
            assert sid >= 700 and sid not in PROCESSES
