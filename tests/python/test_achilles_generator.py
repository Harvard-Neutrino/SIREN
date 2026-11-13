"""Achilles through AchillesTrialGenerator, BatchResampler and ExternalGeneratorCrossSection.

Needs a SIREN built with SIREN_WITH_ACHILLES=ON and two environment variables:

  SIREN_ACHILLES_RUN_CARD   an Achilles run card with a Geometry beam that
                            configures nu_mu CC on Ar40
  SIREN_ACHILLES_DATA_DIR   the Achilles build or install directory (data/)

The checks:

* vertex content: the outgoing muon first, then every non-spectator hadron
  with Achilles' status, exactly one struck nucleon (status 3), Achilles'
  process ID, and baryon number (spectators + final and captured nucleons
  = 40, or 41 after a pion interacted in the cascade, as Achilles reports it);
* lab frame: for a beam along +x, the muons point along +x on average;
* rate: a monoenergetic 2 GeV beam injected in a cylinder inside uniform
  argon gives sum(w) = sigma * n_Ar * V. The cross section in the oracle is
  measured here from plain Achilles trials (envelope * mean trial weight),
  independently of the batch estimator, the injector and the weighter;
* HepMC3 keeps statuses and Achilles' process IDs;
* invalid channels and run cards are refused.
"""
import math
import os

import numpy as np
import pytest

siren = pytest.importorskip("siren")

from siren import dataclasses as dc
from siren import injection as inj
from siren import interactions as xs
from siren import utilities

import test_external_generator as common

P = dc.Particle.ParticleType
RUN_CARD = os.environ.get("SIREN_ACHILLES_RUN_CARD")
DATA_DIR = os.environ.get("SIREN_ACHILLES_DATA_DIR")

pytestmark = pytest.mark.skipif(
    not getattr(xs, "HAS_ACHILLES", False) or not RUN_CARD or not DATA_DIR,
    reason="needs SIREN built with Achilles and SIREN_ACHILLES_RUN_CARD/SIREN_ACHILLES_DATA_DIR")

ENERGY = 2.0
BATCH = 32
ALLOWED_STATUSES = {1, 3, 26, 29}


def channel(target=P.Ar40Nucleus, lepton=P.MuMinus, secondaries=None):
    s = dc.InteractionSignature()
    s.primary_type = P.NuMu
    s.target_type = target
    s.secondary_types = secondaries if secondaries is not None else [lepton, P.Hadrons]
    return s


@pytest.fixture(scope="module")
def achilles():
    cwd = os.getcwd()
    gen = xs.AchillesTrialGenerator(RUN_CARD, DATA_DIR, [channel()])
    assert os.getcwd() == cwd
    return gen


def processes(model, dm, n, seed, direction=(0.0, 0.0, 1.0)):
    """common.processes with a choice of beam direction."""
    from siren import distributions as dist
    from siren import geometry as geo
    from siren.math import Vector3D
    primary = inj.PrimaryInjectionProcess()
    primary.primary_type = P.NuMu
    primary.interactions = xs.InteractionCollection(P.NuMu, [model])
    primary.distributions = [
        dist.PrimaryMass(0), dist.Monoenergetic(ENERGY), dist.FixedDirection(Vector3D(*direction)),
        dist.CylinderVolumePositionDistribution(
            geo.Cylinder(geo.Placement(Vector3D(0, 0, 0)), common.RADIUS, 0.0, common.LENGTH)),
    ]
    physical = inj.PhysicalProcess()
    physical.primary_type = P.NuMu
    physical.interactions = xs.InteractionCollection(P.NuMu, [model])
    physical.distributions = [dist.PrimaryMass(0), dist.Monoenergetic(ENERGY),
                              dist.FixedDirection(Vector3D(*direction))]
    return inj._Injector(n, dm, primary, utilities.SIREN_random(seed)), physical


def _events(gen, n, seed, direction=(0.0, 0.0, 1.0), batch=BATCH):
    model = xs.ExternalGeneratorCrossSection(xs.BatchResampler(gen, batch))
    dm = common.argon_world()
    generator, physical = processes(model, dm, n, seed, direction)
    cwd = os.getcwd()
    events = [generator.GenerateEvent() for _ in range(n)]
    assert os.getcwd() == cwd
    assert all(len(e.tree) for e in events)
    weighter = inj._Weighter([generator], dm, physical)
    return events, weighter, (model, generator, dm)


def test_the_channel_is_opaque_and_matches_the_generator(achilles):
    lo, hi = achilles.EnergyRange()
    assert lo < ENERGY < hi
    assert achilles.Envelope(channel(), ENERGY) > achilles.TotalCrossSection(channel(), ENERGY) > 0
    assert achilles.TotalCrossSection(channel(), hi * 1.01) == 0.0
    model = xs.ExternalGeneratorCrossSection(xs.BatchResampler(achilles, BATCH))
    assert not model.FinalStateDensityEvaluable(channel())
    assert model.FinalStateSamplerKey(channel()).startswith("Achilles|card=")


def test_vertex_content(achilles):
    events, _, keep = _events(achilles, 200, 3)
    nucleons = {P.PPlus, P.Neutron}
    pions = {P.PiPlus, P.PiMinus, P.Pi0}
    excesses = []
    for e in events:
        r = e.tree[0].record
        if r.sampler_factor == 0.0:
            continue
        types = list(r.secondary_types)
        statuses = [r.get_secondary_status(i) for i in range(len(types))]
        assert types[0] == P.MuMinus and statuses[0] == 1
        assert set(statuses) <= ALLOWED_STATUSES, statuses
        assert statuses.count(3) == 1
        assert types[statuses.index(3)] in nucleons
        assert "process_id" in r.interaction_int_parameters
        # Baryon number: the struck nucleon (3) and nucleons that interacted
        # again in the cascade (29) are earlier states of others, so the
        # nucleus' 40 nucleons are the spectators plus the final (1) and
        # captured (26) nucleons. After a pion interacts in the cascade
        # (absorption) Achilles 9d86fb6c reports one nucleon more, which looks
        # like the absorbing nucleon kept among the spectators; SIREN passes
        # Achilles' list through as it is.
        remaining = sum(t in nucleons and s in (1, 26) for t, s in zip(types, statuses))
        excess = r.interaction_int_parameters["achilles_spectators"] + remaining - 40
        pion_interacted = any(t in pions and s == 29 for t, s in zip(types, statuses))
        assert excess == 0 or (excess == 1 and pion_interacted), (excess, list(zip(types, statuses)))
        excesses.append(excess)
        assert r.interaction_parameters["achilles_weight"] >= 1.0
        assert 0 < r.secondary_momenta[0][0] < ENERGY
        assert all(np.all(np.isfinite(p)) for p in r.secondary_momenta)
        assert r.signature == channel()
    # Most events balance exactly.
    assert excesses.count(0) > 0.7 * len(excesses)


def test_events_are_in_the_lab_frame(achilles):
    events, _, keep = _events(achilles, 150, 4, direction=(1.0, 0.0, 0.0))
    u = []
    for e in events:
        r = e.tree[0].record
        if r.sampler_factor == 0.0:
            continue
        p = np.array(r.secondary_momenta[0][1:])
        u.append(p / np.linalg.norm(p))
    u = np.mean(u, axis=0)
    assert u[0] > 0.5 and abs(u[1]) < 0.1 and abs(u[2]) < 0.1, u


def test_rate_matches_flux_times_cross_section_times_targets(achilles):
    # Oracle cross section from plain trials at the beam energy.
    record = dc.InteractionRecord()
    record.signature = channel()
    record.primary_momentum = [ENERGY, 0.0, 0.0, ENERGY]
    record.target_mass = 37.215
    cdr = dc.CrossSectionDistributionRecord(record)
    rand = utilities.SIREN_random(1)
    trials = 40000
    w = np.array([achilles.Trial(cdr, rand) for _ in range(trials)])
    envelope = achilles.Envelope(channel(), ENERGY)
    sigma = envelope * w.mean()
    sigma_err = envelope * w.std(ddof=1) / math.sqrt(trials)
    spline = achilles.TotalCrossSection(channel(), ENERGY)
    assert abs(sigma / spline - 1) < 0.05

    events, weighter, keep = _events(achilles, 1500, 7)
    weights = np.array([weighter.EventWeight(e) for e in events])
    n_ar = common.DENSITY / 39.9623831 * common.AVOGADRO
    volume = math.pi * (100 * common.RADIUS) ** 2 * (100 * common.LENGTH)
    expected = sigma * n_ar * volume * 1e-4
    error = math.hypot(math.sqrt((weights ** 2).sum()), expected * sigma_err / sigma)
    assert abs(weights.sum() - expected) < 4 * error, (weights.sum(), expected, error)


def test_hepmc3_keeps_statuses_and_process_ids(achilles, tmp_path):
    events, weighter, keep = _events(achilles, 40, 9)
    events = [e for e in events if e.tree[0].record.sampler_factor > 0]
    from siren import hepmc3 as sio
    out = str(tmp_path / "achilles.hepmc3")
    sio.SaveInteractionTreesAsHepMC3(events, out)
    back = sio.LoadInteractionTreesFromHepMC3(out)
    for a, b in zip(events, back):
        ra, rb = a.tree[0].record, b.tree[0].record
        assert list(rb.secondary_types) == list(ra.secondary_types)
        assert list(rb.secondary_statuses) == list(ra.secondary_statuses)
        assert rb.signature == ra.signature
        assert rb.interaction_int_parameters["process_id"] == ra.interaction_int_parameters["process_id"]
        assert weighter.EventWeight(b) == pytest.approx(weighter.EventWeight(a), rel=1e-9)
    text = open(out).read()
    for e in events:
        assert ("signal_process_id %d" % e.tree[0].record.interaction_int_parameters["process_id"]) in text


def test_invalid_configurations_are_refused(achilles, tmp_path):
    with pytest.raises(ValueError, match="lepton \\+ Hadrons"):
        xs.AchillesTrialGenerator(RUN_CARD, DATA_DIR, [channel(secondaries=[P.MuMinus, P.PPlus])])
    with pytest.raises(ValueError, match="no process"):
        xs.AchillesTrialGenerator(RUN_CARD, DATA_DIR, [channel(target=P.O16Nucleus)])
    card = tmp_path / "flat.yml"
    card.write_text("Beam:\n  Type: Monochromatic\n  Energy: 1000\n  Species: [14]\n")
    with pytest.raises(ValueError, match="Geometry"):
        xs.AchillesTrialGenerator(str(card), DATA_DIR, [channel()])
    record = dc.InteractionRecord()
    record.signature = channel(target=P.O16Nucleus)
    with pytest.raises(ValueError, match="not one of"):
        achilles.Trial(dc.CrossSectionDistributionRecord(record), utilities.SIREN_random(2))
