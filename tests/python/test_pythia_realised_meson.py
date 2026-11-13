"""PythiaDISCrossSection under the realised-secondaries and opaque-final-state rules.

The model registers one signature per charm species (D0, D+, Ds+ for nu_mu).
The injector selects one by rate; Pythia then produces whichever charm meson
its fragmentation gives. That meson is now recorded as the realised type of
the meson slot (CrossSectionDistributionRecord::SetSecondaryType), while the
signature stays the selected channel. Before, the selected signature itself
was overwritten through a const_cast.

The checks:

* the realised species is Pythia's own: its distribution does not depend on
  which channel was selected, and the signature stays the selected channel;
* the weight is the same whether the signature names the selected species or
  the realised one (the old behaviour), because the species' fragmentation
  fraction enters the rate and the selection probability alike;
* two Pythia objects configured alike count as one final-state sampler, and a
  differently configured one is refused (the final state is opaque without a
  differential spline);
* HepMC3 keeps the realised species and the channel.

Needs SIREN built with Pythia8 and:
  SIREN_PYTHIA_WIDE_SIGMA  total sigma(E) FITS spline (CC nu_mu on p)
  LHAPDF_DATA_PATH         with the PDF set (default LHAPDF6:CT18NLO)
  PYTHIA8DATA              Pythia's xmldoc directory
On macOS Pythia also needs its lib directory on DYLD_LIBRARY_PATH, where it
finds the LHAPDF6 plugin.
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

P = dc.Particle.ParticleType
SIGMA = os.environ.get("SIREN_PYTHIA_WIDE_SIGMA")
PDF = os.environ.get("SIREN_PYTHIA_PDF", "LHAPDF6:CT18NLO")
PYTHIA_DATA = os.environ.get("PYTHIA8DATA", "")
M_N = 0.9314943
ENERGY = 100.0
D_TYPES = (P.D0, P.DPlus, P.DsPlus)

pytestmark = pytest.mark.skipif(
    not hasattr(xs, "PythiaDISCrossSection") or not SIGMA or not os.path.exists(SIGMA or "")
    or not os.environ.get("LHAPDF_DATA_PATH") or not PYTHIA_DATA,
    reason="needs SIREN with Pythia8 and SIREN_PYTHIA_WIDE_SIGMA, LHAPDF_DATA_PATH, PYTHIA8DATA")


def pythia(minimum_Q2=1.0):
    return xs.PythiaDISCrossSection("", SIGMA, 1, M_N, minimum_Q2, [P.NuMu], [P.PPlus],
                                    PYTHIA_DATA, PDF, "cm")


@pytest.fixture(scope="module")
def model():
    return pythia()


def channels(model):
    sigs = list(model.GetPossibleSignaturesFromParents(P.NuMu, P.PPlus))
    assert sorted(int(s.secondary_types[2]) for s in sigs) == sorted(int(t) for t in D_TYPES)
    return {s.secondary_types[2]: s for s in sigs}


def sample(model, signature, n, seed):
    rng = utilities.SIREN_random(seed)
    out = []
    while len(out) < n:
        ir = dc.InteractionRecord()
        ir.signature = signature
        ir.primary_momentum = [ENERGY, 0.0, 0.0, ENERGY]
        ir.target_mass = M_N
        cdr = dc.CrossSectionDistributionRecord(ir)
        model.SampleFinalState(cdr, rng)
        done = dc.InteractionRecord()
        done.signature = ir.signature
        done.primary_momentum = ir.primary_momentum
        cdr.finalize(done)
        out.append(done)
    return out


def test_the_realised_species_is_pythias_and_the_signature_the_channel(model):
    table = []
    for selected, sig in channels(model).items():
        records = sample(model, sig, 40, 100 + int(selected) % 97)
        row = []
        for r in records:
            assert r.signature == sig
            assert r.signature.secondary_types[2] == selected
            assert r.secondary_types[2] in D_TYPES
        for t in D_TYPES:
            row.append(sum(r.secondary_types[2] == t for r in records))
        table.append(row)
    table = np.array(table, float)
    # Independence of the realised species from the selected channel.
    expected = table.sum(1, keepdims=True) * table.sum(0, keepdims=True) / table.sum()
    keep = expected > 0
    chi2 = (((table - expected) ** 2)[keep] / expected[keep]).sum()
    assert chi2 < 18.5, (chi2, table)  # p ~ 0.001 for 4 degrees of freedom
    # Pythia's species differ from the selected one in a good share of events.
    off_diagonal = table.sum() - np.trace(table)
    assert off_diagonal > 0.2 * table.sum()


def argon_world():
    dm = siren.detector.DetectorModel()
    materials = dm.Materials
    materials.AddMaterial("ARGON", {1000180400: 1.0})
    dm.Materials = materials
    sector = siren.detector.DetectorSector()
    sector.name = "argon"
    sector.level = 0
    sector.material_id = materials.GetMaterialId("ARGON")
    sector.geo = siren.geometry.Sphere(siren.geometry.Placement(siren.math.Vector3D(0, 0, 0)), 100.0, 0.0)
    sector.density = siren.detector.ConstantDensityDistribution(1.4)
    dm.AddSector(sector)
    return dm


def processes(gen_model, phys_model):
    from siren import distributions as dist
    from siren import geometry as geo
    from siren.math import Vector3D
    primary = inj.PrimaryInjectionProcess()
    primary.primary_type = P.NuMu
    primary.interactions = xs.InteractionCollection(P.NuMu, [gen_model])
    primary.distributions = [
        dist.PrimaryMass(0), dist.Monoenergetic(ENERGY), dist.FixedDirection(Vector3D(0, 0, 1)),
        dist.CylinderVolumePositionDistribution(geo.Cylinder(geo.Placement(Vector3D(0, 0, 0)), 1.0, 0.0, 2.0)),
    ]
    physical = inj.PhysicalProcess()
    physical.primary_type = P.NuMu
    physical.interactions = xs.InteractionCollection(P.NuMu, [phys_model])
    physical.distributions = [dist.PrimaryMass(0), dist.Monoenergetic(ENERGY),
                              dist.FixedDirection(Vector3D(0, 0, 1))]
    return primary, physical


@pytest.fixture(scope="module")
def injected(model):
    dm = argon_world()
    primary, physical = processes(model, model)
    generator = inj._Injector(12, dm, primary, utilities.SIREN_random(5))
    events = [generator.GenerateEvent() for _ in range(12)]
    assert all(len(e.tree) for e in events)
    return dict(dm=dm, generator=generator, physical=physical, events=events,
                weighter=inj._Weighter([generator], dm, physical))


def test_the_weight_does_not_depend_on_which_species_the_signature_names(model, injected):
    by_species = channels(model)
    weighter = injected["weighter"]
    for e in injected["events"]:
        record = e.tree[0].record
        w = weighter.EventWeight(e)
        assert np.isfinite(w) and w > 0
        # The old behaviour: the signature names the realised species.
        old = dc.InteractionRecord(record)
        old.signature = by_species[record.secondary_types[2]]
        old.secondary_types = list(old.signature.secondary_types)
        tree = dc.InteractionTree()
        tree.add_entry(old, None)
        assert weighter.EventWeight(tree) == pytest.approx(w, rel=1e-12)


def test_alike_pythia_objects_share_a_sampler_and_others_are_refused(model, injected):
    from siren.errors import ConfigurationError
    assert not model.FinalStateDensityEvaluable(channels(model)[P.D0])
    key = model.FinalStateSamplerKey(channels(model)[P.D0])
    assert "Pythia8" in key and PDF in key
    twin = pythia()
    assert twin.FinalStateSamplerKey(channels(twin)[P.D0]) == key
    _, physical = processes(model, twin)
    w = inj._Weighter([injected["generator"]], injected["dm"], physical)
    e = injected["events"][0]
    assert w.EventWeight(e) == pytest.approx(injected["weighter"].EventWeight(e), rel=1e-12)
    other = pythia(minimum_Q2=2.0)
    _, physical = processes(model, other)
    with pytest.raises(ConfigurationError, match="sampler"):
        inj._Weighter([injected["generator"]], injected["dm"], physical)


def test_hepmc3_keeps_the_realised_species_and_the_channel(injected, tmp_path):
    from siren import hepmc3 as sio
    out = str(tmp_path / "pythia.hepmc3")
    sio.SaveInteractionTreesAsHepMC3(injected["events"], out)
    back = sio.LoadInteractionTreesFromHepMC3(out)
    for a, b in zip(injected["events"], back):
        ra, rb = a.tree[0].record, b.tree[0].record
        assert rb.signature == ra.signature
        assert list(rb.secondary_types) == list(ra.secondary_types)
        assert injected["weighter"].EventWeight(b) == pytest.approx(injected["weighter"].EventWeight(a), rel=1e-9)


@pytest.mark.skipif(not os.environ.get("SIREN_PYTHIA_WIDE_DSDXDY"),
                    reason="needs SIREN_PYTHIA_WIDE_DSDXDY (differential spline)")
def test_a_differential_spline_makes_the_final_state_evaluable():
    with_density = xs.PythiaDISCrossSection(os.environ["SIREN_PYTHIA_WIDE_DSDXDY"], SIGMA, 1, M_N, 1.0,
                                            [P.NuMu], [P.PPlus], PYTHIA_DATA, PDF, "cm")
    for signature in with_density.GetPossibleSignaturesFromParents(P.NuMu, P.PPlus):
        assert with_density.FinalStateDensityEvaluable(signature)
    assert not pythia().FinalStateDensityEvaluable(signature)
