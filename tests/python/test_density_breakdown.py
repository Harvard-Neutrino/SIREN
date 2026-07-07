"""MultiChannelPhaseSpace.DensityBreakdown invariant.

DensityBreakdown(record)[i] is channel i's alpha-weighted, common-measure
contribution; the elements must sum to exactly Density(record), and (with a
shared measure, so no conversion) each element equals
weights[i] * channels[i].Density(record).  This is the surface the optimizer
uses instead of re-evaluating every channel from Python.
"""
import math
import os

import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import _util

PT = siren.dataclasses.ParticleType
M_V1 = 0.017
_NuMu = dc.Particle.ParticleType.NuMu


def _make_record(m_chi=0.008, E_V1=0.030, vertex=(0.0, 0.0, 0.0)):
    pz = math.sqrt(max(E_V1 * E_V1 - M_V1 * M_V1, 0.0))
    rec = siren.dataclasses.InteractionRecord()
    rec.signature.primary_type = PT.N4
    rec.signature.target_type = PT.Decay
    rec.signature.secondary_types = [PT.NuLight, PT.Gamma]
    rec.primary_mass = M_V1
    rec.primary_momentum = [E_V1, 0.0, 0.0, pz]
    rec.secondary_masses = [m_chi, m_chi]
    rec.secondary_momenta = [[0, 0, 0, 0], [0, 0, 0, 0]]
    rec.secondary_helicities = [0, 0]
    rec.interaction_vertex = list(vertex)
    rec.primary_initial_position = list(vertex)
    return rec


def test_density_breakdown_sums_to_density():
    target = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 2.0)),
        1.0, 1.0, 1.0)
    iso = siren.injection.Isotropic2BodyChannel(0)
    directed = siren.injection.DetectorDirected2BodyChannel(target, 0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [iso, directed]
    mc.weights = [0.3, 0.7]

    rng = siren.utilities.SIREN_random(7)
    for _ in range(2000):
        r = _make_record()
        mc.Sample(rng, None, r)
        contribs = mc.DensityBreakdown(None, r)
        g = mc.Density(None, r)
        assert len(contribs) == len(mc.channels)
        assert sum(contribs) == g  # exact: same accumulation order as Density()
        for i, ch in enumerate(mc.channels):
            assert contribs[i] == mc.weights[i] * ch.Density(None, r)


# --------------------------------------------------------------------------- #
# VertexWeightFactors.channel_densities contract, as observed through a real  #
# assembled Injector/Weighter pair (not the bare MultiChannelPhaseSpace above).#
#                                                                              #
# ComputeVertexFactors (Weighter.cxx) never writes into channel_densities on   #
# any current code path: the DummyCrossSection assembly below has no vertex   #
# backed by a MultiChannelPhaseSpace, so the field stays at its default-       #
# constructed empty map for every vertex.  This test pins that observed        #
# current contract -- a dict on every vertex, empty here -- rather than        #
# asserting non-empty population that this assembly never produces.           #
# --------------------------------------------------------------------------- #

def _skip_unless_ccm_data():
    """Skip only when the CCM detector data files are absent."""
    try:
        det_dir = _util.get_detector_model_path("CCM")
    except ValueError as e:
        pytest.skip(f"CCM detector model path unavailable: {e}")
    missing = [p for p in (os.path.join(det_dir, "materials.dat"),
                           os.path.join(det_dir, "densities.dat"))
               if not os.path.exists(p)]
    if missing:
        pytest.skip("CCM detector data missing: " + ", ".join(missing))


def _load_ccm_detector():
    dm = detector.DetectorModel()
    det_dir = _util.get_detector_model_path("CCM")
    dm.LoadMaterialModel(os.path.join(det_dir, "materials.dat"))
    dm.LoadDetectorModel(os.path.join(det_dir, "densities.dat"))
    return dm


def _build_chain(detector_model, n_inject, seed):
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])

    primary_inj = injection.PrimaryInjectionProcess()
    primary_inj.primary_type = _NuMu
    primary_inj.interactions = int_col
    primary_inj.distributions = [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(smath.Vector3D(0, 0, 0), 25.0),
    ]

    primary_phys = injection.PhysicalProcess()
    primary_phys.primary_type = _NuMu
    primary_phys.interactions = int_col
    primary_phys.distributions = [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]

    sec_xs = interactions.DummyCrossSection()
    sec_col = interactions.InteractionCollection(_NuMu, [sec_xs])

    sec_inj = injection.SecondaryInjectionProcess()
    sec_inj.primary_type = _NuMu
    sec_inj.interactions = sec_col
    sec_inj.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    sec_phys = injection.PhysicalProcess()
    sec_phys.primary_type = _NuMu
    sec_phys.interactions = sec_col
    sec_phys.distributions = [distributions.SecondaryPhysicalVertexDistribution()]

    rand = utilities.SIREN_random(seed)
    inj = injection._Injector(n_inject, detector_model, primary_inj, [sec_inj], rand)
    inj.SetStoppingCondition(lambda tree, datum, i: datum.depth(tree) >= 1)
    weighter = injection._Weighter([inj], detector_model, primary_phys, [sec_phys])

    keepalive = (xs, int_col, sec_xs, sec_col, primary_inj, primary_phys,
                 sec_inj, sec_phys, rand)
    return inj, weighter, keepalive


def test_vertex_channel_densities_is_an_empty_dict_for_this_assembly():
    """VertexWeightFactors.channel_densities is a dict on every vertex and stays
    empty for the DummyCrossSection assembly, which has no
    MultiChannelPhaseSpace-backed vertex. The finiteness check and the final
    assertion form a tripwire: if a future change starts populating the map,
    this test fails so the check can be tightened to the populated values."""
    _skip_unless_ccm_data()
    detector_model = _load_ccm_detector()
    inj, weighter, _keepalive = _build_chain(detector_model, 2000, 2468)

    events = []
    for _ in range(2000):
        if len(events) >= 20:
            break
        try:
            ev = inj.GenerateEvent()
        except RuntimeError as err:
            if "maximum number of injection attempts" not in str(err):
                raise
            break
        if len(ev.tree) == 0:
            continue
        events.append(ev)
    assert len(events) == 20

    saw_nonempty = False
    for ev in events:
        bd = weighter.EventWeightWithBreakdown(ev)
        for v in bd.vertices:
            cd = v.channel_densities
            assert isinstance(cd, dict)
            if cd:
                saw_nonempty = True
                values = list(cd.values())
                assert all(math.isfinite(x) for x in values)

    # Documents the current, observed contract for this assembly: the
    # DummyCrossSection chain has no MultiChannelPhaseSpace-backed vertex, so
    # ComputeVertexFactors never populates channel_densities.
    assert saw_nonempty is False, (
        "channel_densities was populated by this assembly; the "
        "'observed empty' contract this test documents has changed and the "
        "assertions above should be tightened to check the populated values.")
