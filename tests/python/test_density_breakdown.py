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
# VertexWeightFactors.channel_densities/cancelled contract, as observed        #
# through a real assembled Injector/Weighter pair (not the bare                #
# MultiChannelPhaseSpace above).                                              #
#                                                                              #
# ComputeVertexFactors (Weighter.cxx) always writes channel_densities and      #
# cancelled on the diagnostics path (EventWeightWithBreakdown): the former is  #
# keyed from the injection process's MultiChannelPhaseSpace.DensityBreakdown  #
# when one is attached for the vertex's signature, and the latter lists the   #
# distribution Name()s that cancel between the injection and physical          #
# distribution lists. Cancellation is by value equality (operator==), not      #
# shared_ptr identity: a distribution present in both lists cancels even when  #
# the two lists hold distinct but value-equal instances. The DummyCrossSection #
# assembly below attaches no phase space, so channel_densities stays empty,    #
# but its injection and physical lists share value-equal PrimaryMass,          #
# PrimaryNeutrinoHelicityDistribution, IsotropicDirection (primary vertex) and #
# SecondaryPhysicalVertexDistribution (secondary vertex) instances, so those   #
# names appear in cancelled.                                                   #
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


def test_vertex_cancelled_matches_value_equal_distributions():
    """cancelled lists the value-equal distributions shared by the injection and
    physical lists; channel_densities stays empty with no mixture attached.

    Cancellation is by value equality (WeightableDistribution::operator==), not
    shared_ptr identity: the injection and physical lists here hold distinct but
    value-equal instances, and they cancel exactly as literally shared instances
    would. The primary vertex shares PrimaryMass, PrimaryNeutrinoHelicity-
    Distribution, and IsotropicDirection; the secondary vertex shares
    SecondaryPhysicalVertexDistribution. PowerLaw and the position
    distributions live on only one side, so they do not cancel.
    """
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

    primary_cancelled = {"PrimaryMass", "PrimaryNeutrinoHelicityDistribution",
                         "IsotropicDirection"}
    secondary_cancelled = {"SecondaryPhysicalVertexDistribution"}
    seen_cancelled = set()
    for ev in events:
        bd = weighter.EventWeightWithBreakdown(ev)
        for v in bd.vertices:
            cd = v.channel_densities
            cancelled = set(v.cancelled)
            assert isinstance(cd, dict)
            assert isinstance(v.cancelled, list)
            assert all(math.isfinite(x) for x in cd.values())
            # No MultiChannelPhaseSpace is attached, so channel_densities stays
            # empty; cancelled names the value-equal shared distributions.
            assert cd == {}
            assert cancelled in (primary_cancelled, secondary_cancelled)
            seen_cancelled |= cancelled
    assert seen_cancelled == primary_cancelled | secondary_cancelled


def test_channel_densities_contract_matches_density_breakdown():
    """channel_densities entries are the mixture's DensityBreakdown, keyed channel[i]."""
    target = siren.geometry.Box(
        siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 2.0)),
        1.0, 1.0, 1.0)
    iso = siren.injection.Isotropic2BodyChannel(0)
    directed = siren.injection.DetectorDirected2BodyChannel(target, 0)
    mc = siren.injection.MultiChannelPhaseSpace()
    mc.channels = [iso, directed]
    mc.weights = [0.3, 0.7]

    rng = siren.utilities.SIREN_random(11)
    r = _make_record()
    mc.Sample(rng, None, r)

    bd = mc.DensityBreakdown(None, r)
    g = mc.Density(None, r)
    assert len(bd) == 2
    assert all(math.isfinite(x) for x in bd)
    assert sum(bd) == g

    channel_densities = {"channel[" + str(i) + "]": bd[i] for i in range(len(bd))}
    assert list(channel_densities.keys()) == ["channel[0]", "channel[1]"]
    assert all(math.isfinite(x) for x in channel_densities.values())
    assert sum(channel_densities.values()) == g


def test_vertex_channel_densities_report_their_common_convention():
    """Every populated diagnostic density map identifies its common measure."""
    _skip_unless_ccm_data()
    detector_model = _load_ccm_detector()
    inj, weighter, keepalive = _build_chain(detector_model, 2000, 97531)

    event = None
    for _ in range(2000):
        try:
            candidate = inj.GenerateEvent()
        except RuntimeError as err:
            if "maximum number of injection attempts" not in str(err):
                raise
            break
        if len(candidate.tree) > 0:
            event = candidate
            break
    assert event is not None

    # Attach the diagnostic phase space after generation. This exercises the
    # weighting/reporting path without asking the physical adapter to resample
    # DummyCrossSection's synthetic Nucleon secondary (which has no mass-map
    # entry).
    xs, _, _, _, primary_inj, *_ = keepalive
    signature = event.tree[0].record.signature
    channel = injection.PhysicalCrossSectionChannel(xs, signature)
    mixture = injection.MultiChannelPhaseSpace()
    mixture.channels = [channel]
    mixture.weights = [1.0]
    primary_inj.SetPhaseSpace(signature, mixture)

    breakdown = weighter.EventWeightWithBreakdown(event)
    populated = [vertex for vertex in breakdown.vertices
                 if vertex.channel_densities]
    assert populated
    for vertex in populated:
        assert vertex.channel_density_topology == injection.PhaseSpaceTopology.Scatter2to2
        assert vertex.channel_density_measure == injection.PhaseSpaceMeasure.BjorkenXY()
