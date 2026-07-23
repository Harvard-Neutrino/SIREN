"""MultiChannelPhaseSpace.DensityBreakdown invariant.

DensityBreakdown(record)[i] is channel i's alpha-weighted, common-measure
contribution; the elements must sum to exactly Density(record), and (with a
shared measure, so no conversion) each element equals
weights[i] * channels[i].Density(record).  This is the surface the optimizer
uses instead of re-evaluating every channel from Python.
"""
import math

import siren

PT = siren.dataclasses.ParticleType
M_V1 = 0.017


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
