"""Demo of the weight-decomposition gauge: weighter.explain().

Builds the same small pion-decay chain, generates a few events, weights
them, and prints the per-vertex physical/generation factor breakdown for one
event, plus .culprit() when the total weight is unusable. See
siren.Weighter.explain, siren.Results.explain, and siren.report.WeightBreakdown.
"""
import math
import os
import siren
from siren import _util, distributions, injection
from siren.Injector import Injector
from siren.Weighter import Weighter

proc_dir = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
meson = _util.load_module("meson_diag2", os.path.join(proc_dir, "MesonProduction.py"))
pion_decay = meson.MesonThreeBodySIRENDecay(
    m_mediator=0.017, g_mu=1e-3, pdgid_meson=211, pdgid_lepton=-13,
    pdgid_neutrino=14, pdgid_mediator=5922)
sig = pion_decay.GetPossibleSignatures()[0]
fiducial = siren.geometry.Box(
    siren.geometry.Placement(siren.math.Vector3D(0, 0.3, 0)), 0.6, 0.6, 0.6)

mc = injection.MultiChannelPhaseSpace()
mc.channels = [injection.PhysicalDecayChannel(pion_decay, sig),
               injection.DetectorDirected3BodyChannel(
                   factorization=injection.ThreeBodyMode.Direct,
                   target=fiducial, directed_index=2,
                   topology=injection.PhaseSpaceTopology.Decay3Body)]
mc.weights = [0.5, 0.5]

direction = distributions.FixedDirection(siren.math.Vector3D(0, 0, 1))
energy = distributions.Monoenergetic(0.3)
injector = Injector(
    events=200, seed=1234, primary_type=sig.primary_type,
    detector_model=siren.utilities.load_detector("SBN", detector="SBND"),
    primary_interactions=[pion_decay],
    primary_injection_distributions=[
        distributions.PrimaryMass(0.13957039), energy, direction,
        distributions.PointSourcePositionDistribution([0, 0, 0], 300.0)],
    primary_phase_spaces={sig: mc})
weighter = Weighter(
    injectors=[injector], detector_model=injector.detector_model,
    primary_type=sig.primary_type, primary_interactions=[pion_decay],
    primary_physical_distributions=[energy, direction])

trees = injector.generate(5, on_shortfall="warn")
bd = weighter.explain(trees[0])
print(bd)
if not (math.isfinite(bd.total) and bd.total > 0):
    print("culprit:", bd.culprit())
