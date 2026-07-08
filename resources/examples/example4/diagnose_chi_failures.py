"""Demo of the injection-failure gauge: injector.report().

Builds a small pion-decay injector (directed 3-body channel toward a narrow
off-axis fiducial box), generates events, and prints the per-vertex attrition
table injector.report() renders from the engine's failure ledger. See
siren.Injector.report and siren.report.InjectionReport.
"""
import os
import siren
from siren import _util, distributions, injection
from siren.Injector import Injector

proc_dir = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
meson = _util.load_module("meson_diag1", os.path.join(proc_dir, "MesonProduction.py"))
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

injector = Injector(
    events=200, seed=1234, primary_type=sig.primary_type,
    detector_model=siren.utilities.load_detector("SBN", detector="SBND"),
    primary_interactions=[pion_decay],
    primary_injection_distributions=[
        distributions.PrimaryMass(0.13957039),
        distributions.Monoenergetic(0.3),
        distributions.FixedDirection(siren.math.Vector3D(0, 0, 1)),
        distributions.PointSourcePositionDistribution([0, 0, 0], 300.0)],
    primary_phase_spaces={sig: mc})

injector.generate(20, on_shortfall="warn")
print(injector.report())
