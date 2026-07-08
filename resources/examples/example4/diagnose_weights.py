"""Demo of the sampler-vs-density closure gauge: siren.check_closure().

Loads a DarkNewsTables model and runs check_closure, which checks absolute
normalization and shape flatness of the model's sampler against its own
density. Measures without a self-contained reference density get
shape-only, normalization=nan. See siren.check_closure and
siren.ClosureReport.
"""
import os
import siren
from siren import _util

proc_dir = os.path.join(_util.resource_package_dir(), "processes", "DarkNewsTables")
vp = _util.load_module("vp_diag3", os.path.join(proc_dir, "VectorPortal.py"))

model = vp.DarkPhotonToChiDecay(0.017, 0.008, 1.0, pdgid_V1=5922, pdgid_chi=5917)

rep = siren.check_closure(model, primary_energy=0.5, samples=2000)
print(rep)
