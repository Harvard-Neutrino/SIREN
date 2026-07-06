"""Smoke test run by cibuildwheel against the freshly built+repaired wheel.

Confirms that the siren.hepmc3 extension imports from the installed wheel -- i.e.
that libHepMC3 was actually bundled into the wheel by the repair step -- and that
a tree round-trips through it. This fails hard (no skip) so a missing or
unbundled HepMC3 library turns the wheel red instead of passing silently.
"""
import os
import tempfile

import siren
from siren import hepmc3
from siren import dataclasses as d

rec = d.InteractionRecord()
sig = rec.signature
sig.primary_type = d.ParticleType.NuMu
sig.secondary_types = [d.ParticleType.NuMu]
rec.signature = sig
rec.primary_momentum = [1.0, 0.0, 0.0, 1.0]
rec.secondary_ids = [d.ParticleID(1, 1)]
rec.secondary_momenta = [[1.0, 0.0, 0.0, 1.0]]
rec.secondary_masses = [0.0]
rec.secondary_helicities = [0.0]

tree = d.InteractionTree()
tree.add_entry(rec, None)

path = os.path.join(tempfile.mkdtemp(), "wheel_smoke.hepmc3")
hepmc3.SaveInteractionTreesAsHepMC3([tree], path)
loaded = hepmc3.LoadInteractionTreesFromHepMC3(path)
assert len(loaded) == 1, loaded
assert len(loaded[0].tree) == 1, loaded[0].tree

print("siren.hepmc3 wheel smoke test OK")
