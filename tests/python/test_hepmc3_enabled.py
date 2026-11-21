"""CI guard: with SIREN_TEST_REQUIRE_HEPMC3=1, fail (not skip) if HepMC3 support
is missing, so a find_package regression cannot silently green a HepMC3-enabled
CI job. Unset, this test skips and has no effect on other environments.
"""
import os

import pytest

siren = pytest.importorskip("siren")

REQUIRE_HEPMC3 = os.environ.get("SIREN_TEST_REQUIRE_HEPMC3") == "1"


@pytest.mark.skipif(
    not REQUIRE_HEPMC3,
    reason="SIREN_TEST_REQUIRE_HEPMC3 not set; HepMC3 support not required here",
)
def test_hepmc3_support_is_enabled(tmp_path):
    from siren import dataclasses as dc
    from siren import hepmc3

    record = dc.InteractionRecord()
    sig = record.signature
    sig.primary_type = dc.ParticleType.NuMu
    sig.target_type = dc.ParticleType.EPlus
    sig.secondary_types = [dc.ParticleType.NuMu]
    record.signature = sig
    record.primary_momentum = [10.0, 0.0, 0.0, 10.0]
    record.target_mass = 0.000511
    record.interaction_vertex = [1.0, 2.0, 3.0]
    record.interaction_time = 5.0
    record.secondary_ids = [dc.ParticleID(1, 1)]
    record.secondary_momenta = [[10.0, 0.0, 0.0, 10.0]]
    record.secondary_masses = [0.0]
    record.secondary_helicities = [-1.0]
    record.secondary_times = [5.0]

    tree = dc.InteractionTree()
    tree.add_entry(record, None)
    tree.header.event_number = 1
    tree.header.weights = [1.0]

    out = str(tmp_path / "enabled_check.hepmc3")
    # No try/except here: SIREN_TEST_REQUIRE_HEPMC3=1 means a RuntimeError
    # ("...without HepMC3...") must fail this test, not skip it.
    hepmc3.SaveInteractionTreesAsHepMC3([tree], out)

    loaded = hepmc3.LoadInteractionTreesFromHepMC3(out)
    assert len(loaded) == 1
