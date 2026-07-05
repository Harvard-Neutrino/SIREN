"""Flagship BSM closure: a DarkNews dipole-portal sample written to HepMC3, read
back, and reweighted -- the event weights must survive the HepMC3 round trip.

This exercises the reweighting-critical surfaces (interaction_parameters map,
helicities, momenta, geometry) for a real BSM model, plus the G.R.11 BSM
particle-number declaration (N4 == 5914).

Requires the DarkNews package and its cross-section tables, so it SKIPS cleanly
when DarkNews is unavailable or the tables cannot be built in this environment.
It is modeled on resources/examples/example2/DipolePortal_CCM.py.
"""
import os
import tempfile

import pytest

siren = pytest.importorskip("siren")
pytest.importorskip("DarkNews")


def _build_darknews_ccm(events_to_inject=2):
    """Return (events, gen_times, weighter) for a tiny CCM dipole run.

    Any failure to construct the DarkNews processes/tables raises, and the
    caller turns that into a skip -- the point of the test is the HepMC3 round
    trip, not DarkNews table generation.
    """
    import numpy as np
    from siren import utilities
    from siren._util import GenerateEvents, get_processes_model_path

    model_kwargs = {
        "m4": 0.0235, "mu_tr_mu4": 6e-7, "UD4": 0, "Umu4": 0, "epsilon": 0.0,
        "gD": 0.0, "decay_product": "photon", "noHC": True, "HNLtype": "dirac",
    }
    experiment = "CCM"
    detector_model = utilities.load_detector(experiment)
    primary_type = siren.dataclasses.Particle.ParticleType.NuMu

    table_name = f"DarkNewsTables-v{siren.utilities.darknews_version()}/"
    table_name += "Dipole_M%2.2e_mu%2.2e" % (model_kwargs["m4"], model_kwargs["mu_tr_mu4"])
    table_dir = os.path.join(get_processes_model_path("DarkNewsTables"), table_name)
    os.makedirs(table_dir, exist_ok=True)

    primary_processes, secondary_processes, _, _ = utilities.load_processes(
        "DarkNewsTables", primary_type=primary_type,
        detector_model=detector_model, table_name=table_name, **model_kwargs)

    mass_ddist = siren.distributions.PrimaryMass(0)
    edist = siren.distributions.Monoenergetic(0.02965)
    opening_angle = np.arctan(5 / 23.0)
    lower_target_origin = siren.math.Vector3D(0, 0, -0.241)
    detector_origin = siren.math.Vector3D(23, 0, -0.65)
    lower_dir = detector_origin - lower_target_origin
    lower_dir.normalize()
    lower_inj_ddist = siren.distributions.Cone(lower_dir, opening_angle)
    phys_ddist = siren.distributions.IsotropicDirection()
    max_dist = 25
    lower_pos_dist = siren.distributions.PointSourcePositionDistribution(
        lower_target_origin - detector_origin, max_dist)

    primary_injection_distributions = [mass_ddist, edist, lower_inj_ddist, lower_pos_dist]
    primary_physical_distributions = [edist, phys_ddist]

    fiducial_volume = siren.utilities.get_fiducial_volume(experiment)
    secondary_injection_distributions = {
        st: [siren.distributions.SecondaryBoundedVertexDistribution(fiducial_volume, max_dist)]
        for st in secondary_processes.keys()
    }

    injector = siren.injection.Injector()
    injector.number_of_events = events_to_inject
    injector.detector_model = detector_model
    injector.primary_type = primary_type
    injector.primary_interactions = primary_processes[primary_type]
    injector.primary_injection_distributions = primary_injection_distributions
    injector.secondary_interactions = secondary_processes
    injector.secondary_injection_distributions = secondary_injection_distributions

    events, gen_times = GenerateEvents(injector)

    weighter = siren.injection.Weighter()
    weighter.injectors = [injector]
    weighter.detector_model = detector_model
    weighter.primary_type = primary_type
    weighter.primary_interactions = primary_processes[primary_type]
    weighter.secondary_interactions = secondary_processes
    weighter.primary_physical_distributions = primary_physical_distributions
    weighter.secondary_physical_distributions = {}
    return events, gen_times, weighter


def test_darknews_hepmc3_roundtrip_reweight():
    from siren import hepmc3, _util
    try:
        events, gen_times, weighter = _build_darknews_ccm(events_to_inject=2)
    except Exception as exc:  # DarkNews tables/model unavailable in this environment
        pytest.skip(f"DarkNews setup unavailable: {exc}")

    if len(events) == 0:
        pytest.skip("no DarkNews events generated")

    out = os.path.join(tempfile.mkdtemp(), "dn.hepmc3")
    try:
        hepmc3.SaveInteractionTreesAsHepMC3(events, out)
    except RuntimeError as exc:
        if "without HepMC3" in str(exc):
            pytest.skip("SIREN built without HepMC3 support")
        raise

    # BSM particle-number declaration (G.R.11): N4 == 5914 must be present.
    with open(out) as f:
        text = f.read()
    assert "5914" in text

    loaded = _util.LoadEventsFromHepMC3(out)
    assert len(loaded) == len(events)

    # Reweighting closure: the weight recomputed from the round-tripped tree must
    # match the weight of the original in-memory event.
    for original, roundtripped in zip(events, loaded):
        w0 = weighter(original)
        w1 = weighter(roundtripped)
        assert w0 > 0
        assert abs(w1 - w0) / w0 < 1e-6
