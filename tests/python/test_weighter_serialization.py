"""Weighter archive header, headerless-legacy fallback, and save guard.

Mirrors test_injector_serialization.py for the Weighter. SaveWeighter writes a
magic+version header tied to the class version; LoadWeighter reads it into a
temporary (a failed parse cannot half-mutate the live weighter) and still reads
headerless version-0 archives; a missing file names the path. The Python
Weighter.save round-trips phase-space maps and refuses configurations that
would silently reload with different physics (Python-subclass models and
stopping conditions), matching Injector.save.

Data-free and fixed-seed throughout except the one CCM-guarded weight-match
round-trip, which skips when the CCM detector data files are absent.
"""
import os
import struct

import pytest

import siren
from siren import dataclasses as dc
from siren import injection
from siren import interactions
from siren import distributions
from siren import detector
from siren import math as smath
from siren import utilities
from siren import errors
from siren import _util

_NuMu = dc.Particle.ParticleType.NuMu

# Header word SaveWeighter stamps ("SWGT"), and the class version it is tied to.
_MAGIC = 0x53575447
_VERSION = 0


def _inj_distributions(max_distance):
    return [
        distributions.PrimaryMass(0),
        distributions.PowerLaw(2.0, 0.5, 5.0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
        distributions.PointSourcePositionDistribution(
            smath.Vector3D(0, 0, 0), max_distance),
    ]


def _phys_distributions():
    return [
        distributions.PrimaryMass(0),
        distributions.PrimaryNeutrinoHelicityDistribution(),
        distributions.IsotropicDirection(),
    ]


def _bare_raw_weighter(seed=3):
    """A raw _Injector + _Weighter over a bare detector (no data files).

    Enough to serialize (the empty detector and matched processes archive and
    round-trip); not enough to generate/weight events. Returns
    (detector, injector, weighter, keepalive)."""
    dm = detector.DetectorModel()
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])

    p_inj = injection.PrimaryInjectionProcess()
    p_inj.primary_type = _NuMu
    p_inj.interactions = int_col
    p_inj.distributions = _inj_distributions(0.0)

    p_phys = injection.PhysicalProcess()
    p_phys.primary_type = _NuMu
    p_phys.interactions = int_col
    p_phys.distributions = _phys_distributions()

    inj = injection._Injector(5, dm, p_inj, utilities.SIREN_random(seed))
    weighter = injection._Weighter([inj], dm, p_phys, [])
    return dm, inj, weighter, (xs, int_col, p_inj, p_phys)


# ------------------------------------------------------------------ #
#  magic + version header                                         #
# ------------------------------------------------------------------ #

def test_saved_archive_has_magic_version_header(tmp_path):
    """The first 8 bytes of a .siren_weighter are (magic, version) little-endian,
    with the version tied to the registered class version."""
    _dm, _inj, weighter, _keep = _bare_raw_weighter()
    base = str(tmp_path / "w")
    weighter.SaveWeighter(base)
    path = base + ".siren_weighter"
    assert os.path.exists(path)
    with open(path, "rb") as f:
        head = f.read(8)
    assert len(head) == 8
    magic, version = struct.unpack("<II", head)
    assert (magic, version) == (_MAGIC, _VERSION)


def test_headered_archive_roundtrips(tmp_path):
    """A headered archive reloads through the (injectors, filename) ctor; passing
    an empty injector list keeps the archived injectors, proving the payload
    parsed."""
    _dm, _inj, weighter, _keep = _bare_raw_weighter()
    base = str(tmp_path / "w")
    weighter.SaveWeighter(base)
    reloaded = injection._Weighter([], base)
    assert len(reloaded.GetInjectors()) == 1
    assert reloaded.GetDetectorModel() is not None


# ------------------------------------------------------------------ #
#  headerless legacy fallback                                     #
# ------------------------------------------------------------------ #

def test_headerless_archive_still_loads(tmp_path):
    """Stripping the 8-byte header leaves a legacy version-0 archive that
    LoadWeighter still parses (the headerless fallback)."""
    _dm, _inj, weighter, _keep = _bare_raw_weighter()
    base = str(tmp_path / "w")
    weighter.SaveWeighter(base)
    with open(base + ".siren_weighter", "rb") as f:
        full = f.read()

    legacy_base = str(tmp_path / "legacy")
    with open(legacy_base + ".siren_weighter", "wb") as f:
        f.write(full[8:])  # drop magic+version -> raw version-0 payload

    reloaded = injection._Weighter([], legacy_base)
    assert len(reloaded.GetInjectors()) == 1


# ------------------------------------------------------------------ #
#  missing file names the path                                    #
# ------------------------------------------------------------------ #

def test_load_missing_file_names_path(tmp_path):
    """Weighter.load on a nonexistent archive raises naming the file path,
    instead of a cryptic cereal stream error."""
    missing = str(tmp_path / "does_not_exist")
    w = injection.Weighter()  # no injectors configured -> loads with []
    with pytest.raises(RuntimeError) as exc:
        w.load(missing)
    assert (missing + ".siren_weighter") in str(exc.value)


# ------------------------------------------------------------------ #
#  phase-space round-trip and remaining save guards               #
# ------------------------------------------------------------------ #

def test_save_roundtrips_phase_space_map_injector(tmp_path):
    """Weighter.save preserves a raw injector's phase-space map."""
    dm = detector.DetectorModel()
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])
    p_inj = injection.PrimaryInjectionProcess()
    p_inj.primary_type = _NuMu
    p_inj.interactions = int_col
    p_inj.distributions = _inj_distributions(0.0)

    # Attach a phase-space map to the injection process' signature.
    signature = xs.GetPossibleSignatures()[0]
    mc = injection.MultiChannelPhaseSpace()
    mc.channels = [injection.Isotropic2BodyChannel(0)]
    mc.weights = [1.0]
    p_inj.SetPhaseSpace(signature, mc)
    assert p_inj.HasAnyPhaseSpace()

    inj = injection._Injector(5, dm, p_inj, utilities.SIREN_random(3))
    weighter = injection.Weighter(
        injectors=[inj], detector_model=dm, primary_type=_NuMu,
        primary_interactions=[xs],
        primary_physical_distributions=[distributions.PrimaryMass(0),
                                        distributions.IsotropicDirection()])

    base = str(tmp_path / "w")
    weighter.save(base)
    assert os.path.exists(base + ".siren_weighter")

    reloaded = injection.Weighter()
    reloaded.load(base)
    loaded_injectors = reloaded.engine.GetInjectors()
    assert len(loaded_injectors) == 1
    loaded_process = loaded_injectors[0].GetPrimaryProcess()
    assert loaded_process.HasAnyPhaseSpace()
    loaded_map = loaded_process.GetPhaseSpaceMap()
    assert list(loaded_map) == [signature]
    loaded_mixture = loaded_map[signature]
    assert len(loaded_mixture.channels) == len(mc.channels)
    assert list(loaded_mixture.weights) == list(mc.weights)


def test_save_delegates_to_injector_guard(tmp_path):
    """The weighter's save defers to each injector's own guard: a python Injector
    with a stopping condition is refused, its offenders prefixed 'injector: '."""
    dm = detector.DetectorModel()
    pyinj = injection.Injector(
        detector=dm, number_of_events=2, seed=1, primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_inj_distributions(0.0),
        stopping_condition=lambda datum, i: True)
    pyinj._build()

    weighter = injection.Weighter(
        pyinj, primary_physical=[distributions.PrimaryMass(0),
                                 distributions.IsotropicDirection()])
    base = str(tmp_path / "w")
    with pytest.raises(errors.NotSerializableError) as exc:
        weighter.save(base)
    assert any(o.startswith("injector: ") for o in exc.value.offenders)
    assert not os.path.exists(base + ".siren_weighter")


def test_save_refuses_python_subclass_distribution(tmp_path):
    """A physical distribution written as a python subclass of a pybind base
    does not survive the cereal archive; save refuses it and names the class."""
    dm = detector.DetectorModel()

    class UnitWeight(distributions.WeightableDistribution):
        def __init__(self):
            super().__init__()

        def GenerationProbability(self, detector_model, ic, record):
            return 1.0

    pyinj = injection.Injector(
        detector=dm, number_of_events=2, seed=1, primary_type=_NuMu,
        primary_interactions=[interactions.DummyCrossSection()],
        primary_injection_distributions=_inj_distributions(0.0))
    pyinj._build()

    weighter = injection.Weighter(
        pyinj, primary_physical=[distributions.PrimaryMass(0),
                                 distributions.IsotropicDirection(),
                                 UnitWeight()])
    base = str(tmp_path / "w")
    with pytest.raises(errors.NotSerializableError) as exc:
        weighter.save(base)
    assert any("is a Python subclass" in o and "UnitWeight" in o
               for o in exc.value.offenders)
    assert not os.path.exists(base + ".siren_weighter")


# ------------------------------------------------------------------ #
#  reloaded weights match the original (CCM-guarded)              #
# ------------------------------------------------------------------ #

def _load_ccm_detector():
    """Load the CCM detector model or skip if its data files are missing."""
    try:
        det_dir = _util.get_detector_model_path("CCM")
        materials = os.path.join(det_dir, "materials.dat")
        densities = os.path.join(det_dir, "densities.dat")
        if not (os.path.exists(materials) and os.path.exists(densities)):
            pytest.skip("CCM detector data files not present")
        dm = detector.DetectorModel()
        dm.LoadMaterialModel(materials)
        dm.LoadDetectorModel(densities)
        return dm
    except (ValueError, RuntimeError) as exc:
        pytest.skip("CCM detector model unavailable: {}".format(exc))


def test_reloaded_weights_match_original(tmp_path):
    """A guard-clean weighter saved and reloaded through the header path weights
    a real event identically to the original."""
    dm = _load_ccm_detector()
    xs = interactions.DummyCrossSection()
    int_col = interactions.InteractionCollection(_NuMu, [xs])

    p_inj = injection.PrimaryInjectionProcess()
    p_inj.primary_type = _NuMu
    p_inj.interactions = int_col
    p_inj.distributions = _inj_distributions(25.0)

    p_phys = injection.PhysicalProcess()
    p_phys.primary_type = _NuMu
    p_phys.interactions = int_col
    p_phys.distributions = _phys_distributions()

    inj = injection._Injector(50, dm, p_inj, utilities.SIREN_random(1234))
    weighter = injection._Weighter([inj], dm, p_phys, [])

    event = None
    for _ in range(200):
        ev = inj.GenerateEvent()
        if len(ev.tree) > 0:
            event = ev
            break
    assert event is not None, "no event generated"
    reference = weighter.EventWeight(event)
    assert reference > 0.0

    base = str(tmp_path / "w")
    weighter.SaveWeighter(base)
    reloaded = injection._Weighter([inj], base)  # reuse the live injector
    assert reloaded.EventWeight(event) == pytest.approx(reference, rel=1e-12)
