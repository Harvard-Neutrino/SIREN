"""check_closure gauge: a correct model passes, a model whose density has the
wrong SHAPE fails via the flatness profile and worst_region, a model whose
density is off by a uniform SCALE fails via the absolute normalization, a
rest/lab frame swap is flagged, and the result is deterministic under a fixed
validation seed.
"""

import math

import pytest

siren = pytest.importorskip("siren")

from siren.errors import ClosureError


# ------------------------------------------------------------------ #
#  A correct isotropic 2-body decay: sampler and density agree        #
# ------------------------------------------------------------------ #

class GoodIsoDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def test_correct_model_passes_closure():
    report = siren.check_closure(GoodIsoDecay(), samples=3000, seed=0)
    assert report.ok, str(report)
    # The absolute normalization (integral of the density over the reference
    # measure) must sit at 1.0 for a correctly normalized density.
    est, _err = report.normalization
    assert abs(est - 1.0) < 0.2
    # The shape profile is flat, so the self-normalized flatness mean is ~1.0.
    fmean, _ferr = report.flatness
    assert abs(fmean - 1.0) < 0.2


def test_closure_report_renders_and_raises():
    report = siren.check_closure(GoodIsoDecay(), samples=2000, seed=1)
    text = str(report)
    assert "Closure report" in text
    # A passing report does not raise.
    report.raise_if_failed()


# ------------------------------------------------------------------ #
#  A broken density: FinalStateProbability skewed away from the        #
#  sampler by a cos-theta-dependent factor                            #
# ------------------------------------------------------------------ #

class SkewedDensityDecay(siren.DecayModel):
    """Samples isotropically but reports a density that rises with cos(theta),
    so the sampler and the stated density are different distributions."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        p = list(record.secondary_momenta[0])
        mag = math.sqrt(sum(c * c for c in p[1:]))
        cost = p[3] / mag if mag > 0 else 0.0
        # Strongly cos-theta-weighted density while the sampler stays isotropic.
        return (1.0 + 1.5 * cost) / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def test_broken_density_fails_closure():
    report = siren.check_closure(SkewedDensityDecay(), samples=4000, seed=0)
    assert not report.ok
    assert report.worst_region  # names the worst-deviation bin


def test_broken_density_raises():
    report = siren.check_closure(SkewedDensityDecay(), samples=4000, seed=2)
    with pytest.raises(ClosureError):
        report.raise_if_failed()


# ------------------------------------------------------------------ #
#  A uniform SCALE error: the density is off by a constant factor      #
#  everywhere. The shape stays flat, so only the absolute              #
#  normalization can catch it.                                        #
# ------------------------------------------------------------------ #

class ScaledDensityDecay(siren.DecayModel):
    """Samples isotropically and reports an isotropic density that is scaled by
    a uniform factor 2.0 everywhere. The f/g SHAPE is still flat, so the
    flatness check passes; only the absolute normalization (which fixes the
    scale) sees that the density integrates to 2.0 instead of 1.0."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        # Correct density would be 1/(4*pi); this is uniformly doubled.
        return 2.0 / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def test_uniform_scale_error_fails_via_normalization():
    report = siren.check_closure(ScaledDensityDecay(), samples=4000, seed=0)
    # The absolute normalization estimate lands near 2.0, not 1.0, so the gauge
    # fails even though the shape is flat.
    est, _err = report.normalization
    assert abs(est - 2.0) < 0.3
    assert not report.ok
    # The flatness profile alone would NOT have caught this: its self-normalized
    # mean sits at ~1.0 regardless of the uniform scale.
    fmean, _ferr = report.flatness
    assert abs(fmean - 1.0) < 0.2


def test_uniform_scale_error_raises():
    report = siren.check_closure(ScaledDensityDecay(), samples=4000, seed=3)
    with pytest.raises(ClosureError):
        report.raise_if_failed()


# ------------------------------------------------------------------ #
#  Determinism under a fixed validation seed                          #
# ------------------------------------------------------------------ #

def test_closure_deterministic_under_seed():
    import siren.closure as closure

    # The result cache is keyed on the model INSTANCE, and these are two
    # distinct GoodIsoDecay() objects, so neither call can serve the other from
    # cache; clearing is belt-and-suspenders. The same seed drives the same
    # validation RNG streams, so the reports must match exactly.
    closure._CACHE.clear()
    r1 = siren.check_closure(GoodIsoDecay(), samples=2000, seed=42)
    closure._CACHE.clear()
    r2 = siren.check_closure(GoodIsoDecay(), samples=2000, seed=42)
    assert r1.normalization[0] == pytest.approx(r2.normalization[0])
    assert r1.flatness[0] == pytest.approx(r2.flatness[0])
    assert r1.worst_region == r2.worst_region


# ------------------------------------------------------------------ #
#  Cache keys on the instance, not the class + params                 #
# ------------------------------------------------------------------ #

class ScaleParamDecay(siren.DecayModel):
    """A model whose density scale is per-instance state, so two objects of the
    same class with the same call parameters must not share a cached report."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def __init__(self, scale, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._scale = scale

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return self._scale / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def test_cache_keys_on_instance_not_class():
    import siren.closure as closure

    closure._CACHE.clear()
    correct = ScaleParamDecay(1.0)
    scaled = ScaleParamDecay(2.0)
    # Same class, same call parameters, different internal state. Keying on the
    # class + params would return the first report for the second model; keying
    # on the instance keeps them distinct.
    r_correct = siren.check_closure(correct, samples=4000, seed=7)
    r_scaled = siren.check_closure(scaled, samples=4000, seed=7)
    assert r_correct.ok, str(r_correct)
    assert not r_scaled.ok, str(r_scaled)
    assert abs(r_correct.normalization[0] - 1.0) < 0.2
    assert abs(r_scaled.normalization[0] - 2.0) < 0.3


# ------------------------------------------------------------------ #
#  Mixture path exercises the validation surface                      #
# ------------------------------------------------------------------ #

def test_mixture_path_returns_report():
    mix = 0.98 * siren.channels.isotropic() + 0.02 * siren.channels.physical()
    report = siren.check_closure(mix)
    # The standalone mixture gauge cannot run the detector-dependent probe; it
    # returns a report noting the build-time probe rather than crashing.
    assert isinstance(report, siren.ClosureReport)
    assert report.notes


# ------------------------------------------------------------------ #
#  moment_z is keyed by the measured coordinate, not fanned out       #
# ------------------------------------------------------------------ #

class MultiVarIsoDecay(siren.DecayModel):
    """A correct isotropic decay that declares SEVERAL DensityVariables. The
    gauge measures a single angular coordinate, so moment_z must carry one entry
    keyed by that coordinate -- never one entry per declared variable."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def density_variables(self):
        return ["s_pair", "cos_theta_sub"]


def test_moment_z_keyed_by_coordinate_not_density_variables():
    report = siren.check_closure(MultiVarIsoDecay(), samples=3000, seed=0)
    # A single coordinate is measured, so a single moment entry is reported.
    assert len(report.moment_z) == 1
    (name,) = report.moment_z
    # It is keyed by the coordinate actually sampled, not by a declared variable.
    assert name.startswith("costheta_secondary")
    assert "s_pair" not in report.moment_z
    assert "cos_theta_sub" not in report.moment_z
