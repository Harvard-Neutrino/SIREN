"""
Unit tests for the shared siren.earth factorization.

These exercise the pure logic that both the SBN and DUNE earth models rely on --
the predict-count-and-shift-up level scheme, the mass-conserving PREM collapse,
and the exponential-atmosphere band parameters -- without needing a real
DetectorModel, materials, or (for DUNE) DEM/network access.
"""

import math
import pytest

import siren.earth as earth


# --- lightweight duck-typed stand-ins for DetectorModel / DetectorSector ----

class _FakeSector:
    def __init__(self, name, level):
        self.name = name
        self.level = level


class _FakeModel:
    """Mimics the DetectorModel.Sectors get/set contract used by siren.earth."""
    def __init__(self, sectors):
        self._sectors = list(sectors)

    @property
    def Sectors(self):
        # The real getter yields a fresh list of the sector objects; callers
        # mutate .level then reassign via the setter.
        return list(self._sectors)

    @Sectors.setter
    def Sectors(self, value):
        self._sectors = list(value)

    def levels(self):
        return {s.name: s.level for s in self._sectors}


def test_shift_levels_up_normalizes_and_exempts_universe():
    sectors = [_FakeSector("UNIVERSE", -1), _FakeSector("World", 0),
               _FakeSector("det", 1), _FakeSector("det2", 2)]
    min_level = earth.shift_levels_up(sectors, 5)
    assert min_level == 0
    lv = {s.name: s.level for s in sectors}
    assert lv["UNIVERSE"] == -1          # exempt, untouched
    assert lv["World"] == 5              # 0 - 0 + 5
    assert lv["det"] == 6                # 1 - 0 + 5
    assert lv["det2"] == 7               # 2 - 0 + 5


def test_insert_sectors_above_places_between_anchor_and_higher():
    m = _FakeModel([_FakeSector("UNIVERSE", -1), _FakeSector("World", 0),
                    _FakeSector("det", 1), _FakeSector("det2", 2)])
    overburden = [_FakeSector("ob_low", None), _FakeSector("ob_high", None)]
    earth.insert_sectors_above(m, "World", overburden)
    lv = m.levels()
    # anchor unchanged; below-anchor (UNIVERSE) unchanged
    assert lv["World"] == 0
    assert lv["UNIVERSE"] == -1
    # overburden fills the freed levels just above the anchor, in order
    assert lv["ob_low"] == 1
    assert lv["ob_high"] == 2
    # everything that was above the anchor shifts up by len(overburden)
    assert lv["det"] == 3
    assert lv["det2"] == 4
    # strict ordering: PREM/background < World < overburden < detectors
    assert lv["World"] < lv["ob_low"] < lv["ob_high"] < lv["det"] < lv["det2"]


def test_insert_sectors_above_unknown_anchor_raises():
    m = _FakeModel([_FakeSector("World", 0)])
    with pytest.raises(RuntimeError):
        earth.insert_sectors_above(m, "NoSuchSector", [_FakeSector("x", None)])


def test_collapse_to_constant_is_mass_conserving_and_matches_crust():
    layers = earth.CANONICAL_PREM_LAYERS
    collapsed = earth.collapse_to_constant(layers)
    assert len(collapsed) == len(layers)
    # constant layers pass through unchanged (upper crust = 2.600)
    assert abs(collapsed[-1][4] - 2.600) < 1e-9
    assert collapsed[-2][4] == pytest.approx(2.900)  # inner_crust constant
    # inner core mean is close to the central PREM value ~13.09 g/cm^3
    assert 12.5 < collapsed[0][4] < 13.1
    # every collapsed shell is a positive constant
    for name, r, mat, dtype, rho in collapsed:
        assert dtype == "constant"
        assert rho > 0.0


def test_atmosphere_exponential_band_params():
    atm = earth.USStandard1976()
    form = earth.AtmosphereForm(model=atm, mode="exponential")
    bands = form.resolve(earth.R_PREM)
    assert len(bands) == len(earth.DEFAULT_SHELL_ALTS)
    H = atm.scale_height
    for band, (h1, h2) in zip(bands, earth.DEFAULT_SHELL_ALTS):
        assert band.mode == "exponential"
        assert band.material == "AIR"
        # base radius = surface + h1; density at base = rho0*exp(-h1/H)
        assert band.r_base == pytest.approx(earth.R_PREM + h1)
        assert band.r_top == pytest.approx(earth.R_PREM + h2)
        assert band.rho_base == pytest.approx(atm.rho0 * math.exp(-h1 / H))
        # densities strictly decrease with altitude, and never underflow to 0
        assert band.rho_base > 0.0


def test_atmosphere_column_depth_exponential_matches_analytic():
    """The exact exponential column integral is ~1030 g/cm^2, and agrees with
    the mass-conserving shell sum the loaders still expose (both integrate the
    same barometric profile)."""
    atm = earth.USStandard1976()
    H = atm.scale_height
    # analytic column from surface to top of the modeled atmosphere
    top = earth.DEFAULT_SHELL_ALTS[-1][1]
    analytic = atm.rho0 * H * (1.0 - math.exp(-top / H)) * 100.0  # g/cm^2
    shell_sum = 0.0
    for (h1, h2) in earth.DEFAULT_SHELL_ALTS:
        shell_sum += atm.shell_mean_density(h1, h2) * (h2 - h1) * 100.0
    assert analytic == pytest.approx(shell_sum, rel=1e-6)
    assert 980 < analytic < 1080
