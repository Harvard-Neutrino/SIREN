"""Smoke tests for the Constants pybinding submodule."""
import pytest

siren = pytest.importorskip("siren")


@pytest.fixture(scope="module")
def Constants():
    from siren.utilities import Constants
    return Constants


class TestConstantsBinding:
    """Verify that the Constants submodule is importable and exposes all
    expected attributes with reasonable values."""

    def test_import(self):
        from siren.utilities import Constants
        assert Constants is not None

    def test_mathematical_constants(self, Constants):
        assert Constants.pi == pytest.approx(3.14159265358979, rel=1e-10)
        assert Constants.tau == pytest.approx(2 * Constants.pi, rel=1e-14)
        assert Constants.degrees == pytest.approx(Constants.pi / 180.0, rel=1e-14)

    def test_particle_masses_positive(self, Constants):
        for name in ("protonMass", "neutronMass", "electronMass",
                     "muonMass", "tauMass", "wMass", "zMass"):
            val = getattr(Constants, name)
            assert val > 0, f"{name} should be positive, got {val}"

    def test_proton_mass_value(self, Constants):
        assert Constants.protonMass == pytest.approx(0.938272, rel=1e-4)

    def test_electron_mass_value(self, Constants):
        assert Constants.electronMass == pytest.approx(0.000511, rel=1e-2)

    def test_fermi_constant_positive(self, Constants):
        assert Constants.FermiConstant > 0

    def test_speed_of_light(self, Constants):
        assert Constants.c > 0

    def test_unit_conversions(self, Constants):
        # GeV is the base energy unit
        assert Constants.GeV > 0
        assert Constants.MeV == pytest.approx(Constants.GeV * 1e-3, rel=1e-10)
        assert Constants.TeV == pytest.approx(Constants.GeV * 1e3, rel=1e-10)

    def test_gweak_exposed(self, Constants):
        assert hasattr(Constants, "gweak"), "gweak constant is missing from pybinding"
        assert Constants.gweak == pytest.approx(0.64, rel=1e-6)

    def test_ckm_elements_exposed(self, Constants):
        """All nine CKM matrix elements should be bound."""
        ckm_names = ("Vud", "Vus", "Vub",
                     "Vcd", "Vcs", "Vcb",
                     "Vtd", "Vts", "Vtb")
        for name in ckm_names:
            val = getattr(Constants, name, None)
            assert val is not None, f"CKM element {name} is missing"
            assert 0 < val <= 1, f"CKM element {name}={val} out of [0,1] range"

    def test_ckm_unitarity_first_row(self, Constants):
        """First row of CKM matrix should satisfy unitarity: sum |V|^2 ~ 1."""
        row_sum = Constants.Vud**2 + Constants.Vus**2 + Constants.Vub**2
        assert row_sum == pytest.approx(1.0, abs=0.01)

    def test_weinberg_angle(self, Constants):
        # sin^2(theta_W) ~ 0.231
        assert 0.1 < Constants.thetaWeinberg < 0.4

    def test_isoscalar_mass(self, Constants):
        """isoscalarMass is used by processes.py at runtime; verify it exists."""
        assert hasattr(Constants, "isoscalarMass")
        assert Constants.isoscalarMass > 0
