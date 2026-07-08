"""particles.define registers BSM names with a collision policy."""
import pytest

import siren
from siren import particles
from siren.errors import ConfigurationError

# A PDG code well outside the built-in enum range, used only by these tests.
_PDG = 9000042
_PDG_OTHER = 9000043


def test_define_registers_resolvable_name():
    pt = particles.define("TestChiA", _PDG, 0.05)
    assert particles.resolve("TestChiA") == pt
    assert particles.resolve("TestChiA") == siren.dataclasses.ParticleType(_PDG)


def test_define_is_idempotent_on_identical_values():
    pt1 = particles.define("TestChiB", _PDG_OTHER, 0.10)
    pt2 = particles.define("TestChiB", _PDG_OTHER, 0.10)
    assert pt1 == pt2


def test_define_name_collision_with_different_values_raises():
    particles.define("TestChiC", 9000050, 0.20)
    with pytest.raises(ConfigurationError):
        particles.define("TestChiC", 9000050, 0.99)  # same name, different mass


def test_define_pdg_collision_under_new_name_raises():
    particles.define("TestChiD", 9000060, 0.30)
    with pytest.raises(ConfigurationError):
        particles.define("TestChiE", 9000060, 0.30)  # same pdg, different name


def test_define_builtin_name_with_own_pdg_is_noop():
    numu_pdg = int(siren.dataclasses.ParticleType.NuMu)
    result = particles.define("NuMu", numu_pdg, 0.0)
    assert result == siren.dataclasses.ParticleType.NuMu


def test_define_builtin_name_with_different_pdg_raises():
    with pytest.raises(ConfigurationError):
        particles.define("NuMu", 9000099, 0.0)


def test_define_new_name_reusing_builtin_pdg_raises():
    numu_pdg = int(siren.dataclasses.ParticleType.NuMu)
    with pytest.raises(ConfigurationError):
        particles.define("MyNuMuAlias", numu_pdg, 0.0)


def test_define_builtin_name_records_mass_and_rejects_mismatch():
    """A built-in name's first define() records its mass; a later different
    mass with the same pdg is rejected, not silently dropped."""
    gamma_pdg = int(siren.dataclasses.ParticleType.Gamma)
    result = particles.define("Gamma", gamma_pdg, 0.0)
    assert result == siren.dataclasses.ParticleType.Gamma
    # The supplied mass is now on record for the built-in name.
    assert particles._name_to_mass.get("Gamma") == 0.0
    # Re-registering the same pdg with a different mass is a collision.
    with pytest.raises(ConfigurationError):
        particles.define("Gamma", gamma_pdg, 5.0)
