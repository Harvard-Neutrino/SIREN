"""Authoring-base contracts: derived signatures, closure-by-construction
default sampler, the recursion guard, near-miss override rejection, and the
override audit at Injector build.
"""

import math

import pytest

siren = pytest.importorskip("siren")

import siren.models as models
from siren.errors import ConfigurationError


P = siren.dataclasses.ParticleType


def _template_record(signature, energy=0.05, primary_mass=0.02):
    from siren import _validation
    rec = _validation.build_template_record(
        signature, primary_mass=primary_mass, energy=energy)
    rec.secondary_masses = [0.0 for _ in signature.secondary_types]
    return rec


def _template_csdr(signature, energy=0.05, primary_mass=0.02):
    # Keep the source record alive alongside the CSDR: the CSDR holds a
    # reference to it, so the record must outlive the CSDR.
    rec = _template_record(signature, energy=energy, primary_mass=primary_mass)
    csdr = siren.dataclasses.CrossSectionDistributionRecord(rec)
    return csdr, rec


# ------------------------------------------------------------------ #
#  A minimal 2-body decay on the default (SolidAngleRest) sampler      #
# ------------------------------------------------------------------ #

class IsoDecay(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleRest()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def density_variables(self):
        return "cost"


def test_derived_signature_and_topology():
    m = IsoDecay()
    sigs = m.GetPossibleSignatures()
    assert len(sigs) == 1
    sig = sigs[0]
    assert sig.primary_type == P.N4
    assert sig.target_type == P.Decay
    assert list(sig.secondary_types) == [P.NuLight, P.Gamma]
    assert m.Topology() == siren.Topology.Decay2Body
    assert m.GetPossibleSignaturesFromParent(P.N4)
    assert m.GetPossibleSignaturesFromParent(P.EMinus) == []


def test_density_variables_normalized_to_list():
    assert IsoDecay().DensityVariables() == ["cost"]


def test_width_overload_pair_and_fsp():
    m = IsoDecay()
    assert m.TotalDecayWidth(P.N4) == 1.0
    assert m.TotalDecayWidth(P.EMinus) == 0.0
    assert m.TotalDecayWidthAllFinalStates(P.N4) == 1.0
    rec = _template_record(m.GetPossibleSignatures()[0])
    # FinalStateProbability = differential / total.
    assert m.FinalStateProbability(rec) == pytest.approx(1.0 / (4.0 * math.pi))


def test_default_sampler_writes_secondaries():
    m = IsoDecay()
    csdr, _rec = _template_csdr(m.GetPossibleSignatures()[0])
    rng = siren.utilities.SIREN_random(1234)
    m.SampleFinalState(csdr, rng)
    secs = csdr.get_secondary_particle_records()
    assert len(secs) == 2
    # Both daughters got a nonzero four-momentum from the engine channel.
    for spr in secs:
        p = list(spr.four_momentum)
        assert any(abs(c) > 0 for c in p)


def test_recursion_canary_completes():
    """The default sampler samples the engine channel, never recursing back
    into SampleFinalState through PhysicalDecayChannel."""
    m = IsoDecay()
    sig = m.GetPossibleSignatures()[0]
    rng = siren.utilities.SIREN_random(7)
    for _ in range(500):
        csdr, _rec = _template_csdr(sig)
        m.SampleFinalState(csdr, rng)


# ------------------------------------------------------------------ #
#  A declared measure with no self-contained channel                  #
# ------------------------------------------------------------------ #

class LabDecayNoSampler(siren.DecayModel):
    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleLab()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)


def test_default_sampler_rejects_unsupported_measure():
    m = LabDecayNoSampler()
    csdr, _rec = _template_csdr(m.GetPossibleSignatures()[0])
    rng = siren.utilities.SIREN_random(3)
    with pytest.raises(ConfigurationError):
        m.SampleFinalState(csdr, rng)


class LabDecayWithSampler(siren.DecayModel):
    """A non-default measure with an explicit sample() override; the engine
    entry point SampleFinalState must reach the override."""

    parent = "N4"
    daughters = ("NuLight", "Gamma")
    measure = siren.Measure.SolidAngleLab()

    def total_width(self):
        return 1.0

    def differential_width(self, record):
        return 1.0 / (4.0 * math.pi)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sampled = False

    def sample(self, record, random):
        self.sampled = True
        secs = record.get_secondary_particle_records()
        for spr in secs:
            spr.four_momentum = [0.01, 0.0, 0.0, 0.01]
            spr.mass = 0.0


def test_sample_override_reached_through_engine_entry_point():
    m = LabDecayWithSampler()
    csdr, _rec = _template_csdr(m.GetPossibleSignatures()[0])
    rng = siren.utilities.SIREN_random(5)
    # SampleFinalState (the C++ virtual entry point) dispatches to sample().
    m.SampleFinalState(csdr, rng)
    assert m.sampled


def test_audit_accepts_overridden_sampler_for_unsupported_measure():
    from siren import _validation
    _validation.audit_overrides([LabDecayWithSampler()])


# ------------------------------------------------------------------ #
#  __init_subclass__ near-miss override rejection                     #
# ------------------------------------------------------------------ #

def test_near_miss_override_rejected_at_definition():
    with pytest.raises(ConfigurationError):
        class Typo(siren.DecayModel):
            parent = "N4"
            daughters = ("NuLight", "Gamma")
            measure = siren.Measure.SolidAngleRest()

            def total_width(self):
                return 1.0

            def diferential_width(self, record):   # noqa: intentional typo
                return 1.0


# ------------------------------------------------------------------ #
#  Cross-section base                                                  #
# ------------------------------------------------------------------ #

class ElasticXS(siren.CrossSectionModel):
    primary = "NuMu"
    target = "PPlus"
    finals = ("NuMu", "PPlus")
    measure = siren.Measure.SolidAngleRest()

    def total_xs(self, record):
        return 1e-38

    def differential_xs(self, record):
        return 1e-38 / (4.0 * math.pi)

    def density_variables(self):
        return ["cost"]


def test_cross_section_signature_and_targets():
    m = ElasticXS()
    assert m.GetPossiblePrimaries() == [P.NuMu]
    assert m.GetPossibleTargets() == [P.PPlus]
    assert m.GetPossibleTargetsFromPrimary(P.NuMu) == [P.PPlus]
    assert m.GetPossibleTargetsFromPrimary(P.NuE) == []
    assert m.Topology() == siren.Topology.Scatter2to2
    sig = m.GetPossibleSignatures()[0]
    assert list(sig.secondary_types) == [P.NuMu, P.PPlus]


# ------------------------------------------------------------------ #
#  audit_overrides on a legacy-style typo'd direct subclass           #
# ------------------------------------------------------------------ #

class _LegacyBrokenDecay(siren.interactions.Decay):
    """A direct pyDecay subclass that misspells DifferentialDecayWidth."""

    def __init__(self):
        siren.interactions.Decay.__init__(self)
        self._sig = siren.dataclasses.InteractionSignature()
        self._sig.primary_type = P.N4
        self._sig.target_type = P.Decay
        self._sig.secondary_types = [P.NuLight, P.Gamma]

    def equal(self, other):
        return self is other

    def GetPossibleSignatures(self):
        return [self._sig]

    def GetPossibleSignaturesFromParent(self, primary_type):
        return [self._sig]

    def TotalDecayWidth(self, arg):
        return 1.0

    def TotalDecayWidthAllFinalStates(self, arg):
        return 1.0

    # DifferentialDecayWidth intentionally MISSING (the silent pure-virtual case)

    def FinalStateProbability(self, record):
        return 1.0

    def DensityVariables(self):
        return ["cost"]

    def SampleFinalState(self, csdr, random):
        pass


def test_audit_overrides_catches_legacy_typo():
    from siren import _validation
    with pytest.raises(ConfigurationError):
        _validation.audit_overrides([_LegacyBrokenDecay()])


def test_audit_overrides_passes_complete_model():
    from siren import _validation
    _validation.audit_overrides([IsoDecay()])
    _validation.audit_overrides([ElasticXS()])


def test_audit_flags_default_sampler_without_channel():
    from siren import _validation
    with pytest.raises(ConfigurationError):
        _validation.audit_overrides([LabDecayNoSampler()])


# ------------------------------------------------------------------ #
#  Base-parameterizable factory over the DarkNews C++ base            #
# ------------------------------------------------------------------ #

def test_decay_model_base_over_darknews_base():
    DarkNewsDecay = pytest.importorskip("siren.interactions").DarkNewsDecay
    Base = models.decay_model_base(base=DarkNewsDecay)

    class DNStyle(Base):
        parent = "N4"
        daughters = ("NuLight", "Gamma")
        measure = siren.Measure.SolidAngleRest()

        def total_width(self):
            return 1.0

        def differential_width(self, record):
            return 1.0 / (4.0 * math.pi)

    m = DNStyle()
    assert isinstance(m, DarkNewsDecay)
    assert isinstance(m, siren.interactions.Decay)
