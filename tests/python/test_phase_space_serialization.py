"""Raw-binding pickle coverage for phase-space channels and mixtures."""

import pickle

import pytest


siren = pytest.importorskip("siren")


def _box():
    placement = siren.geometry.Placement(siren.math.Vector3D(0.0, 0.0, 50.0))
    return siren.geometry.Box(placement, 2.0, 2.0, 2.0)


def _round_trip(value):
    payload = pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
    loaded = pickle.loads(payload)
    assert pickle.dumps(loaded, protocol=pickle.HIGHEST_PROTOCOL) == payload
    return loaded


def test_concrete_phase_space_channels_pickle():
    target = _box()
    inner = siren.injection.MultiChannelPhaseSpace(
        [siren.injection.Isotropic2BodyChannel(1)], [1.0]
    )
    nested = siren.injection.NestedMixtureChannel(inner)
    nested.label = "directed group"

    decay = siren.interactions.CharmMesonDecay(
        siren.dataclasses.ParticleType.D0
    )
    cross_section = siren.interactions.DummyCrossSection()

    channels = [
        siren.injection.Isotropic2BodyChannel(1),
        siren.injection.DetectorDirected2BodyChannel(
            target, 0, siren.injection.DirectedMode.Volume, 8.0
        ),
        siren.injection.DetectorDirectedAngularSectorChannel(
            target, 0.8, 1.0, 0.0, 1.0, 0
        ),
        siren.injection.DetectorDirected3BodyChannel(
            factorization=siren.injection.ThreeBodyMode.Direct,
            target=target,
            directed_index=0,
            mass_mode=siren.injection.InvariantMassMode.PowerLaw,
            power_law_nu=0.7,
            power_law_offset=0.05,
            volume=8.0,
        ),
        siren.injection.DetectorDirectedScatteringChannel(
            target,
            directed_index=1,
            variable=siren.injection.ScatteringVariable.Q2,
            mode=siren.injection.DirectedMode.Volume,
            q2_mode=siren.injection.ScatteringQ2Mode.Propagator,
            mediator_mass=0.2,
            volume=8.0,
        ),
        siren.injection.PhysicalDecayChannel(decay),
        siren.injection.PhysicalCrossSectionChannel(cross_section),
        nested,
    ]

    for channel in channels:
        loaded = _round_trip(channel)
        assert loaded.Name() == channel.Name()
        assert loaded.Topology() == channel.Topology()
        assert loaded.Measure() == channel.Measure()

    loaded_nested = _round_trip(nested)
    assert loaded_nested.label == "directed group"
    assert loaded_nested.mixture.weights == [1.0]
    assert loaded_nested.mixture.channels[0].Name() == "Isotropic2Body"


def test_multi_channel_phase_space_pickles_optimizer_state():
    mixture = siren.injection.MultiChannelPhaseSpace(
        [
            siren.injection.Isotropic2BodyChannel(0),
            siren.injection.DetectorDirected2BodyChannel(
                _box(), 0, siren.injection.DirectedMode.Volume, 8.0
            ),
        ],
        [0.25, 0.75],
    )
    mixture.kp_accumulator = [1.5, 2.5]
    mixture.kp_count = 7
    mixture.kp_succ_select = [3.0, 4.0]
    mixture.kp_fail_select = [0.5, 1.5]

    loaded = _round_trip(mixture)

    assert loaded.weights == [0.25, 0.75]
    assert loaded.kp_accumulator == [1.5, 2.5]
    assert loaded.kp_count == 7
    assert loaded.kp_succ_select == [3.0, 4.0]
    assert loaded.kp_fail_select == [0.5, 1.5]
    assert [channel.Name() for channel in loaded.channels] == [
        "Isotropic2Body",
        "DetectorDirected2Body",
    ]
    assert loaded.CommonTopology() == siren.injection.PhaseSpaceTopology.Decay2Body
