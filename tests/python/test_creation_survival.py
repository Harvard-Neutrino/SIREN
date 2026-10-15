"""Independent exponential survival checks on the native API."""

import math
import pickle

import pytest
import siren
from siren import dataclasses as dc, distributions as dist, injection as inj
from siren import interactions as xs, utilities
from siren.math import Vector3D
from test_decay_channels import channel, P

HBARC = 1.973e-16


def setup(tmp_path, depth, full, z0=0.0, branch=1.0, n=8, volume=False, selected=False):
    width = depth / 10.0 * 3.0 * HBARC
    model = channel([P.NuMu, P.Gamma], branch * width, width)
    models = [model]
    if selected:
        models.append(channel([P.NuE, P.Gamma], 0.3 * width, width))
    collection = xs.InteractionCollection(P.N4, models)
    if selected:
        collection.SetDecayChannels(model.GetPossibleSignatures())
    detector = siren.detector.DetectorModel()
    path = tmp_path / "source.csv"
    path.write_text(f"E,px,py,pz,x0,y0,z0,m\n{math.sqrt(10)},0,0,3,0,0,{z0},1\n")
    sphere = siren.geometry.Sphere(
        siren.geometry.Placement(Vector3D(0, 0, 11)), 1.0, 0.0
    )
    process = inj.PrimaryInjectionProcess(P.N4, collection)
    process.distributions = [
        dist.PrimaryExternalDistribution(str(path)),
        (dist.SphereVolumePositionDistribution(sphere) if volume else
         dist.PrimaryBoundedVertexDistribution(sphere, 100.0)),
    ]
    mode = (
        inj.VertexWeightingMode.PropagatedFromCreation()
        if full
        else inj.VertexWeightingMode.Propagated()
    )
    process.weighting_mode = mode
    physical = inj.PhysicalProcess(P.N4, collection)
    physical.weighting_mode = mode
    generator = inj._Injector(n, detector, process, utilities.SIREN_random(31))
    return generator, inj._Weighter([generator], detector, physical)


@pytest.mark.parametrize("depth", [0.5, 5.0, 50.0, 300.0, 700.0])
@pytest.mark.parametrize("branch", [1.0, 1e-6])
def test_upstream_survival_is_charged_exactly_once(tmp_path, depth, branch):
    generator, weighter = setup(tmp_path, depth, True, branch=branch)
    plain, conditional = setup(tmp_path, depth, False, branch=branch)
    events = [generator.GenerateEvent() for _ in range(8)]
    other = [plain.GenerateEvent() for _ in range(8)]
    assert all(len(e.tree) == 1 for e in events + other)
    inside = -math.expm1(-0.2 * depth)
    survival = math.exp(-depth)
    for a, b in zip(events, other):
        assert list(a.tree[0].record.interaction_vertex) == list(
            b.tree[0].record.interaction_vertex
        )
        expected = branch * inside / 8
        assert conditional.EventWeight(b) == pytest.approx(expected, rel=2e-12, abs=0)
        if depth == 700.0 and branch == 1e-6:
            with pytest.raises(RuntimeError, match="inverse weight overflow"):
                weighter.EventWeight(a)
            continue
        weight = weighter.EventWeight(a)
        assert weight == pytest.approx(expected * survival, rel=1e-11, abs=0)
        breakdown = weighter.EventWeightWithBreakdown(a)
        assert breakdown.total == weight
        assert breakdown.vertices[0].survival == pytest.approx(
            survival, rel=1e-11, abs=0
        )
        assert conditional.EventWeightWithBreakdown(b).vertices[0].survival == 1.0


def test_creation_inside_volume_has_no_upstream_survival(tmp_path):
    generator, weighter = setup(tmp_path, 20.0, True, z0=11.0)
    plain, conditional = setup(tmp_path, 20.0, False, z0=11.0)
    events = [generator.GenerateEvent() for _ in range(8)]
    other = [plain.GenerateEvent() for _ in range(8)]
    assert [weighter.EventWeight(e) for e in events] == [
        conditional.EventWeight(e) for e in other
    ]


@pytest.mark.parametrize("secondary", [False, True])
def test_each_process_uses_its_own_creation_point(secondary):
    model = channel([P.NuMu, P.Gamma], 3 * HBARC, 3 * HBARC)
    collection = xs.InteractionCollection(P.N4, [model])
    physical = inj.PhysicalProcess(P.N4, collection)
    physical.weighting_mode = inj.VertexWeightingMode.PropagatedFromCreation()
    process_type = (
        inj.SecondaryInjectionProcess if secondary else inj.PrimaryInjectionProcess
    )
    weighter_type = (
        inj.SecondaryProcessWeighter if secondary else inj.PrimaryProcessWeighter
    )
    process = process_type(P.N4, collection)
    process.weighting_mode = physical.weighting_mode
    weighter = weighter_type(physical, process, siren.detector.DetectorModel())
    record = dc.InteractionRecord()
    record.signature = model.GetPossibleSignatures()[0]
    record.primary_mass = 1.0
    record.primary_momentum = [math.sqrt(10), 0, 0, 3]
    record.secondary_masses = [0.0, 0.0]
    record.secondary_momenta = [[0.0] * 4 for _ in range(2)]
    record.secondary_helicities = [0.0, 0.0]
    record.primary_initial_position = [0, 0, 7]
    record.interaction_vertex = [0, 0, 10.3]
    bounds = (Vector3D(0, 0, 10), Vector3D(0, 0, 12))
    assert weighter.PhysicalProbability(bounds, record) == pytest.approx(
        math.exp(-3.3) / (4 * math.pi), rel=2e-12, abs=0
    )


def test_survival_underflow_remains_a_visible_weight_error(tmp_path):
    generator, weighter = setup(tmp_path, 800.0, True, n=1)
    event = generator.GenerateEvent()
    assert len(event.tree) == 1
    with pytest.raises(RuntimeError, match="CreationSurvival"):
        weighter.EventWeight(event)
    with pytest.raises(RuntimeError, match="CreationSurvival"):
        weighter.EventWeightWithBreakdown(event)


@pytest.mark.parametrize(
    "preset", [inj.VertexWeightingMode.Fixed, inj.VertexWeightingMode.ExternalBounds]
)
def test_invalid_mode_does_not_replace_process_state(preset):
    mode = preset()
    mode.survival_from_creation = True
    with pytest.raises(RuntimeError, match="survival_from_creation"):
        mode.Validate()
    process = inj.PhysicalProcess()
    with pytest.raises(RuntimeError, match="survival_from_creation"):
        process.weighting_mode = mode
    assert process.weighting_mode == inj.VertexWeightingMode.Propagated()


def test_native_archives_retain_opt_in_mode(tmp_path):
    from test_native_phase_space import record_state

    model = channel([P.NuMu, P.Gamma], 0.1e-16, 1e-16)
    collection = xs.InteractionCollection(P.N4, [model])
    detector = siren.detector.DetectorModel()
    process = inj.PrimaryInjectionProcess(P.N4, collection)
    process.weighting_mode = inj.VertexWeightingMode.PropagatedFromCreation()
    process.distributions = [
        dist.PrimaryMass(1.0),
        dist.Monoenergetic(2.0),
        dist.FixedDirection(Vector3D(0, 0, 1)),
        dist.PointSourcePositionDistribution(Vector3D(), 10.0),
    ]
    generator = inj._Injector(8, detector, process, utilities.SIREN_random(3))
    physical = inj.PhysicalProcess(P.N4, collection)
    physical.weighting_mode = process.weighting_mode
    weighter = inj._Weighter([generator], detector, physical)
    generator.GenerateEvent()
    path = str(tmp_path / "survival.injector")
    generator.SaveInjector(path)
    loaded = inj._Injector(8, path, utilities.SIREN_random(1))
    copied = pickle.loads(pickle.dumps(generator))
    weight_path = str(tmp_path / "survival")
    weighter.SaveWeighter(weight_path)
    restored = inj._Weighter([], weight_path)
    for _ in range(4):
        event = generator.GenerateEvent()
        assert (
            record_state(event)
            == record_state(loaded.GenerateEvent())
            == record_state(copied.GenerateEvent())
        )
        assert restored.EventWeight(event) == pytest.approx(
            weighter.EventWeight(event) * generator.InjectionAttempts(),
            rel=2e-12,
            abs=0,
        )
    assert loaded.GetPrimaryProcess().weighting_mode.survival_from_creation
    assert restored.GetPrimaryPhysicalProcess().weighting_mode.survival_from_creation


@pytest.mark.parametrize("depth", [0.5, 50.0])
@pytest.mark.parametrize("full", [False, True])
def test_selected_propagated_channel_uses_full_lifetime(tmp_path, depth, full):
    generator, weighter = setup(tmp_path, depth, full, branch=0.2, n=32, selected=True)
    events = [generator.GenerateEvent() for _ in range(32)]
    expected = 0.2 * -math.expm1(-0.2 * depth) / 32
    if full:
        expected *= math.exp(-depth)
    for event in events:
        assert list(event.tree[0].record.signature.secondary_types) == [P.NuMu, P.Gamma]
        assert weighter.EventWeight(event) == pytest.approx(expected, rel=1e-11, abs=0)
        assert weighter.EventWeightWithBreakdown(event).total == weighter.EventWeight(event)


@pytest.mark.parametrize("injection_full", [False, True])
@pytest.mark.parametrize("physical_full", [False, True])
def test_physical_process_controls_creation_survival(tmp_path, injection_full, physical_full):
    generator, weighter = setup(tmp_path, 5.0, True)
    mode = lambda full: (inj.VertexWeightingMode.PropagatedFromCreation() if full
                         else inj.VertexWeightingMode.Propagated())
    generator.GetPrimaryProcess().weighting_mode = mode(injection_full)
    weighter.GetPrimaryPhysicalProcess().weighting_mode = mode(physical_full)
    events = [generator.GenerateEvent() for _ in range(8)]
    expected = -math.expm1(-1.0) / 8 * (math.exp(-5.0) if physical_full else 1.0)
    for event in events:
        assert weighter.EventWeight(event) == pytest.approx(expected, rel=1e-11, abs=0)


def test_volume_sampler_replaces_creation_with_bounds_entry(tmp_path):
    generator, weighter = setup(tmp_path, 50.0, True, volume=True)
    plain, conditional = setup(tmp_path, 50.0, False, volume=True)
    events = [generator.GenerateEvent() for _ in range(8)]
    other = [plain.GenerateEvent() for _ in range(8)]
    for event, reference in zip(events, other):
        record = event.tree[0].record
        assert list(record.primary_initial_position) != [0.0, 0.0, 0.0]
        assert weighter.EventWeightWithBreakdown(event).vertices[0].survival == pytest.approx(1.0)
        assert weighter.EventWeight(event) == pytest.approx(conditional.EventWeight(reference), rel=1e-12, abs=0)
