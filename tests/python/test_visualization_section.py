"""Filled sections, display presets and the side panel's state."""
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

from siren import _visualization_presets as presets
from siren import _visualization_section as section
from siren._visualization_render import SceneRenderer, feature_lines
from siren._visualization_scene import mesh_arrays, section_holes

pg = pytest.importorskip('pyg4ometry')
vtk = pytest.importorskip('vtk')


def rectangle(a, b, centre=(0, 0)):
    x, y = centre
    p = np.array([(x - a/2, y - b/2), (x + a/2, y - b/2), (x + a/2, y + b/2), (x - a/2, y + b/2)])
    return np.stack((p, np.roll(p, -1, axis=0)), axis=1)


def solid(name, *args):
    reg = pg.geant4.Registry()
    return mesh_arrays(getattr(pg.geant4.solid, name)('s', *args, reg).mesh())


def scene_of(*rows):
    """rows: (name, material, density, matrix, parent)."""
    return dict(prototypes={n: dict(name=n, material=m, density=d, role='') for n, m, d, _, _ in rows},
                instances=[dict(prototype=n, name='world/' + n, matrix=np.asarray(x).tolist(), parent=p)
                           for n, _, _, x, p in rows])


def shifted(x=0, y=0, z=0):
    m = np.eye(4)
    m[:3, 3] = x, y, z
    return m


def viewer():
    return SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, cutterNormals={},
                           cutterColors={}, bClipper=False)


def test_fill_removes_touching_and_overlapping_holes():
    outer = rectangle(10, 6)
    assert section.area(section.fill(outer)) == pytest.approx(60)
    assert section.area(section.fill(outer, [rectangle(2, 2)])) == pytest.approx(56)
    touching = [rectangle(5, 2, (-2.5, 0)), rectangle(5, 2, (2.5, 0))]
    assert section.area(section.fill(outer, touching)) == pytest.approx(40)
    overlapping = [rectangle(6, 2, (-2, 0)), rectangle(6, 2, (2, 0))]
    assert section.area(section.fill(outer, overlapping)) == pytest.approx(40)
    assert len(section.fill(outer, [rectangle(10, 6)])) == 0
    assert section.area(section.fill(outer, [rectangle(1, 1, (50, 0))])) == pytest.approx(60)
    assert len(section.fill(np.empty((0, 2, 2)))) == 0


def test_fill_does_not_need_shared_vertices():
    # A long side next to two short ones, as a Boolean mesh can leave behind.
    outline = np.array([[(0, 0), (4, 0)], [(4, 0), (4, 1)], [(4, 1), (4, 3)], [(4, 3), (0, 3)], [(0, 3), (0, 0)]], float)
    assert section.area(section.fill(outline[[3, 0, 4, 2, 1]])) == pytest.approx(12)


@pytest.mark.parametrize('height,expected', [(0, 8), (2.999, 8), (3, 8), (-3, 0), (3.001, 0)])
def test_cut_through_faces_and_vertices_stays_closed(height, expected):
    # At z = 3 and z = -3 the plane holds a whole face of the box.
    meshes = {'box': solid('Box', 2, 4, 6)}
    scene = scene_of(('box', 'Steel', 7.8, np.eye(4), -1))
    cut = section.SceneSections(scene['instances'], meshes).cut((0, 0, height), (0, 0, 1))
    assert sum(section.area(q) for q in cut.pieces.values()) == pytest.approx(expected)


def test_section_follows_rotated_scaled_and_reflected_placements():
    m = np.eye(4)
    m[:3, :3] = [[0, -2, 0], [1, 0, 0], [0, 0, 3]]
    mirrored = np.diag([-1., 1., 1., 1.])
    scene = scene_of(('a', 'Steel', 7.8, m, -1), ('b', 'Steel', 7.8, mirrored @ shifted(30), -1))
    box = solid('Box', 2, 4, 6)
    cut = section.SceneSections(scene['instances'], dict(a=box, b=box)).cut((0, 0, 1), (0, 0, 2))
    assert section.area(cut.pieces[0]) == pytest.approx(8 * 2)
    assert section.area(cut.pieces[1]) == pytest.approx(8)
    assert cut.find((-30, 0, 1)) == 1 and cut.find((0, 0, 1)) == 0 and cut.find((15, 0, 1)) is None


def test_round_section_area_and_tilted_plane():
    scene = scene_of(('tube', 'Steel', 7.8, np.eye(4), -1))
    reg = pg.geant4.Registry()
    tube = mesh_arrays(pg.geant4.solid.Tubs('t', 3, 5, 20, 0, 2*np.pi, reg, nslice=200).mesh())
    cut = section.SceneSections(scene['instances'], dict(tube=tube)).cut((0, 0, 0), (0, 1, 1))
    expected = np.pi * (25 - 9) * np.sqrt(2)
    assert section.area(cut.pieces[0]) == pytest.approx(expected, rel=1e-3)


def test_daughters_are_taken_out_and_levels_decide_overlaps():
    rows = [('outer', 'Concrete', 2.3, np.eye(4), -1), ('inner', 'Steel', 7.8, np.eye(4), 0)]
    meshes = dict(outer=solid('Box', 10, 10, 10), inner=solid('Box', 4, 4, 4))
    nested = scene_of(*rows)
    cut = section.SceneSections(nested['instances'], meshes).cut((0, 0, 0), (0, 0, 1))
    assert section.area(cut.pieces[0]) == pytest.approx(84) and section.area(cut.pieces[1]) == pytest.approx(16)
    # Not drawing the daughter leaves its opening.
    cut = section.SceneSections(nested['instances'], meshes).cut((0, 0, 0), (0, 0, 1), drawn={'outer'})
    assert set(cut.pieces) == {0} and section.area(cut.pieces[0]) == pytest.approx(84)
    # The overlapping sectors of a model: the higher level is the one there.
    flat = scene_of(rows[0], ('inner', 'Steel', 7.8, np.eye(4), -1))
    for levels, areas in (((1, 2), (84, 16)), ((2, 1), (100, 0))):
        for row, level in zip(flat['instances'], levels):
            row['level'] = level
        cut = section.SceneSections(flat['instances'], meshes).cut((0, 0, 0), (0, 0, 1))
        assert [section.area(cut.pieces.get(i, np.empty((0, 4, 2)))) for i in (0, 1)] == pytest.approx(areas)


def test_a_gap_in_one_strip_does_not_spread():
    # Earth-sized outlines round to slightly different ends: here the left side
    # starts 2e-4 above the right one, leaving a strip with a single crossing.
    outline = np.array([[(0, 2e-4), (0, 5)], [(4, 0), (4, 5)], [(0, 5), (4, 5)], [(4, 0), (0, 2e-4)]], float)
    outline[3] = [(4, 0), (4, 0)]  # the bottom is lost
    assert section.area(section.fill(outline)) == pytest.approx(20, abs=.01)


def test_exact_shapes_replace_or_remove_instances():
    scene = scene_of(('box', 'Steel', 7.8, shifted(5), -1), ('shell', 'Steel', 7.8, np.eye(4), -1))
    meshes = dict(box=solid('Box', 2, 2, 2), shell=solid('Box', 20, 20, 20))
    for row, level in zip(scene['instances'], (1, 2)):
        row['level'] = level
    cut = section.SceneSections(scene['instances'], meshes).cut((0, 0, 0), (0, 0, 1))
    assert set(cut.pieces) == {1}  # the shell's box covers everything
    left_out = section.SceneSections(scene['instances'], meshes, {1: lambda: None}).cut((0, 0, 0), (0, 0, 1))
    assert section.area(left_out.pieces[0]) == pytest.approx(4) and 1 not in left_out.pieces
    exact = {1: lambda: solid('Box', 2, 2, 2)}
    replaced = section.SceneSections(scene['instances'], meshes, exact).cut((0, 0, 0), (0, 0, 1))
    assert [section.area(replaced.pieces[i]) for i in (0, 1)] == pytest.approx([4, 4])
    # A hollow shell as a Boolean: its inside is open again.
    shell = {1: lambda: ('subtraction', solid('Box', 20, 20, 20), solid('Box', 18, 18, 18))}
    hollow = section.SceneSections(scene['instances'], meshes, shell).cut((0, 0, 0), (0, 0, 1))
    assert [section.area(hollow.pieces[i]) for i in (0, 1)] == pytest.approx([4, 76])


def test_boolean_cuts_combine_in_the_plane():
    a, b = rectangle(4, 4), rectangle(4, 4, (2, 2))
    for operation, expected in (('union', 28), ('intersection', 4), ('subtraction', 12)):
        quads = section.combine((operation, a, b))
        assert section.area(quads) == pytest.approx(expected)
        assert section.area(section.fill(section.sides(quads))) == pytest.approx(expected)
    flush = ('subtraction', ('union', a, rectangle(4, 4, (4, 0))), rectangle(8, 2, (2, -1)))
    assert section.area(section.combine(flush)) == pytest.approx(16)
    nothing = np.empty((0, 2, 2))
    assert len(section.combine(('intersection', a, nothing))) == 0
    assert section.area(section.combine(('subtraction', a, nothing))) == pytest.approx(16)


@pytest.mark.parametrize('detector,half', [('IceCube/IceCube-v1', 7e6), ('MINERvA/MINERvA-v1', 4.),
                                           ('ND280/ND280-v1', 6.), ('CCM/CCM-v1', 15.)])
def test_cuts_through_models_agree_with_siren(detectors_dir, tmp_path, detector, half):
    """Earth shells, sectors the GDML holds as boxes, and overlapping levels."""
    from siren import _visualization_scene as scenes, visualization
    from siren._visualization_view import ViewSession
    from siren.detector import DetectorModel, GeometryPosition
    from siren.math import Vector3D
    pytest.importorskip('pyvista')
    model = DetectorModel()
    model.LoadMaterialModel(str(detectors_dir / detector / 'materials.dat'))
    model.LoadDetectorModel(str(detectors_dir / detector / 'densities.dat'))
    session = ViewSession.__new__(ViewSession)
    session.model, session.vis = model, visualization
    stats = visualization.to_gdml(model, str(tmp_path / 'm.gdml'))
    session.boxed = set(stats['bbox_sectors'])
    events = []
    scenes.prepare_scene(tmp_path / 'm.gdml', tmp_path / 'out', overlapping=True, cache=False,
                         emit=events.append)
    scene = events[0]['scene']
    session._set_levels(scene)
    meshes = {e['prototype']: scenes.read_mesh(e['path']) for e in events if e['kind'] == 'mesh'}
    cutter = section.SceneSections(scene['instances'], meshes, session._exact_shapes(scene['instances']))
    origin = model.DetectorOrigin.get()
    centre = np.array([origin.GetX(), origin.GetY(), origin.GetZ()])
    names = {s.level: s.name for s in model.Sectors}
    rng = np.random.default_rng(1)
    for normal in ((1, 0, 0), (0, 0, 1), (.3, .5, .8)):
        offset = half * np.array([.0007, .0011, .0013])  # off any exact face
        cut = cutter.cut((centre + offset) * 1000, normal, drawn=set(scene['selected']))
        differ = 0
        points = cut.to_world((rng.random((300, 2)) - .5) * 2000 * half)
        for p in points:
            place = model.GeoPositionToDetPosition(GeometryPosition(Vector3D(*(p / 1000))))
            sector = model.GetContainingSector(place)
            want = sector.name if visualization._sector_density(sector) > .05 else None
            index = cut.find(p)
            differ += want != (None if index is None else names[scene['instances'][index]['level']])
        assert differ <= 6, normal  # facets of round surfaces


def test_undrawn_daughters_of_solids_become_openings():
    scene = scene_of(('world', 'Air', .001, np.eye(4), -1), ('wall', 'Concrete', 2.3, np.eye(4), 0),
                     ('duct', 'Air', .001, np.eye(4), 1), ('loose', 'Air', .001, np.eye(4), 0))
    assert section_holes(scene, ['wall']) == ['duct']
    assert section_holes(scene, ['wall', 'duct']) == []
    assert section_holes(scene, ['wall'], overlapping=True) == ['duct', 'loose']


def test_caps_clipping_layers_and_picking():
    scene = scene_of(('outer', 'Concrete', 2.3, np.eye(4), -1), ('inner', 'Steel', 7.8, np.eye(4), 0))
    renderer = SceneRenderer(viewer(), presets=dict(layers=[dict(label='Core', match='inner')]))
    renderer.set_structure(scene)
    renderer.add_mesh('outer', *solid('Box', 10, 10, 10))
    renderer.add_mesh('inner', *solid('Box', 4, 4, 4))
    renderer.finish_stage()
    renderer.set_section((0, 0, 1), (0, 0, 2))
    # vtkPlane keeps its positive side, so it points away from the viewer.
    assert renderer.plane.GetNormal() == (0, 0, -1)
    assert all(a.GetMapper().GetNumberOfClippingPlanes() == 1 for a in renderer.actor_instances)
    assert len(renderer.cap_actors) == 2
    hit = renderer.pick_ray([1, 1, 50], [1, 1, -50])
    assert hit['name'] == 'world/inner' and hit['position'] == pytest.approx([1, 1, 1]) and hit['section']
    assert renderer.pick_ray([4, 4, 50], [4, 4, -50])['name'] == 'world/outer'
    # Refinement keeps the section on the new actors.
    renderer.add_mesh('outer', *solid('Box', 10, 10, 10))
    renderer.finish_stage()
    assert all(a.GetMapper().GetNumberOfClippingPlanes() == 1 for a in renderer.actor_instances)
    renderer.layers['Core'] = False
    renderer.refresh_overlays()
    assert [a.GetVisibility() for a in renderer.cap_actors if renderer.actor_layer[a] == 'Core'] == [0]
    renderer.refresh_overlays(frozenset(['Concrete']))
    assert not any(a.GetVisibility() for a in renderer.cap_actors)
    renderer.clear_section()
    assert not renderer.cap_actors
    assert all(a.GetMapper().GetNumberOfClippingPlanes() == 0 for a in renderer.actor_instances)


def test_openings_are_cut_but_not_drawn():
    scene = scene_of(('wall', 'Concrete', 2.3, np.eye(4), -1), ('duct', 'Air', .001, np.eye(4), 0))
    scene.update(selected=['wall'], holes=['duct'])
    renderer = SceneRenderer(viewer())
    renderer.set_structure(scene)
    renderer.add_mesh('wall', *solid('Box', 10, 10, 10))
    renderer.add_mesh('duct', *solid('Box', 4, 4, 20))
    renderer.finish_stage()
    renderer.set_section((0, 0, 0), (0, 0, 1))
    assert len(renderer.viewer.actors) == 1 and len(renderer.cap_actors) == 1
    assert section.area(renderer.section.pieces[0]) == pytest.approx(84)
    # Asking for it later draws it.
    renderer.set_structure(dict(scene, selected=['duct'], holes=[]))
    renderer.add_mesh('duct', *solid('Box', 4, 4, 20))
    assert len(renderer.viewer.actors) == 2


def test_one_volume_can_be_placed_in_two_layers():
    scene = scene_of(('box', 'Steel', 7.8, np.eye(4), -1))
    scene['instances'] = [dict(scene['instances'][0], name='world/' + n, matrix=shifted(x).tolist())
                          for n, x in (('left', -10), ('right0x1a2b', 10))]
    renderer = SceneRenderer(viewer(), presets=dict(layers=[dict(label='Right', match='right', shown=False)]))
    renderer.set_structure(scene)
    renderer.add_mesh('box', *solid('Box', 2, 4, 6))
    assert {renderer.actor_layer[a]: a.GetVisibility() for a in renderer.actor_instances} == {
        presets.OTHER: 1, 'Right': 0}
    assert renderer.pick_ray([10, 0, 50], [10, 0, -50]) is None


def test_daughters_inherit_their_mothers_layer_but_not_the_worlds():
    scene = scene_of(('world', 'Air', .001, np.eye(4), -1), ('hall', 'Concrete', 2.3, np.eye(4), 0),
                     ('bolt', 'Steel', 7.8, np.eye(4), 1), ('loose', 'Steel', 7.8, np.eye(4), 0))
    layers = presets.load_presets(dict(layers=[dict(label='All', match='world'),
                                               dict(label='Hall', match='ha*')]))['layers']
    assert presets.assign_layers(scene, layers) == ['All', 'Hall', 'Hall', presets.OTHER]


def test_feature_lines_skip_facets_of_round_surfaces():
    assert len(feature_lines(*solid('Box', 2, 4, 6))) == 12
    reg = pg.geant4.Registry()
    tube = mesh_arrays(pg.geant4.solid.Tubs('t', 0, 5, 20, 0, 2*np.pi, reg, nslice=48).mesh())
    assert len(feature_lines(*tube)) == 96  # the two rims, none along the side


def test_materials_keep_density_colours_unless_a_preset_names_them():
    from siren.visualization import _material_vis_options
    registry = SimpleNamespace(materialDict={n: SimpleNamespace(density=d) for n, d in
                                             [('AIR', .001), ('MANTLE', 4.), ('ROCK', 2.65), ('CORE', 12.)]})
    plain = _material_vis_options(registry, pg.visualisation.VisualisationOptions)
    assert len({tuple(v.colour) for v in plain.values()}) == 4
    assert [plain[n].alpha for n in ('AIR', 'ROCK')] == [.12, .9]
    rows = presets.load_presets(dict(materials=[dict(label='Rock', match='RO*', colour='#ff0000')]))['materials']
    named = _material_vis_options(registry, pg.visualisation.VisualisationOptions, rows, 1.)
    assert named['ROCK'].colour == [1, 0, 0] and named['ROCK'].alpha == 1
    assert named['MANTLE'].colour == plain['MANTLE'].colour


def test_presets_override_colours(tmp_path):
    rows = dict(title='T', materials=[dict(label='Mine', match=['*Fe'], colour='#102030')],
                sections=[dict(label='Middle', origin=[0, 0, 1], normal=[0, 0, 1])])
    path = tmp_path / 'p.json'
    path.write_text(json.dumps(rows))
    loaded = presets.load_presets(path)
    row = presets.material_row('G4_Fe0xabc', loaded['materials'])
    assert row['label'] == 'Mine' and row['colour'] == pytest.approx((16/255, 32/255, 48/255))
    assert presets.material_row('G4_Al', loaded['materials']) is None
    assert presets.load_presets(loaded) == loaded  # already normalized


@pytest.mark.parametrize('bad,message', [
    (dict(colour='red'), 'unknown preset'), (dict(up=[0, 0, 0]), 'nonzero'), (dict(opacity=0), 'opacity'),
    (dict(sections=[dict(origin=[0, 0, 0])]), 'normal'),
    (dict(sections=[dict(origin=[0, 0, 0], normal=[0, 0, 0])]), 'nonzero'),
    (dict(materials=[dict(label='a', match='*', colour='red')]), 'rrggbb'),
    (dict(layers=[dict(label='a', match='*'), dict(label='a', match='b')]), 'different labels')])
def test_bad_presets_are_refused(bad, message):
    with pytest.raises(ValueError, match=message):
        presets.load_presets(bad)


native_window = pytest.mark.skipif(os.environ.get('SIREN_TEST_VTK_RENDER') != '1',
                                   reason='set SIREN_TEST_VTK_RENDER=1 with a native display')


@native_window
def test_native_panel_clicks_cut_hide_and_restore(tmp_path):
    from siren._visualization_scene import prepare_scene
    from siren._visualization_view import ViewSession
    from test_visualization_scene import fixture_gdml
    rows = dict(title='Fixture', layers=[dict(label='Turned', match='rotated')],
                sections=[dict(label='Middle', origin=[0, 0, 0], normal=[0, 0, 1])])
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, presets=rows,
                          cache_dir=tmp_path / 'cache')
    try:
        prepare_scene(session.path, tmp_path / 'out', emit=session.receive, **session.options)
        iren, panel = session.viewer.iren, session.panel

        def click(label):
            actor = next(a for a in panel._props if a.IsA('vtkTextActor') and a.GetInput() == label)
            x, y = actor.GetPosition()
            iren.SetEventInformation(int(x) + 2, int(y), 0, 0)
            iren.LeftButtonPressEvent()
            iren.LeftButtonReleaseEvent()

        click('Middle')
        assert session.current_section == 'Middle' and len(session.renderer.section.pieces) == 2
        click('Turned')
        assert [a.GetVisibility() for a in session.renderer.actor_instances] .count(0) == 1
        click('Steel')
        assert not any(a.GetVisibility() for a in session.renderer.cap_actors)
        click('Steel')
        click('Turned')
        click('x')
        slider = next(t for t in panel.targets if getattr(t[4], 'draggable', False))
        iren.SetEventInformation(int(slider[2]) - 8, int((slider[1] + slider[3]) / 2), 0, 0)
        iren.LeftButtonPressEvent()
        iren.LeftButtonReleaseEvent()
        low, high = session.slider_range()
        assert session.slider_value > .5 * (low + high) and session.current_section == ''
        before = list(session.renderer.section_normal)
        click('Other side')
        assert session.renderer.section_normal == [-c for c in before]
        camera = session.viewer.ren.GetActiveCamera().GetPosition()
        iren.SetEventInformation(10, 10, 0, 0)  # an empty part of the panel
        iren.LeftButtonPressEvent()
        iren.SetEventInformation(60, 60, 0, 0)
        iren.MouseMoveEvent()
        iren.LeftButtonReleaseEvent()
        assert session.viewer.ren.GetActiveCamera().GetPosition() == camera
        click('Whole model')
        assert session.renderer.plane is None and not session.renderer.cap_actors
    finally:
        session.viewer.renWin.Finalize()
