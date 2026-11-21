"""Prototype/transform checks plus opt-in real native-window acceptance."""
import os
import time
from types import SimpleNamespace

import numpy as np
import pytest

from siren._visualization_render import SceneRenderer, glyph_components
from siren._visualization_scene import mesh_arrays, normalize_regions
from test_visualization_scene import fixture_gdml, REGIONS

pg = pytest.importorskip('pyg4ometry')
vtk = pytest.importorskip('vtk')


def renderer_fixture(instancing=True):
    viewer = SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, bClipper=False)
    renderer = SceneRenderer(viewer, instancing=instancing)
    transforms = []
    for x in (-10, 10):
        m = np.eye(4)
        m[:3, :3] = [[0, -2, 0], [1, 0, 0], [0, 0, 3]]
        m[0, 3] = x
        transforms.append(m.tolist())
    # A sheared placement must retain the full matrix rather than approximate
    # it as rotation plus per-axis scaling.
    m = np.eye(4)
    m[0, 1] = .3
    transforms.append(m.tolist())
    metadata = dict(prototypes={'box':dict(name='box', material='Steel', density=7.8, role='')},
                    instances=[dict(prototype='box', name='copy%d' % i, matrix=m)
                               for i, m in enumerate(transforms)])
    renderer.set_structure(metadata)
    reg = pg.geant4.Registry()
    renderer.add_mesh('box', *mesh_arrays(pg.geant4.solid.Box('b', 2, 4, 6, reg).mesh()))
    return renderer, transforms


def test_glyph_quaternion_scale_and_exact_matrix_fallback():
    renderer, transforms = renderer_fixture()
    assert renderer.stats == dict(prototype_faces=12, placed_faces=36,
                                 glyph_instances=2, matrix_instances=1)
    glyph = next(a for a in renderer.viewer.actors.values() if a.GetMapper().IsA('vtkGlyph3DMapper'))
    data = glyph.GetMapper().GetInput()
    assert data.GetPoints().GetDataType() == vtk.VTK_DOUBLE
    q = data.GetPointData().GetArray('rotation').GetTuple(0)
    scale = data.GetPointData().GetArray('scale').GetTuple(0)
    rotation = [[0.] * 3 for _ in range(3)]
    vtk.vtkMath.QuaternionToMatrix3x3(q, rotation)
    assert np.asarray(rotation) @ np.diag(scale) == pytest.approx(np.asarray(transforms[0])[:3, :3])
    exact = next(a for a in renderer.viewer.actors.values() if not a.GetMapper().IsA('vtkGlyph3DMapper'))
    assert [[exact.GetMatrix().GetElement(i, j) for j in range(4)] for i in range(4)] == pytest.approx(np.asarray(transforms[2]))
    reflected = np.diag([-1., 1., 1., 1.])
    assert glyph_components(reflected) is None


def test_refinement_replaces_actors_and_keeps_material_groups():
    renderer, _ = renderer_fixture()
    old = list(renderer.viewer.actors.values())
    renderer.add_mesh('box', *(a.copy() for a in renderer.meshes['box']))
    assert len(renderer.viewer.actors) == 2
    assert len(renderer.material_actors['Steel']) == 2
    assert not any(a in renderer.viewer.actors.values() for a in old)
    assert renderer.stats['glyph_instances'] == 2


def test_large_scene_statistics_do_not_rescan_after_each_mesh(monkeypatch):
    viewer = SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, bClipper=False)
    renderer = SceneRenderer(viewer)
    names = ['volume_%d' % i for i in range(200)]
    renderer.set_structure(dict(
        prototypes={n:dict(name=n, material='Steel', density=7.8, role='') for n in names},
        instances=[dict(prototype=n, name=n, matrix=np.eye(4).tolist()) for n in names]))
    reg = pg.geant4.Registry()
    arrays = mesh_arrays(pg.geant4.solid.Box('b', 2, 4, 6, reg).mesh())
    scans = []
    original = renderer.update_stats
    def scan():
        scans.append(len(renderer.meshes))
        original()
    monkeypatch.setattr(renderer, 'update_stats', scan)
    for name in names:
        renderer.add_mesh(name, *arrays)
    assert scans == []
    assert renderer.stats == dict(prototype_faces=2400, placed_faces=2400,
                                 glyph_instances=200, matrix_instances=0)
    assert scans == [200]
    assert renderer.stats['glyph_instances'] == 200
    assert scans == [200]
    renderer.add_mesh(names[0], *(a.copy() for a in arrays))
    assert renderer.stats['prototype_faces'] == 2400
    assert scans == [200, 200]


def test_solid_aliases_instance_without_losing_material_or_pick_identity():
    viewer = SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, bClipper=False)
    renderer = SceneRenderer(viewer)
    prototypes = {n:dict(name=n, material=mat, density=density, role='', mesh_key='solid:box')
                  for n, mat, density in [('left','Steel',7.8),('right','Steel',7.8),('other','Al',2.7)]}
    instances = []
    for n, x in [('left', -10), ('right', 10), ('other', 30)]:
        m = np.eye(4)
        m[0, 3] = x
        instances.append(dict(prototype=n, name='world/'+n, matrix=m.tolist()))
    metadata = dict(prototypes=prototypes, instances=instances, selected=['left', 'other'])
    renderer.set_structure(metadata)
    reg = pg.geant4.Registry()
    arrays = mesh_arrays(pg.geant4.solid.Box('b', 2, 4, 6, reg).mesh())
    renderer.add_mesh('left', *arrays)
    renderer.add_mesh('other', *arrays)
    assert len(viewer.actors) == 2 and 'right' not in renderer.meshes
    assert renderer.pick_ray([10, 0, 10], [10, 0, -10]) is None
    # Deferred named alias joins the existing material's glyph actor.
    renderer.set_structure(dict(metadata, selected=['right']))
    renderer.add_mesh('right', *arrays)
    assert len(viewer.actors) == 2
    assert renderer.stats['glyph_instances'] == 3
    assert renderer.meshes['left'][0] is renderer.meshes['right'][0]
    for name, x, material in [('left',-10,'Steel'),('right',10,'Steel'),('other',30,'Al')]:
        hit = renderer.pick_ray([x, 0, 10], [x, 0, -10])
        assert hit['name'] == 'world/'+name and hit['material'] == material
        assert hit['position'] == pytest.approx([x, 0, 3])
    assert len(renderer.locators) == 1  # one geometry BVH across materials
    # Repeated delivery must not rebuild an already-instanced group.
    old = list(viewer.actors.values())
    renderer.add_mesh('right', *arrays)
    assert list(viewer.actors.values()) == old


def test_worker_burst_yields_and_keeps_mesh_files_until_consumed(monkeypatch):
    from collections import deque
    from siren import _visualization_view as view
    session = view.ViewSession.__new__(view.ViewSession)
    now, consumed, frames, closed = [0.], [], [], []
    events = [dict(kind='mesh', index=i) for i in range(100)] + [dict(kind='done')]
    def poll():
        batch = events[:]
        events.clear()
        return batch
    session.job = SimpleNamespace(finished=True, poll=poll, close=lambda:closed.append(True))
    session._events = deque()
    session.pending, session.loaded = set(), set()
    session.error = None
    session._closing = lambda:False
    def receive(e):
        assert not closed, 'worker mesh files were removed before the queue drained'
        consumed.append(e)
        now[0] += .02
    session.receive = receive
    session._render = lambda:frames.append(True)
    monkeypatch.setattr(view.time, 'perf_counter', lambda:now[0])
    session.tick()
    assert 1 <= len(consumed) <= 3
    assert session.job is not None and not closed and len(frames) == 1
    for _ in range(100):
        if session.job is None:
            break
        session.tick()
    assert [e['index'] for e in consumed if e['kind'] == 'mesh'] == list(range(100))
    assert consumed[-1]['kind'] == 'done'
    assert closed == [True] and session.job is None and not session._events


def test_direct_mesh_actor_is_pickable_without_expanding_prototypes():
    from siren._visualization_render import polydata
    renderer, _ = renderer_fixture()
    vertices = np.array([[49., -1., 5.], [51., -1., 5.], [50., 1., 5.]])
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata(vertices, np.array([[0, 1, 2]], dtype=np.int64)))
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer.external_actors['direct'] = actor
    assert renderer.pick_ray([50, 0, 10], [50, 0, -10])['position'] == pytest.approx([50, 0, 5])
    actor.SetVisibility(False)
    assert renderer.pick_ray([50, 0, 10], [50, 0, -10]) is None


def test_exterior_legend_matches_region_colours():
    viewer = SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, bClipper=False)
    renderer = SceneRenderer(viewer, display='exterior', regions=normalize_regions(REGIONS))
    renderer.set_structure(dict(prototypes={
        'front': dict(name='front', material='Glass', density=2.2, role='photocathode_inner_surface'),
        'rear': dict(name='rear', material='Glass', density=2.2, role='reflector_inner_surface'),
    }, instances=[]))
    assert 'Glass' not in viewer.legend_options
    assert viewer.legend_options['PMT photocathode'].colour == renderer.vis_option('front').colour
    assert viewer.legend_options['PMT reflector'].colour == renderer.vis_option('rear').colour


native_window = pytest.mark.skipif(os.environ.get('SIREN_TEST_VTK_RENDER') != '1',
                                   reason='set SIREN_TEST_VTK_RENDER=1 with a native display')


@native_window
@pytest.mark.parametrize('coloured', [True, False])
def test_native_glyph_pick_and_controls(tmp_path, coloured):
    from siren._visualization_view import ViewSession
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, legend=False,
                          interactive=False, coloured=coloured, cache_dir=tmp_path / 'cache')
    from siren._visualization_scene import prepare_scene
    try:
        prepare_scene(session.path, tmp_path / 'out', emit=session.receive, **session.options)
        v = session.viewer
        v.ren.SetWorldPoint(-10, 0, 3, 1)
        v.ren.WorldToDisplay()
        x, y, _ = v.ren.GetDisplayPoint()
        picked = session.renderer.pick(round(x), round(y))
        assert picked is not None and picked['name'].endswith('/other')
        assert picked['material'] == 'Steel'
        assert picked['position'] == pytest.approx([-10, 0, 3], abs=.1)
        # Change opacity, refine, then restore: replacement actors retain state.
        v.iren.SetKeyEventInformation(0, 0, 'n', 0, 'n')
        v.interactorStyle.InvokeEvent('KeyPressEvent')
        opacity = session.controls['opacity']
        for name, arrays in list(session.renderer.meshes.items()):
            session.renderer.add_mesh(name, *arrays)
        session.controls['refresh']()
        assert all(a.GetProperty().GetOpacity() == pytest.approx(opacity)
                   for a in session.renderer.material_actors['Steel'])
        v.iren.SetKeyEventInformation(0, 0, 'v', 0, 'v')
        v.interactorStyle.InvokeEvent('KeyPressEvent')
        assert all(a.GetProperty().GetOpacity() == pytest.approx(.9 if coloured else .5)
                   for a in session.renderer.material_actors['Steel'])
    finally:
        if session.job is not None:
            session.job.close()
        session.viewer.renWin.Finalize()


@native_window
@pytest.mark.parametrize('cutter,clipper', [(True, False), (False, True), (True, True)])
def test_native_sections_and_clipping_use_full_meshes(tmp_path, cutter, clipper):
    from siren._visualization_view import ViewSession
    from siren._visualization_scene import prepare_scene
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, legend=False,
                          interactive=False, cutter=cutter, clipper=clipper,
                          display='exterior', regions=REGIONS, cache_dir=tmp_path / 'cache')
    try:
        prepare_scene(session.path, tmp_path / 'out', emit=session.receive, **session.options)
        assert session.options['display'] == 'full'
        assert session.renderer.stats['glyph_instances'] == 0
        if cutter:
            assert {'xy', 'xz', 'yz'} <= session.viewer.cutters.keys()
            assert session.viewer.getCutterPolydata('xy').GetNumberOfPoints() > 0
        if clipper:
            assert session.viewer.clippers and session.controls['clipper_widget'] is not None
            for clip in session.viewer.clippers:
                clip.Update()
                # ClippedOutput retains unused source points; inspect cells.
                output = clip.GetClippedOutput()
                assert output.GetNumberOfCells() > 0
                assert all(output.GetCell(i).GetBounds()[1] <= 1e-8
                           for i in range(output.GetNumberOfCells()))
            session.viewer.setClipper([0, 0, 0], [0, 1, 0])
            assert tuple(session.viewer.clippers[0].GetClipFunction().GetNormal()) == (0, 1, 0)
    finally:
        session.viewer.renWin.Finalize()


@native_window
def test_native_progress_refinement_deferred_loading_and_camera(tmp_path):
    from siren._visualization_view import ViewSession
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, legend=False,
                          cache_dir=tmp_path / 'cache')
    receive = session.receive
    witness = {}
    def event(e):
        receive(e)
        if e['kind'] == 'stage_done' and e['stage'] == 'preview':
            session.viewer.iren.SetKeyEventInformation(0, 0, ' ', 0, 'Left')
            session.viewer.interactorStyle.InvokeEvent('KeyPressEvent')
            assert session.camera_dirty
            witness['camera'] = session.viewer.ren.GetActiveCamera().GetPosition()
    session.receive = event
    # Each worker process pays the pyg4ometry import (several seconds on some
    # installations), so the budget is per phase rather than per test.
    deadline = [time.monotonic() + 20]
    def stop(*args):
        if time.monotonic() > deadline[0]:
            witness['timeout'] = True
            session.viewer.iren.TerminateApp()
        elif session.detail_ready and 'gas' not in witness:
            witness['detail_camera'] = session.viewer.ren.GetActiveCamera().GetPosition()
            session.controls['gas_visible'] = True
            session.request_geometry('gas')
            witness['gas'] = True
            deadline[0] = time.monotonic() + 20
        elif session.detail_ready and session.job is None and {'world', 'parent'} <= session.loaded:
            session.viewer.iren.TerminateApp()
    session.viewer.iren.AddObserver('TimerEvent', stop)
    session.run()
    assert 'timeout' not in witness
    assert witness['camera'] == pytest.approx(witness['detail_camera'])
    assert {'world', 'parent', 'cube'} <= session.loaded
    assert session.metrics['window_frame_seconds'] < session.metrics['preview_frame_seconds']
    assert session.metrics['preview_frame_seconds'] <= session.metrics['detail_frame_seconds']
    assert session.job is None


@native_window
def test_native_early_close_cancels_worker(tmp_path):
    from siren._visualization_view import ViewSession
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, legend=False,
                          cache_dir=tmp_path / 'cache')
    jobs = []
    def close(*args):
        if session.job is not None:
            jobs.append((session.job.process, session.job.directory.name))
        session.viewer.iren.TerminateApp()
    session.viewer.iren.AddObserver('TimerEvent', close)
    session.run()
    assert jobs and session.metrics['cancelled']
    assert all(p.poll() is not None and not os.path.exists(d) for p, d in jobs)


@native_window
def test_native_siren_model_and_temporary_export_cache(tmp_path, detectors_dir):
    from siren.detector import DetectorModel
    from siren.visualization import view
    model = DetectorModel()
    base = detectors_dir / 'CCM' / 'CCM-v1'
    model.LoadMaterialModel(str(base / 'materials.dat'))
    model.LoadDetectorModel(str(base / 'densities.dat'))
    first, second = {}, {}
    for stats in (first, second):
        result = view(model, screenshot=tmp_path / 'model.png', axes=False,
                      legend=False, cache_dir=tmp_path / 'cache', timings=stats)
        assert result == tmp_path / 'model.png'
        assert (tmp_path / 'model.png').stat().st_size > 1000
        assert stats['prototype_faces'] > 0
    assert first['stages'][0]['cache_misses'] > 0
    assert second['stages'][0]['cache_misses'] == 0
    assert second['stages'][0]['cache_hits'] == first['stages'][0]['cache_misses']


def test_direct_mesh_actors_do_not_need_pyvista_without_mesh_sectors(monkeypatch):
    import sys
    from siren.visualization import _add_mesh_actors
    monkeypatch.setitem(sys.modules, 'pyvista', None)  # import raises ImportError
    model = SimpleNamespace(Sectors=[SimpleNamespace(geo=SimpleNamespace(), material_id=0)])
    assert _add_mesh_actors(vtk, vtk.vtkRenderer(), {}, model, None, SimpleNamespace(materialDict={})) == 0


@native_window
def test_native_direct_mesh_failure_keeps_viewer(tmp_path, detectors_dir, monkeypatch):
    from siren.detector import DetectorModel
    from siren import visualization
    model = DetectorModel()
    base = detectors_dir / 'CCM' / 'CCM-v1'
    model.LoadMaterialModel(str(base / 'materials.dat'))
    model.LoadDetectorModel(str(base / 'densities.dat'))
    def broken(*args, **kwargs):
        raise ImportError("No module named 'pyvista'")
    monkeypatch.setattr(visualization, '_add_mesh_actors', broken)
    stats = {}
    with pytest.warns(RuntimeWarning, match='pyvista'):
        result = visualization.view(model, screenshot=tmp_path / 'model.png', axes=False,
                                    legend=False, cache_dir=tmp_path / 'cache', timings=stats)
    assert result == tmp_path / 'model.png' and stats['prototype_faces'] > 0
    assert 'pyvista' in stats['direct_mesh_error']


@native_window
def test_native_plain_viewer_has_no_legend_and_point_only_picking(tmp_path):
    from siren._visualization_view import ViewSession
    from siren._visualization_scene import prepare_scene
    session = ViewSession(fixture_gdml(tmp_path / 'a.gdml'), axes=False, legend=True,
                          interactive=False, coloured=False, picker=False,
                          cache_dir=tmp_path / 'cache')
    try:
        prepare_scene(session.path, tmp_path / 'out', emit=session.receive, **session.options)
        assert session.controls['legend'] is None
        assert session.viewer.pick_scene is None
    finally:
        session.viewer.renWin.Finalize()


def test_pick_orders_hits_by_world_distance_across_scaled_placements():
    viewer = SimpleNamespace(ren=vtk.vtkRenderer(), actors={}, cutterOrigins={}, bClipper=False)
    renderer = SceneRenderer(viewer)
    near, far = np.eye(4), np.eye(4)
    near[2, 3] = 20                      # unit box centred at z=20
    far[:3, :3] = np.diag([1, 1, 10])    # stretched box centred at z=0, top at z=+30
    proto = {'box': dict(name='box', material='Steel', density=7.8, role='')}
    renderer.set_structure(dict(prototypes=proto, instances=[
        dict(prototype='box', name='far', matrix=far.tolist()),
        dict(prototype='box', name='near', matrix=near.tolist())]))
    reg = pg.geant4.Registry()
    renderer.add_mesh('box', *mesh_arrays(pg.geant4.solid.Box('b', 2, 2, 6, reg).mesh()))
    hit = renderer.pick_ray([0, 0, 100], [0, 0, -100])
    assert hit['name'] == 'far' and hit['position'] == pytest.approx([0, 0, 30])
    hit = renderer.pick_ray([0, 0, 25], [0, 0, -100])
    assert hit['name'] == 'near' and hit['position'] == pytest.approx([0, 0, 23])


def test_session_surfaces_worker_warnings(monkeypatch):
    from siren import _visualization_view as view
    session = view.ViewSession.__new__(view.ViewSession)
    session.metrics, session._closing = {}, lambda: False
    with pytest.warns(RuntimeWarning, match='cache disabled'):
        session.receive(dict(kind='warning', message='mesh cache disabled: test'))
    assert session.metrics['warnings'] == ['mesh cache disabled: test']


def test_refinement_releases_preview_polydata():
    renderer, _ = renderer_fixture()
    preview = renderer.meshes['box']
    assert len(renderer._source_data) == 1
    renderer.add_mesh('box', *(a.copy() for a in preview))
    assert len(renderer._source_data) == 1
    assert (id(preview[0]), id(preview[1])) not in renderer._source_data
    # Redelivering the same arrays is a no-op and keeps one live source.
    renderer.add_mesh('box', *renderer.meshes['box'])
    assert len(renderer._source_data) == 1 and renderer._source_users[next(iter(renderer._source_data))] == 1


def test_finish_stage_does_not_accumulate_pipeline_bookkeeping():
    from pyg4ometry.visualisation import VtkViewerNew
    viewer = VtkViewerNew(defaultCutters=True, axisCubeWidget=False)
    renderer = SceneRenderer(viewer)
    reg = pg.geant4.Registry()
    arrays = mesh_arrays(pg.geant4.solid.Box('b', 2, 4, 6, reg).mesh())
    renderer.set_structure(dict(prototypes={'box': dict(name='box', material='Steel', density=7.8, role='')},
                                instances=[dict(prototype='box', name='c%d' % i, matrix=np.eye(4).tolist()) for i in range(3)]))
    for _ in range(3):  # preview, detail, deferred load
        renderer.add_mesh('box', *(a.copy() for a in arrays))
        renderer.finish_stage()
    assert len(viewer.instanceNameDict) == 3 and len(viewer.polydata) == 1
    assert len(renderer.material_actors['Steel']) == 1
    viewer.renWin.Finalize()
