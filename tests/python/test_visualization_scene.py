"""Independent geometry/cache witnesses for selective viewer preparation."""
import hashlib
import os
from pathlib import Path
import time

import numpy as np
import pytest

from siren import _visualization_scene as scene

pg = pytest.importorskip("pyg4ometry")

# Example surface-region mapping for display='exterior' (CCM's PMT shells).
# SIREN itself carries no detector-specific region names; callers supply them.
REGIONS = {
    "auxtype": "pmt_region",
    "styles": {
        "external_tpb": {"label": "PMT TPB", "colour": [.95, .94, .82], "alpha": 1.},
        "photocathode_inner_surface": {"label": "PMT photocathode", "colour": [.65, .42, .16], "alpha": 1.},
        "reflector_inner_surface": {"label": "PMT reflector", "colour": [.8, .82, .85], "alpha": 1.},
        "bare_transparent_glass": {"label": "PMT bare glass", "colour": [.72, .88, .94], "alpha": .25},
    },
    "hidden": ["vacuum"],
}


def fixture_gdml(path, body=None):
    structure = body or '''
      <volume name="cube"><materialref ref="Steel"/><solidref ref="cube_s"/></volume>
      <volume name="parent"><materialref ref="Air"/><solidref ref="parent_s"/>
       <physvol name="rotated"><volumeref ref="cube"/>
        <position x="10" y="0" z="0"/><rotation x="0" y="0" z="90" unit="deg"/>
       </physvol>
      </volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>
       <physvol name="translated"><volumeref ref="parent"/><position x="0" y="20" z="0"/></physvol>
       <physvol name="other"><volumeref ref="cube"/><position x="-10" y="0" z="0"/></physvol>
      </volume>'''
    path.write_text('''<?xml version="1.0"?>
    <gdml><define/><materials>
     <material name="Air" Z="1"><D value="0.001" unit="g/cm3"/><atom value="1"/></material>
     <material name="Steel" Z="1"><D value="7.8" unit="g/cm3"/><atom value="1"/></material>
    </materials><solids>
     <box name="cube_s" x="2" y="4" z="6" lunit="mm"/>
     <box name="parent_s" x="50" y="50" z="50" lunit="mm"/>
     <box name="world_s" x="100" y="100" z="100" lunit="mm"/>
    </solids><structure>''' + structure + '''</structure>
    <setup name="Default" version="1.0"><world ref="world"/></setup></gdml>''')
    return path


def prepare(path, tmp_path, **kwargs):
    events = []
    result, stats = scene.prepare_scene(path, tmp_path / 'out', emit=events.append,
                                        cache_dir=tmp_path / 'cache', **kwargs)
    meshes = {e['prototype']: scene.read_mesh(e['path']) for e in events if e['kind'] == 'mesh'}
    return result, stats, meshes


def volume(vertices, faces):
    triangles = vertices[faces]
    return np.einsum('ij,ij->i', triangles[:, 0],
                     np.cross(triangles[:, 1], triangles[:, 2])).sum() / 6


def test_deferred_parents_preserve_nested_transform(tmp_path, monkeypatch):
    calls = []
    original = pg.geant4.solid.Box.mesh
    def mesh(solid):
        calls.append(solid.name)
        return original(solid)
    monkeypatch.setattr(pg.geant4.solid.Box, 'mesh', mesh)
    s, stats, meshes = prepare(fixture_gdml(tmp_path / 'a.gdml'), tmp_path)
    assert calls == ['cube_s']  # Air parents are traversed without meshing.
    assert set(meshes) == {'cube'}
    assert stats['cache_misses'] == 1
    assert volume(*meshes['cube']) == pytest.approx(48)
    placements = [i for i in s['instances'] if i['prototype'] == 'cube']
    # GDML passive +90 degree rotation is an active -90 degree rotation.
    assert np.asarray(placements[0]['matrix']) @ [1, 0, 0, 1] == pytest.approx([10, 19, 0, 1])
    assert np.asarray(placements[1]['matrix']) @ [1, 0, 0, 1] == pytest.approx([-9, 0, 0, 1])


def test_logical_aliases_share_mesh_files_but_keep_names_and_selection(tmp_path, monkeypatch):
    path = fixture_gdml(tmp_path / 'aliases.gdml', '''
      <volume name="left"><materialref ref="Steel"/><solidref ref="cube_s"/></volume>
      <volume name="right"><materialref ref="Steel"/><solidref ref="cube_s"/></volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>
       <physvol name="first"><volumeref ref="left"/><position x="-10" y="0" z="0"/></physvol>
       <physvol name="second"><volumeref ref="right"/><position x="10" y="0" z="0"/></physvol>
      </volume>''')
    events = []
    metadata, stats = scene.prepare_scene(path, tmp_path/'out', cache_dir=tmp_path/'cache', emit=events.append)
    meshes = [e for e in events if e['kind'] == 'mesh']
    assert {e['prototype'] for e in meshes} == {'left', 'right'}
    assert len({e['path'] for e in meshes}) == 1
    assert metadata['prototypes']['left']['mesh_key'] == metadata['prototypes']['right']['mesh_key']
    assert stats['cache_misses'] == 1 and stats['shared_mesh_reuses'] == 1
    def forbidden(*args):
        raise AssertionError('shared cached solid must not be remeshed')
    monkeypatch.setattr(pg.geant4.solid.Box, 'mesh', forbidden)
    assert scene.cached_scene_available(path, cache_dir=tmp_path/'cache')
    events.clear()
    _, warm = scene.prepare_scene(path, tmp_path/'filtered', cache_dir=tmp_path/'cache',
                                  hidden_volumes=['right'], emit=events.append)
    assert warm['cache_hits'] == 1 and warm['cache_misses'] == 0
    assert [e['prototype'] for e in events if e['kind'] == 'mesh'] == ['left']


def test_warm_cache_never_calls_mesher_and_preserves_bits(tmp_path, monkeypatch):
    path = fixture_gdml(tmp_path / 'a.gdml')
    first, _, a = prepare(path, tmp_path)
    def forbidden(*args):
        raise AssertionError('warm cache must not mesh')
    monkeypatch.setattr(pg.geant4.solid.Box, 'mesh', forbidden)
    second, stats, b = prepare(path, tmp_path)
    assert first['cache_key'] == second['cache_key']
    assert stats['cache_hits'] == 1 and stats['cache_misses'] == 0
    assert all(np.array_equal(x, y) for x, y in zip(a['cube'], b['cube']))
    # A regenerated temporary file with identical contents is the same input.
    copied = tmp_path / 'copy.gdml'
    copied.write_bytes(path.read_bytes())
    _, stats, _ = prepare(copied, tmp_path)
    assert stats['cache_hits'] == 1


def test_cache_invalidates_content_and_resolution(tmp_path):
    path = fixture_gdml(tmp_path / 'a.gdml')
    s1, _, _ = prepare(path, tmp_path, mesh_slices=16)
    s2, stats, _ = prepare(path, tmp_path, mesh_slices=24)
    assert s1['cache_key'] != s2['cache_key'] and stats['cache_hits'] == 0
    path.write_text(path.read_text().replace('name="cube_s" x="2"', 'name="cube_s" x="3"'))
    s3, stats, meshes = prepare(path, tmp_path, mesh_slices=24)
    assert s2['cache_key'] != s3['cache_key'] and stats['cache_hits'] == 0
    assert volume(*meshes['cube']) == pytest.approx(72)


@pytest.mark.parametrize('damage', ['truncated', 'coordinates', 'indices'])
def test_corrupt_cache_is_rebuilt(tmp_path, damage):
    path = fixture_gdml(tmp_path / 'a.gdml')
    s, _, meshes = prepare(path, tmp_path)
    token = s['prototypes']['cube']['mesh_key'].encode()
    disk = tmp_path / 'cache' / s['cache_key'] / (hashlib.sha256(token).hexdigest() + '.npz')
    if damage == 'truncated':
        disk.write_bytes(b'PK broken archive')
    else:
        v, f = [a.copy() for a in meshes['cube']]
        if damage == 'coordinates':
            v[0, 0] = np.nan
        else:
            f[0, 0] = len(v)
        scene.write_mesh(disk, v, f)
    _, stats, rebuilt = prepare(path, tmp_path)
    assert stats['cache_hits'] == 0 and stats['cache_misses'] == 1
    assert volume(*rebuilt['cube']) == pytest.approx(48)


def test_disabled_cache_and_hidden_volume_do_no_unrequested_work(tmp_path):
    path = fixture_gdml(tmp_path / 'a.gdml')
    _, stats, meshes = prepare(path, tmp_path, cache=False, hidden_volumes=['cu*'])
    assert not meshes and stats['cache_misses'] == 0
    assert not (tmp_path / 'cache').exists()
    _, stats, meshes = prepare(path, tmp_path, cache=False, show_gas=True, only=['parent'])
    assert set(meshes) == {'parent'} and stats['cache_misses'] == 1
    assert not (tmp_path / 'cache').exists()


def test_settings_restore_success_and_failure(tmp_path, monkeypatch):
    old_slices = pg.config.SolidDefaults.Sphere.nslice
    old_stack = pg.config.SolidDefaults.Tubs.nslice
    old_meshing = pg.config.doMeshing
    with pytest.raises(RuntimeError):
        with scene.mesh_settings(pg, 20, meshing=False):
            assert pg.config.SolidDefaults.Sphere.nslice == 20
            assert not pg.config.doMeshing
            raise RuntimeError('bad input')
    assert pg.config.SolidDefaults.Sphere.nslice == old_slices
    assert pg.config.SolidDefaults.Tubs.nslice == old_stack
    assert pg.config.doMeshing is old_meshing
    prepare(fixture_gdml(tmp_path / 'a.gdml'), tmp_path, mesh_slices=12)
    assert pg.config.SolidDefaults.Sphere.nslice == old_slices
    assert pg.config.doMeshing is old_meshing


def test_dependency_hash_and_cycle_detection(tmp_path):
    child = fixture_gdml(tmp_path / 'child.gdml')
    parent = fixture_gdml(tmp_path / 'parent.gdml', '''
     <volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>
      <physvol name="include"><file name="%s"/></physvol>
     </volume>''' % child)
    a = scene.cache_key(scene.input_manifest(parent), pg)
    child.write_text(child.read_text().replace('name="cube_s" x="2"', 'name="cube_s" x="3"'))
    b = scene.cache_key(scene.input_manifest(parent), pg)
    assert a != b
    fixture_gdml(child, '''<volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>
      <physvol name="cycle"><file name="%s"/></physvol></volume>''' % parent)
    with pytest.raises(ValueError, match='cyclic'):
        scene.input_manifest(parent)


def test_input_change_during_load_does_not_publish_cache(tmp_path, monkeypatch):
    path = fixture_gdml(tmp_path / 'a.gdml')
    original = pg.geant4.solid.Box.mesh
    def mesh(solid):
        path.write_text(path.read_text().replace('name="cube_s" x="2"', 'name="cube_s" x="3"'))
        return original(solid)
    monkeypatch.setattr(pg.geant4.solid.Box, 'mesh', mesh)
    with pytest.raises(RuntimeError, match='changed while loading'):
        prepare(path, tmp_path)
    assert not (tmp_path / 'cache').exists()


def test_region_outer_faces_exclude_inner_wall_and_region_endcaps():
    regions = scene.normalize_regions(REGIONS)
    reg = pg.geant4.Registry()
    tube = pg.geant4.solid.Tubs('glass', 9, 10, 20, 0, 2 * np.pi, reg, nslice=24)
    v, f = scene.mesh_arrays(tube.mesh())
    outer = scene.exterior_faces(v, f, 'bare_transparent_glass', regions)
    assert 0 < len(outer) < len(f)
    assert np.linalg.norm(v[outer][:, :, :2], axis=2) == pytest.approx(np.full((len(outer), 3), 10))
    assert np.array_equal(scene.exterior_faces(v, f, 'unknown', regions), f)
    assert np.array_equal(scene.exterior_faces(v, f, 'bare_transparent_glass', None), f)
    metadata = dict(prototypes={'tube':dict(name='tube', density=2, role='bare_transparent_glass'),
                               'vacuum':dict(name='vacuum', density=0, role='vacuum')})
    assert scene.selected_prototypes(metadata, show_gas=True, display='exterior',
                                     hidden_roles=regions['hidden']) == ['tube']
    assert scene.selected_prototypes(metadata, show_gas=True, display='exterior') == ['tube', 'vacuum']


def test_exterior_display_reads_caller_regions_from_auxiliary(tmp_path):
    path = fixture_gdml(tmp_path / 'regions.gdml', '''
      <volume name="cube"><materialref ref="Steel"/><solidref ref="cube_s"/>
       <auxiliary auxtype="pmt_region" auxvalue="vacuum"/></volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>
       <physvol name="placed"><volumeref ref="cube"/></physvol></volume>''')
    with pytest.raises(ValueError, match='regions'):
        scene.prepare_scene(path, tmp_path / 'out', cache=False, display='exterior')
    with pytest.raises(ValueError, match='regions'):
        scene.normalize_regions({'styles': {}})
    s, _, meshes = prepare(path, tmp_path, cache=False, display='exterior', regions=REGIONS)
    assert s['prototypes']['cube']['role'] == 'vacuum' and not meshes
    _, _, meshes = prepare(path, tmp_path, cache=False, display='full', regions=REGIONS)
    assert set(meshes) == {'cube'}
    s, _, _ = prepare(path, tmp_path, cache=False)
    assert s['prototypes']['cube']['role'] == ''


def test_replica_and_division_compatibility(tmp_path):
    for kind in ('replica', 'division'):
        child = ('<replicavol number="2"><volumeref ref="cube"/><replicate_along_axis>'
                 '<direction x="1"/><width value="50" unit="mm"/><offset value="0" unit="mm"/>'
                 '</replicate_along_axis></replicavol>' if kind == 'replica' else
                 '<divisionvol axis="kXAxis" number="2" width="50" unit="mm"><volumeref ref="cube"/></divisionvol>')
        path = fixture_gdml(tmp_path / (kind + '.gdml'), '''
          <volume name="cube"><materialref ref="Steel"/><solidref ref="cube_s"/></volume>
          <volume name="world"><materialref ref="Air"/><solidref ref="world_s"/>%s</volume>''' % child)
        s, stats, meshes = prepare(path, tmp_path, show_gas=True)
        assert stats['eager_fallback']
        copies = [i for i in s['instances'] if '#' in i['prototype']]
        assert [np.asarray(i['matrix'])[0, 3] for i in copies] == [-25, 25]
        expected = 48 if kind == 'replica' else 50 * 100 * 100
        assert all(volume(*meshes[i['prototype']]) == pytest.approx(expected) for i in copies)


def test_partial_cache_files_are_removed_on_failure_and_when_stale(tmp_path, monkeypatch):
    cache = tmp_path / 'cache'
    cache.mkdir()
    v = np.zeros((3, 3)); f = np.array([[0, 1, 2]], dtype=np.int64)
    # An exception during the write never leaves the temporary behind.
    monkeypatch.setattr(np, 'savez_compressed', lambda *a, **k: (_ for _ in ()).throw(OSError('disk full')))
    with pytest.raises(OSError):
        scene.write_mesh(cache / 'a.npz', v, f)
    assert list(cache.iterdir()) == []
    monkeypatch.undo()
    # A worker killed mid-write leaves one; it is swept once it is clearly abandoned.
    stale, fresh = cache / (scene._PARTIAL_PREFIX + 'old.npz'), cache / (scene._PARTIAL_PREFIX + 'new.npz')
    stale.write_bytes(b'x'); fresh.write_bytes(b'x')
    old = time.time() - 2 * scene._STALE_PARTIAL_SECONDS
    os.utime(stale, (old, old))
    assert scene.sweep_partial_files(cache) == 1
    assert not stale.exists() and fresh.exists()
    fresh.unlink()
    # prepare_scene sweeps the per-key directory it is about to write into.
    path = fixture_gdml(tmp_path / 'a.gdml')
    prepare(path, tmp_path)  # prepare() uses tmp_path / 'cache'
    key_dir, = cache.iterdir()
    for cached in key_dir.glob('*.npz'):
        cached.unlink()  # force a rewrite
    stale = key_dir / (scene._PARTIAL_PREFIX + 'old.npz')
    stale.write_bytes(b'x'); os.utime(stale, (old, old))
    prepare(path, tmp_path)
    assert not stale.exists() and list(key_dir.glob('*.npz'))
    assert scene.write_mesh(cache / 'ok.npz', v, f) is None and (cache / 'ok.npz').exists()


def test_corrupt_cache_entry_does_not_skip_preview(tmp_path):
    path = fixture_gdml(tmp_path / 'a.gdml')
    s, _, _ = prepare(path, tmp_path)
    assert scene.cached_scene_available(path, cache_dir=tmp_path / 'cache')
    token = s['prototypes']['cube']['mesh_key'].encode()
    disk = tmp_path / 'cache' / s['cache_key'] / (hashlib.sha256(token).hexdigest() + '.npz')
    disk.write_bytes(b'PK broken archive')
    assert not scene.cached_scene_available(path, cache_dir=tmp_path / 'cache')


def test_input_change_during_cache_write_discards_staged_files(tmp_path, monkeypatch):
    path = fixture_gdml(tmp_path / 'a.gdml')
    original = scene.stage_mesh
    def stage(disk, v, f):
        if 'cache' in Path(disk).parts:  # only the persistent cache, not transfer files
            path.write_text(path.read_text().replace('name="cube_s" x="2"', 'name="cube_s" x="3"'))
        return original(disk, v, f)
    monkeypatch.setattr(scene, 'stage_mesh', stage)
    with pytest.warns(RuntimeWarning, match='not published'):
        _, stats, meshes = prepare(path, tmp_path)
    assert stats['cache_misses'] == 1 and volume(*meshes['cube']) == pytest.approx(48)
    assert all(not any(d.iterdir()) for d in (tmp_path / 'cache').iterdir())
