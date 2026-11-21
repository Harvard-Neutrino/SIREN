"""Viewer startup: independent meshes, bounded cache lifetime, requested cutters."""
from collections import Counter
from types import SimpleNamespace

import pytest

from siren import visualization


class Mesh:
    def __init__(self, position=0):
        self.position = position

    def clone(self):
        return Mesh(self.position)


class Solid:
    def __init__(self, position=0):
        self.position = position
        self.calls = 0

    def mesh(self):
        self.calls += 1
        return Mesh(self.position)


def test_cache_protects_both_first_result_and_later_clones():
    solid = Solid(2)
    registry = SimpleNamespace(solidDict={'solid': solid, 'alias': solid})
    with visualization._cached_solid_meshes(registry):
        first = solid.mesh()
        first.position += 100  # MultiUnion transforms even its first operand.
        second = solid.mesh()
        assert second.position == 2
        second.position = -30
        assert solid.mesh().position == 2
        assert solid.calls == 1
    assert 'mesh' not in vars(solid)
    solid.position = 7  # Edits after loading must re-mesh with current parameters.
    assert solid.mesh().position == 7
    assert solid.calls == 2


def test_distinct_solids_do_not_share_cached_results():
    a, b = Solid(2), Solid(7)
    with visualization._cached_solid_meshes(SimpleNamespace(solidDict={'a': a, 'b': b})):
        assert a.mesh().position == 2
        assert b.mesh().position == 7
        assert a.mesh().position == 2
    assert a.calls == b.calls == 1


def test_failure_is_not_cached_and_overrides_are_restored():
    solid = Solid()
    calls = []

    def existing_override():
        calls.append(None)
        if len(calls) == 1:
            raise ValueError('meshing failed')
        return Mesh(5)

    solid.mesh = existing_override
    with pytest.raises(RuntimeError, match='reader failed'):
        with visualization._cached_solid_meshes(SimpleNamespace(solidDict={'s': solid})):
            with pytest.raises(ValueError, match='meshing failed'):
                solid.mesh()
            assert solid.mesh().position == 5
            assert solid.mesh().position == 5
            assert len(calls) == 2
            raise RuntimeError('reader failed')
    assert solid.mesh is existing_override


def test_reader_failure_restores_mesh_methods_without_global_patch():
    solid = Solid()
    registry = SimpleNamespace(solidDict={'s': solid})

    class Reader:
        def __init__(self, path):
            self.parseStructure(SimpleNamespace(getElementsByTagName=lambda name: []))

        def getRegistry(self):
            return registry

        def parseStructure(self, xmldoc):
            solid.mesh()
            raise ValueError('invalid structure')

    original = Reader.parseStructure
    pg = SimpleNamespace(gdml=SimpleNamespace(Reader=Reader))
    with pytest.raises(ValueError, match='invalid structure'):
        visualization._read_pyg4ometry_registry('broken.gdml', pg)
    assert Reader.parseStructure is original
    assert 'mesh' not in vars(solid)


def gdml(structure):
    return '''<?xml version="1.0"?>
<gdml>
 <define/>
 <materials><material name="Air" Z="1"><D value="0.001" unit="g/cm3"/><atom value="1"/></material></materials>
 <solids>
  <box name="shared" x="2" y="4" z="6" lunit="mm"/>
  <box name="world_box" x="100" y="100" z="100" lunit="mm"/>
  <multiUnion name="pair_x">
   <multiUnionNode name="x1"><solid ref="shared"/><position x="10" y="0" z="0" unit="mm"/></multiUnionNode>
   <multiUnionNode name="x2"><solid ref="shared"/><position x="-10" y="0" z="0" unit="mm"/></multiUnionNode>
  </multiUnion>
  <multiUnion name="pair_y">
   <multiUnionNode name="y1"><solid ref="shared"/><position x="0" y="20" z="0" unit="mm"/><rotation x="0" y="0" z="90" unit="deg"/></multiUnionNode>
   <multiUnionNode name="y2"><solid ref="shared"/><position x="0" y="-20" z="0" unit="mm"/><rotation x="0" y="0" z="90" unit="deg"/></multiUnionNode>
  </multiUnion>
 </solids>
 <structure>''' + structure + '''</structure>
 <setup name="Default" version="1.0"><world ref="world"/></setup>
</gdml>'''


def test_shared_multiunion_operands_match_analytic_geometry(tmp_path, monkeypatch):
    pg = pytest.importorskip('pyg4ometry')
    np = pytest.importorskip('numpy')
    path = tmp_path / 'shared.gdml'
    path.write_text(gdml('''
      <volume name="x"><materialref ref="Air"/><solidref ref="pair_x"/></volume>
      <volume name="y"><materialref ref="Air"/><solidref ref="pair_y"/></volume>
      <volume name="cube"><materialref ref="Air"/><solidref ref="shared"/></volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_box"/>
       <physvol name="placed_x"><volumeref ref="x"/></physvol>
       <physvol name="placed_y"><volumeref ref="y"/></physvol>
       <physvol name="placed_cube"><volumeref ref="cube"/></physvol>
      </volume>'''))
    calls = Counter()
    original = pg.geant4.solid.Box.mesh

    def counted(solid):
        calls[solid.name] += 1
        return original(solid)

    monkeypatch.setattr(pg.geant4.solid.Box, 'mesh', counted)
    ordinary = pg.gdml.Reader(str(path)).getRegistry()
    assert calls['shared'] == 5
    calls.clear()
    _, cached, _, _ = visualization._load_registry(str(path), None)
    assert calls['shared'] == 1
    expected = {'x': ([-11, -2, -3], [11, 2, 3], 96),
                'y': ([-2, -21, -3], [2, 21, 3], 96),
                'cube': ([-1, -2, -3], [1, 2, 3], 48)}
    for reg in (ordinary, cached):
        for name, (lo, hi, volume) in expected.items():
            mesh = reg.logicalVolumeDict[name].mesh.localmesh
            vertices = np.array(mesh.toVerticesAndPolygons()[0])
            assert vertices.min(0) == pytest.approx(lo)
            assert vertices.max(0) == pytest.approx(hi)
            assert mesh.volume() == pytest.approx(volume)
            assert mesh.isClosed() and mesh.isOutwardOriented()
    assert all('mesh' not in vars(s) for s in cached.solidDict.values())
    # A new load and an explicit edit both construct current geometry.
    _, fresh, _, _ = visualization._load_registry(str(path), None)
    assert calls['shared'] == 2
    fresh.solidDict['shared'].pX = 8
    assert fresh.solidDict['shared'].mesh().volume() == pytest.approx(8 * 4 * 6)


def test_replicas_keep_eager_logical_volume_meshes(tmp_path):
    pytest.importorskip('pyg4ometry')
    path = tmp_path / 'replicas.gdml'
    path.write_text(gdml('''
      <volume name="cube"><materialref ref="Air"/><solidref ref="shared"/></volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_box"/>
       <replicavol number="3"><volumeref ref="cube"/><replicate_along_axis>
        <direction x="1"/><width value="2" unit="mm"/><offset value="0" unit="mm"/>
       </replicate_along_axis></replicavol>
      </volume>'''))
    _, reg, world, _ = visualization._load_registry(str(path), None)
    replica = world.daughterVolumes[0]
    assert len(replica.meshes) == 3
    assert [t[1][0] for t in replica.transforms] == [-2, 0, 2]
    assert all(m.localmesh.volume() == pytest.approx(48) for m in replica.meshes)


def test_divisions_preserve_reader_behavior_when_dimensions_mutate(tmp_path):
    pg = pytest.importorskip('pyg4ometry')
    path = tmp_path / 'divisions.gdml'
    path.write_text(gdml('''
      <volume name="cube"><materialref ref="Air"/><solidref ref="shared"/></volume>
      <volume name="divided"><materialref ref="Air"/><solidref ref="world_box"/>
       <divisionvol axis="kXAxis" number="2" width="50" unit="mm">
        <volumeref ref="cube"/>
       </divisionvol>
      </volume>
      <volume name="world"><materialref ref="Air"/><solidref ref="world_box"/>
       <physvol name="placed_divided"><volumeref ref="divided"/></physvol>
      </volume>'''))
    ordinary = pg.gdml.Reader(str(path)).getRegistry()
    _, cached, _, _ = visualization._load_registry(str(path), None)
    # Some pyg4ometry versions alter the mother's expressions during division.
    # A reused solid must reflect the same changes as the ordinary reader.
    for name in ('cube', 'divided', 'world'):
        assert cached.logicalVolumeDict[name].mesh.localmesh.volume() == pytest.approx(
            ordinary.logicalVolumeDict[name].mesh.localmesh.volume())
    division = cached.logicalVolumeDict['divided'].daughterVolumes[0]
    assert [t[1][0] for t in division.transforms] == [-25, 25]
    assert all(m.localmesh.volume() == pytest.approx(50 * 100 * 100)
               for m in division.meshes)
