"""Signature metadata must advertise the same daughters used by generation."""
from types import SimpleNamespace
import pytest
import siren
pytest.importorskip('DarkNews')
from siren.SIREN_DarkNews import PyDarkNewsCrossSection

@pytest.mark.parametrize('target_pdg,target_type',[(1000060120,siren.dataclasses.Particle.ParticleType(1000060120)),(2212,siren.dataclasses.Particle.ParticleType.HNucleus)])
@pytest.mark.parametrize('surface',['public','resource'])
def test_darknews_signature_metadata_matches_generation(target_pdg,target_type,surface):
    case=SimpleNamespace(nu_projectile=SimpleNamespace(pdgid=14),
        nu_upscattered=SimpleNamespace(pdgid=5914),nuclear_target=SimpleNamespace(pdgid=target_pdg))
    cls = PyDarkNewsCrossSection if surface == 'public' else siren._util.import_resource(
        'processes','DarkNewsTables').PyDarkNewsCrossSection
    model=cls(case)
    model.is_configured=True
    model.target_type=target_type
    advertised=model.GetPossibleSignatures()
    generated=model.GetPossibleSignaturesFromParents(siren.particles.NuMu,target_type)
    assert advertised==generated
    assert list(advertised[0].secondary_types)==[siren.particles.N4,siren.dataclasses.Particle.ParticleType(target_pdg)]
    vertex=siren.Vertex('NuMu',model)
    assert siren.particles.N4 in vertex.as_vertex_spec().secondary_types
