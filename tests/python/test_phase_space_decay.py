"""Native physical-channel adapter, including fresh archive ownership."""
import math
import pickle

import pytest
import siren


def decay(partial=2e-9, total=8e-9):
    signature = siren.dataclasses.InteractionSignature()
    signature.primary_type = siren.particles.Pi0
    signature.target_type = siren.particles.Decay
    signature.secondary_types = [siren.particles.Gamma] * 2
    return siren.injection.PhaseSpaceDecay(
        signature, [0., 0.], partial, total, siren.injection.Isotropic2BodyChannel())


def test_native_decay_physical_density_widths_sample_and_pickle():
    model = decay()
    restored = pickle.loads(pickle.dumps(model))
    assert model == restored
    assert model != decay(1e-9)
    record = siren.dataclasses.InteractionRecord()
    record.signature = model.GetPossibleSignatures()[0]
    record.primary_mass = .1349768
    record.primary_momentum = [record.primary_mass, 0., 0., 0.]
    record.secondary_masses = [0., 0.]
    record.secondary_momenta = [[0.] * 4 for _ in range(2)]
    record.secondary_helicities = [0., 0.]
    physical = siren.injection.PhysicalDecayChannel(restored)
    physical.Sample(siren.utilities.SIREN_random(44), None, record)
    assert record.secondary_momenta[0][0] == pytest.approx(record.primary_mass / 2)
    assert physical.Density(None, record) == pytest.approx(1 / (4 * math.pi))
    assert restored.TotalDecayWidth(record) == 2e-9
    assert restored.TotalDecayWidthAllFinalStates(record) == 2e-9
    assert restored.ParentDecayWidth(record) == 8e-9
    assert restored.DifferentialDecayWidth(record) == pytest.approx(2e-9 / (4 * math.pi))
    record.signature.primary_type = siren.particles.NuMu
    assert restored.TotalDecayWidth(record) == 0
    assert restored.FinalStateProbability(record) == 0


@pytest.mark.parametrize('partial,total', [(0., 1.), (-1., 2.), (2., 1.), (1., math.inf)])
def test_invalid_native_widths_rejected(partial, total):
    with pytest.raises(ValueError):
        decay(partial, total)


def make_record(model, momentum=0.2):
    r = siren.dataclasses.InteractionRecord()
    r.signature = model.GetPossibleSignatures()[0]
    r.primary_mass = .1349768
    r.primary_momentum = [math.hypot(r.primary_mass, momentum), 0, 0, momentum]
    r.secondary_masses = model.SecondaryMasses(r.signature.secondary_types)
    r.secondary_momenta = [[0.] * 4 for _ in r.secondary_masses]
    r.secondary_helicities = [0.] * len(r.secondary_masses)
    return r


@pytest.mark.parametrize('branch', [1e-6, 0.1, 1.0])
@pytest.mark.parametrize('propagated', [False, True])
def test_native_branching_and_lifetime_are_separate(branch, propagated):
    # An independent exponential flight integral, with the angular density 1/4pi.
    width = 1e-16
    model = decay(branch * width, width)
    r = make_record(model)
    collection = siren.interactions.InteractionCollection(siren.particles.Pi0, [model])
    physical = siren.injection.PhysicalProcess(siren.particles.Pi0, collection)
    injection = siren.injection.PrimaryInjectionProcess(siren.particles.Pi0, collection)
    mode = (siren.injection.VertexWeightingMode.Propagated() if propagated
            else siren.injection.VertexWeightingMode.Fixed())
    physical.weighting_mode = mode
    injection.weighting_mode = mode
    weighter = siren.injection.PrimaryProcessWeighter(physical, injection, siren.detector.DetectorModel())
    bounds = (siren.math.Vector3D(0, 0, 0), siren.math.Vector3D(0, 0, 2))
    r.interaction_vertex = [0, 0, .3]
    # SIREN's documented rounded hbar*c in GeV metres; independent length formula.
    length = .2 / r.primary_mass * 1.973e-16 / width
    assert collection.TotalDecayWidthAllFinalStates(r) == width
    assert collection.TotalDecayLengthAllFinalStates(r) == pytest.approx(length, rel=1e-9)
    expected = branch / (4*math.pi)
    if propagated:
        expected *= math.exp(-.3/length)/length
    assert weighter.PhysicalProbability(bounds, r) == pytest.approx(expected, rel=1e-8, abs=0)


def test_shared_parent_width_once_and_conflicting_declarations():
    a = decay(.7e-16, 1e-16)
    sig = a.GetPossibleSignatures()[0]
    sig.secondary_types = [siren.particles.EMinus, siren.particles.EPlus]
    b = siren.injection.PhaseSpaceDecay(sig, [0., 0.], .3e-16, 1e-16,
                                        siren.injection.Isotropic2BodyChannel())
    collection = siren.interactions.InteractionCollection(siren.particles.Pi0, [a, b])
    r = make_record(a)
    assert collection.TotalDecayWidthAllFinalStates(r) == 1e-16
    single = siren.interactions.InteractionCollection(siren.particles.Pi0, [decay(1e-16, 1e-16)])
    assert collection.TotalDecayLengthAllFinalStates(r) == single.TotalDecayLengthAllFinalStates(r)
    bad = siren.interactions.InteractionCollection(siren.particles.Pi0, [a, decay(.1e-16, 2e-16)])
    with pytest.raises(ValueError, match='Conflicting'):
        bad.TotalDecayWidthAllFinalStates(r)
    bad = siren.interactions.InteractionCollection(siren.particles.Pi0, [a, a])
    with pytest.raises(ValueError, match='exceed'):
        bad.TotalDecayWidthAllFinalStates(r)


def test_native_width_combines_with_legacy_owned_width():
    class GammaGamma(siren.DecayModel):
        parent='Pi0'; daughters=('Gamma','Gamma')
        measure=siren.Measure.SolidAngleRest()
        def total_width(self): return .9e-16
        def differential_width(self, record): return self.total_width()/(4*math.pi)
        def sample(self, record, random): self.sample_isotropic(record,random)
    legacy=GammaGamma()
    a=decay(.1e-16,1e-16)
    r=make_record(a)
    collection=siren.interactions.InteractionCollection(siren.particles.Pi0,[a,legacy])
    assert collection.ParentDecayWidth(r)==1e-16
    assert collection.TotalDecayWidthAllFinalStates(r)==1e-16
    physical=siren.injection.PhysicalProcess(siren.particles.Pi0,collection)
    injection=siren.injection.PrimaryInjectionProcess(siren.particles.Pi0,collection)
    physical.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    injection.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    w=siren.injection.PrimaryProcessWeighter(physical,injection,siren.detector.DetectorModel())
    bounds=(siren.math.Vector3D(),siren.math.Vector3D(0,0,2))
    assert w.PhysicalProbability(bounds,r)==pytest.approx(1/(4*math.pi))


def test_python_parent_width_override_dispatches_through_collection():
    class RareGammaGamma(siren.DecayModel):
        parent='Pi0'; daughters=('Gamma','Gamma')
        measure=siren.Measure.SolidAngleRest()
        def total_width(self): return .1e-16
        def ParentDecayWidth(self,record): return 1e-16
        def differential_width(self,record): return self.total_width()/(4*math.pi)
        def sample(self,record,random): self.sample_isotropic(record,random)
    model=RareGammaGamma()
    r=make_record(model)
    collection=siren.interactions.InteractionCollection(siren.particles.Pi0,[model])
    assert collection.TotalDecayWidthAllFinalStates(r)==1e-16


def test_legacy_fixed_model_does_not_require_a_parent_lifetime_query():
    class LegacyFixed(siren.DecayModel):
        parent='Pi0'; daughters=('Gamma','Gamma')
        measure=siren.Measure.SolidAngleRest()
        def total_width(self): return 1e-16
        def TotalDecayWidthAllFinalStates(self,record):
            raise AssertionError('Legacy Fixed process did not request a lifetime')
        def differential_width(self,record): return self.total_width()/(4*math.pi)
        def sample(self,record,random): self.sample_isotropic(record,random)
    model=LegacyFixed();r=make_record(model)
    collection=siren.interactions.InteractionCollection(siren.particles.Pi0,[model])
    assert collection.ParentDecayWidth(r)==0
    physical=siren.injection.PhysicalProcess(siren.particles.Pi0,collection)
    injection=siren.injection.PrimaryInjectionProcess(siren.particles.Pi0,collection)
    physical.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    injection.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    w=siren.injection.PrimaryProcessWeighter(physical,injection,siren.detector.DetectorModel())
    assert w.PhysicalProbability((siren.math.Vector3D(),siren.math.Vector3D(0,0,2)),r)==pytest.approx(1/(4*math.pi))


@pytest.mark.parametrize('propagated',[False,True])
def test_omitted_decay_width_with_material_competition(tmp_path,propagated):
    # Competing scatter and decay hazards have different units at their model
    # boundary; evaluate both as inverse metres in this independent oracle.
    materials=tmp_path/'materials.dat';densities=tmp_path/'densities.dat'
    materials.write_text('ARGON 1\n1000180400 1\n')
    densities.write_text('object sphere 0 0 0 0 0 0 100 active ARGON constant 1.4\ndetector 0 0 0\n')
    detector=siren.detector.DetectorModel();detector.LoadMaterialModel(str(materials));detector.LoadDetectorModel(str(densities))
    model=decay(.1e-16,1e-16);r=make_record(model);r.interaction_vertex=[0,0,.3]
    argon=siren.dataclasses.ParticleType.Ar40Nucleus;sigma=1e-24
    xs=siren.interactions.TrivialCrossSection(sigma,[siren.particles.Pi0],[argon])
    collection=siren.interactions.InteractionCollection(siren.particles.Pi0,[xs],[model])
    physical=siren.injection.PhysicalProcess(siren.particles.Pi0,collection)
    injection=siren.injection.PrimaryInjectionProcess(siren.particles.Pi0,collection)
    mode=siren.injection.VertexWeightingMode.Propagated() if propagated else siren.injection.VertexWeightingMode.Fixed()
    physical.weighting_mode=mode;injection.weighting_mode=mode
    w=siren.injection.PrimaryProcessWeighter(physical,injection,detector)
    nd=detector.GetParticleDensity(siren.detector.DetectorPosition(siren.math.Vector3D()),argon)
    decay_rate=1e-16*r.primary_mass/.2/1.973e-16
    total_rate=100*nd*sigma+decay_rate
    expected=.1*decay_rate/total_rate/(4*math.pi)
    if propagated: expected*=total_rate*math.exp(-.3*total_rate)
    assert w.PhysicalProbability((siren.math.Vector3D(),siren.math.Vector3D(0,0,2)),r)==pytest.approx(expected,rel=1e-10,abs=0)


def test_stationary_fixed_rare_decay_retains_branch():
    model=decay(.1e-16,1e-16);r=make_record(model,0.)
    collection=siren.interactions.InteractionCollection(siren.particles.Pi0,[model])
    physical=siren.injection.PhysicalProcess(siren.particles.Pi0,collection)
    injection=siren.injection.PrimaryInjectionProcess(siren.particles.Pi0,collection)
    physical.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    injection.weighting_mode=siren.injection.VertexWeightingMode.Fixed()
    w=siren.injection.PrimaryProcessWeighter(physical,injection,siren.detector.DetectorModel())
    assert w.PhysicalProbability((siren.math.Vector3D(),siren.math.Vector3D(0,0,2)),r)==pytest.approx(.1/(4*math.pi))


def test_darknews_base_parent_width_override_reaches_collection():
    # The pyDarkNewsDecay trampoline must forward the optional hook like pyDecay.
    Base = siren.decay_model_base(base=siren.interactions.DarkNewsDecay)

    class RareGammaGamma(Base):
        parent = 'Pi0'
        daughters = ('Gamma', 'Gamma')
        measure = siren.Measure.SolidAngleRest()
        def total_width(self): return .1e-16
        def ParentDecayWidth(self, record): return 1e-16
        def differential_width(self, record): return self.total_width()/(4*math.pi)
        def SecondaryMasses(self, types): return [0., 0.]
        def sample(self, record, random): self.sample_isotropic(record, random)

    model = RareGammaGamma()
    r = make_record(model)
    collection = siren.interactions.InteractionCollection(siren.particles.Pi0, [model])
    assert collection.ParentDecayWidth(r) == 1e-16
    assert collection.TotalDecayWidthAllFinalStates(r) == 1e-16
    physical = siren.injection.PhysicalProcess(siren.particles.Pi0, collection)
    injection = siren.injection.PrimaryInjectionProcess(siren.particles.Pi0, collection)
    physical.weighting_mode = siren.injection.VertexWeightingMode.Fixed()
    injection.weighting_mode = siren.injection.VertexWeightingMode.Fixed()
    w = siren.injection.PrimaryProcessWeighter(physical, injection, siren.detector.DetectorModel())
    bounds = (siren.math.Vector3D(), siren.math.Vector3D(0, 0, 2))
    assert w.PhysicalProbability(bounds, r) == pytest.approx(.1/(4*math.pi))


@pytest.mark.parametrize('name', ['ParentDecayWidht', 'ParentDecayWidths', 'ParentWidth'])
def test_parent_width_misspelling_is_rejected(name):
    body = {'parent': 'Pi0', 'daughters': ('Gamma', 'Gamma'),
            'measure': siren.Measure.SolidAngleRest(),
            'total_width': lambda self: 1e-16,
            name: lambda self, record: 1e-16}
    with pytest.raises(siren.errors.ConfigurationError, match='ParentDecayWidth'):
        type('Misspelled', (siren.DecayModel,), body)
