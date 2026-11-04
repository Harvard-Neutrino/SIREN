"""Bounded translation of deprecated Simulation convenience arguments."""
import copy
import warnings

from . import distributions as d, particles
from .vertex import Vertex
from ._util import load_processes

SLOTS = frozenset(('interactions', 'energy', 'direction', 'position',
    'injection_energy', 'physical_energy', 'injection_direction',
    'physical_direction', 'injection_distributions', 'physical_distributions',
    'flux', 'mass', 'targets', 'process', 'isoscalar',
    'secondary_interactions', 'secondary_position', 'weighting'))


def translate(primary, secondaries, detector, options, slots):
    warnings.warn('Simulation convenience slots are deprecated; declare primary=Vertex(...) '
                  'and secondaries=[Vertex(...)].', DeprecationWarning, stacklevel=3)
    slots = {key: value for key, value in slots.items() if value is not None}
    metadata = ()
    if isinstance(primary, Vertex):
        vertex = copy.copy(primary)
        if 'interactions' in slots or options:
            raise TypeError('A primary Vertex already declares its interactions')
        injection = list(vertex.distributions)
        physical = list(vertex.physical)
    else:
        models = slots.pop('interactions', None)
        if models is None:
            raise TypeError('Simulation requires interactions for a legacy primary')
        if isinstance(models, str):
            kwargs = dict(options)
            kwargs.setdefault('primary_types', [particles.resolve(primary)])
            for key, target in (('targets', 'target_types'), ('process', 'process_types'),
                                ('isoscalar', 'isoscalar')):
                if key in slots:
                    if target in kwargs:
                        raise TypeError('Duplicate process option: ' + target)
                    value = slots.pop(key)
                    if key == 'targets':
                        if value == 'detector':
                            from ._util import get_detector_model_targets
                            value = get_detector_model_targets(detector)
                        elif isinstance(value, str):
                            value = [particles.resolve(value)]
                        else:
                            value = [particles.resolve(v) if isinstance(v, str) else v for v in value]
                    elif key == 'process':
                        value = [value] if isinstance(value, str) else value
                    kwargs[target] = value
            bundle = load_processes(models, **kwargs)
            models = bundle.primary[particles.resolve(primary)]
            metadata = bundle.metadata
            if bundle.secondary and 'secondary_interactions' not in slots:
                slots['secondary_interactions'] = bundle.secondary
        elif options:
            raise TypeError('process_options requires a named interaction loader')
        vertex = Vertex(primary, models)
        injection, physical = [], []
    for option in ('targets', 'process', 'isoscalar'):
        if option in slots:
            raise TypeError(option + ' requires a named interaction loader')
    energy = slots.pop('energy', None)
    inj_energy = slots.pop('injection_energy', None)
    phys_energy = slots.pop('physical_energy', None)
    flux = slots.pop('flux', None)
    if flux is not None and phys_energy is not None:
        raise TypeError('Cannot specify both flux and physical_energy')
    phys_energy = flux if flux is not None else phys_energy
    if energy is not None:
        if inj_energy is not None:
            raise TypeError('Cannot specify both energy and injection_energy')
        inj_energy = energy
        phys_energy = energy if phys_energy is None else phys_energy
    if phys_energy is None and inj_energy is not None:
        phys_energy = inj_energy
    direction = slots.pop('direction', None)
    inj_direction = slots.pop('injection_direction', None)
    phys_direction = slots.pop('physical_direction', None)
    if direction is not None:
        if inj_direction is not None or phys_direction is not None:
            raise TypeError('Cannot mix direction with split direction slots')
        inj_direction = phys_direction = direction
    if phys_direction is None and inj_direction is not None:
        phys_direction = inj_direction
    injection += [x for x in (inj_energy, inj_direction, slots.pop('position', None))
                  if x is not None]
    injection += list(slots.pop('injection_distributions', ()) or ())
    for replacement, base in ((phys_energy, d.PrimaryEnergyDistribution),
                              (phys_direction, d.PrimaryDirectionDistribution)):
        if replacement is not None:
            physical = [item for item in physical if not isinstance(item, base)]
            physical.append(replacement)
    physical += list(slots.pop('physical_distributions', ()) or ())
    mass = slots.pop('mass', None)
    has_mass = any(d.DistributionVariable.PrimaryMass in x.ProvidedVariables()
                   for x in injection)
    if mass is not None and has_mass:
        raise TypeError('mass and an explicit mass distribution conflict')
    if not has_mass:
        injection.insert(0, d.PrimaryMass(0.0 if mass is None else mass))
    vertex.distributions = injection
    vertex.physical = physical
    if 'weighting' in slots:
        vertex.weighting = slots.pop('weighting')
    secondary_models = slots.pop('secondary_interactions', None)
    positions = slots.pop('secondary_position', None)
    if secondary_models is not None:
        if secondaries:
            raise TypeError('Cannot mix secondaries with secondary_interactions')
        if positions is None:
            raise TypeError('secondary_interactions requires secondary_position')
        secondaries = [Vertex(particles.resolve(kind) if isinstance(kind, str) else kind,
                              models, position=positions[kind] if isinstance(positions, dict)
                              else positions)
                       for kind, models in secondary_models.items()]
    elif positions is not None:
        raise TypeError('secondary_position requires secondary_interactions')
    if slots:
        raise TypeError('Unused Simulation arguments: ' + ', '.join(sorted(slots)))
    return vertex, tuple(secondaries), metadata
