# Decay channels and propagation

A `DecayModel` describes one physical final state. Its `total_width()` is that
channel's partial width in GeV. Put every competing decay model in the vertex's
`interactions` list; their widths determine the propagation lifetime.

Use `decay_channels` to generate only a chosen subset:

```python
vertex = siren.Vertex(
    "V", interactions=[visible_decay, invisible_decay],
    decay_channels=visible_decay,
    position=position_distribution,
)
```

The full model list remains available to the position distribution and physical
weighter. The engine samples only the selected decay signatures and divides the
physical density by that restricted generation density. For two channels with
widths `Gamma_visible` and `Gamma_invisible`, forcing the visible channel gives
the branching factor

```
Gamma_visible / (Gamma_visible + Gamma_invisible)
```

once in the event weight. The mean flight distance remains
`(|p| / mass) * hbarc / (Gamma_visible + Gamma_invisible)` in metres.
Do not add another branching-fraction normalization or override each channel's
all-final-state width with the shared total. That would count the same physical
width more than once.

`decay_channels` accepts a decay model, an `InteractionSignature`, or a list of
either. A model selects all its signatures for this vertex's particle. Selection
is by signature, so every model contribution to the same signature participates.
Native models that cover several final states remain supported: each model's
all-final-state width must count precisely the physical channels it owns.

`None` leaves generation unrestricted. Empty lists, duplicates, foreign parents,
and signatures missing from the decay models raise `ConfigurationError`.
Selection order has no effect: equivalent lists compare equal and are stored
in canonical signature order. Assigning `Vertex.decay_channels` uses the same
validation as construction; a rejected assignment keeps the old selection.
Closed selected channels have zero generation density and are never sampled;
if no interaction remains, strict generation reports the failed attempt.
Every model's all-final-state width and every advertised partial width must be
finite and nonnegative, including unselected channels used for propagation.
Their total must remain finite. Zero widths remain valid closed channels.
Cross-section candidates remain eligible under a decay restriction. With
competing scattering, channel probabilities use the combined allowed interaction
rates. Selected collections also reject negative or nonfinite scattering rates
and an overflowing sum of rates; zero-rate candidates are never sampled.

The same declaration works for primary and secondary vertices, including
`Simulation`. Phase-space biasing compiles only the selected decay signatures;
its density is still included in the weight. `Fixed()` omits flight probabilities
but retains channel branching. Pure decays of a stationary parent use finite
width ratios for both channel sampling and physical channel probabilities, with
or without a selection and including when several channels compete; `Fixed()`
is the meaningful mode there. The final-state proposal density is still
included. Moving parents keep inverse flight lengths, whose common boost factor
cancels from the ratios. A parent moving so slowly that an inverse flight length
overflows is rejected with `ConfigurationError` rather than treated as stationary.
For lifetime reweighting, supply the full physical
model set at the new parameters through `physical_interactions` or a Weighter
override; leave the generation models unchanged.

For direct process construction, call
`InteractionCollection.SetDecayChannels([signature, ...])`. The lower-level
Python `Injector` constructor accepts signature lists through
`primary_decay_channels` and a particle-keyed `secondary_decay_channels` mapping.
`GetDecayChannels()` returns a copy; use the setter to change the collection.
Replacing `Injector.secondary_interactions` validates every replacement before
changing any native process, so a rejected replacement preserves all processes.

Selections survive native save/load and supported pickle paths. Older
InteractionCollection archives load with unrestricted generation. The existing
serialization restrictions on Python models, callbacks, and distributions still
apply. Readers predating archive version 1 cannot read the new collections.
Version-1 collections written before selection-order normalization still load;
their selections are canonicalized on read.
`Results.merge` requires matching selections, as part of its requirement
that pooled runs have the same configuration. A Weighter combining injectors
still requires each injector to support the events being weighted.
For truly disjoint channel selections with consistent physical normalization,
sum the independently normalized per-run channel estimates. Keep those results
separate from `Results.merge` or pooled Weighter operations. Overlapping channel
selections cannot be added this way without counting their overlap twice.

This feature controls the channels of a process. Secondary routing still allows
one process per particle type; it does not choose different channel restrictions
for repeated occurrences of that particle in a chain. Model-specific particle
aliases used for such routing need a separate migration decision. Vacuum flight
beyond detector geometry also remains a separate distribution contract.
