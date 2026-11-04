# Configuration compatibility

U04 introduces the final authoring interface on the reviewed native contracts.
These observable differences are intentional and covered by the Python tests.
The shared installed runtime and downstream BSM adapters are separate migrations.

| Contract | Behavior and reason |
|---|---|
| Repeated generation | Attempts and failures accumulate until `reset()`. Weights use the combined attempts, including failures. Recompute weights after generating all batches. This repairs the initial draft's denominator reset. |
| Requested count | Retry budgets do not replace the configured event count. Saving and loading retain that count and the RNG stream. The constructor's `max_attempts` setting limits each generation call and survives Python pickle; native archives restore the default retry policy. |
| Efficiency floor | A run with no successes stops after the 50-attempt warmup when below the requested floor. Budgets are bounded by the native counter's remaining capacity. |
| Mixed interfaces | Vertex declarations cannot be combined with legacy process keywords. Inherited Weighter construction rejects legacy keywords and unknown overrides instead of dropping them. |
| Physical target | The first pooled injector supplies the shared physical target unless explicitly overridden. Proposal models may differ; the complete weighting modes must agree for each process type. |
| Distribution updates | Generation distributions cannot be replaced through a built Python Injector. Create a new injector so sampling and weighting bounds remain consistent. |
| Required variables | Transport-dependent position samplers require energy and mass before sampling. Incorrectly ordered configurations now fail during construction. |
| Secondary control | Configured secondaries require explicit expansion or a legacy stopping condition; silently pruning an intended chain is an error. The legacy DipolePortal example now supplies its existing stop function. |
| Expansion | Empty lists terminate. Repeated primary/secondary particle types keep distinct rules and predicates. Invalid rules, absent daughters and unreachable secondary declarations fail early. |
| Native expansion archives | SecondaryExpansion version 1 adds a separate primary policy. Real version-0 expansion archives load with their original policy. Injector remains version 2. |
| Callback state | Native archives refuse callbacks. Python Injector pickle supports picklable stopping functions and predicates; lambdas and other non-picklable state fail explicitly. |
| Weighter pickle | Native generation and resolved physical processes are stored together. Python subclasses' additional attributes remain Python pickle state. Pickling builds configured weighters first: incomplete configurations and live injector callbacks fail explicitly even before first use. Empty weighters remain picklable. Python physics trampolines and `event_factor` still cannot be serialized as native physics. |
| Archive loading | Class loading constructs an injector; legacy instance loading replaces that instance. Duplicate secondary process types are rejected, including ambiguous old archives, instead of silently overwriting one entry. |
| Stopping inspection | Native `GetStoppingCondition()` supplies a stop-all callable when no callback exists. Use `HasStoppingCondition()` or the facade's `stopping_condition` to inspect whether one was supplied. |
| Channel sums | Addition preserves written coefficients across more than two terms; normalization occurs for the resulting mixture. This repairs proposal shares in the initial draft. |
| Probe diagnostics | Density-probe diagnostics and skipped templates are warnings, not discarded results or closure certificates. |
| DarkNews | `always_interpolate=False` retains exact differential evaluation. The authoring port replaces process exits with exceptions and validates numeric results; it does not include the separate DarkNews 0.4.9 API update. |

The retained main/U02 archive and scattering witnesses cover the unchanged native
paths. These checks do not establish full DarkNews cascade or downstream model
acceptance.
