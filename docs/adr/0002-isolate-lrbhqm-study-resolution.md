# Isolate LRBHQM study resolution

LRBHQM study generation uses a dedicated `WormWeaver.LRBHQM.generate_study` entry point that resolves model-specific parameters before delegating to WormWeaver's internal generic generation pipeline. It accepts the ordinary study schema, derives `t/U₁₁` from each expanded configuration's physical `model.hopping.J` and onsite interaction matrix, and resolves `model.onsite.mu`; ordinary `generate_study` retains the configured chemical potential. This keeps generic generation unaware of LRBHQM physics while avoiding duplicated sweep, manifest, provenance, and publication logic.

Unit-density lookup uses crossings labeled with sampling mode `none` by default, even when the simulation schedule enables G2 entry. Callers may enable mode matching to derive the lookup label from the actual schedule; the generated manifest records which policy and lookup mode were used.

Each LRBHQM study names one bundled curated JLD2 dataset through `unit_density_data`. The resolver loads only that file and verifies its model parameters against every expanded configuration; studies spanning different curated model families remain separate.
